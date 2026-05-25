#!/usr/bin/env python

import sys
from argparse import ArgumentParser

from genologics.config import BASEURI, PASSWORD, USERNAME
from genologics.entities import Process
from genologics.lims import Lims

DESC = """EPP used to create a dilution CSV file to drive the fragment analyzer.

The CSV filename includes the instrument name (fergie or fatboy).
Per-sample dilution volumes are calculated from prior quantification (Qubit /
Quant-IT) stored on the input artifacts.  If no concentration is available, a
blanket dilution factor set on the step UDF is used instead.

Required step UDFs:
  - Instrument Used               : must contain 'Fatboy' or 'Fergie'
Optional step UDFs:
  - Analysis Type                 : Smear Analysis (default) | Single Peak
  - Blanket Dilution Factor       : numeric, used when no concentration is found
  - Force Blanket Dilution        : checkbox; when ticked, blanket factor is used
                                    for ALL samples, ignoring measured concentrations
"""

# Target concentrations (ng/uL) for each analysis mode
ANALYSIS_TYPES = {
    "Smear Analysis": 1.0,   # Range 0.5 – 5 ng/uL, ideal 1 ng/uL
    "Single Peak": 0.05,     # Range 5 – 500 pg/uL (0.005 – 0.5 ng/uL), ideal ~50 pg/uL
}
DEFAULT_ANALYSIS_TYPE = "Smear Analysis"

MIN_DILUTION_FACTOR = 5.0   # Always at least 1:5 dilution
MIN_SAMPLE_VOL = 2.0        # Minimum sample volume (uL)
TARGET_TOTAL_VOL = 10.0     # Starting target for total volume (uL)


def _to_ng_ul(conc, units):
    """Convert concentration value to ng/uL."""
    if units in ("ng/ul", "ng/uL"):
        return conc
    elif units == "ng/mL":
        return conc / 1000.0
    elif units in ("pg/uL", "pg/ul"):
        return conc / 1000.0
    return conc  # assume ng/uL if unit is unrecognised


def get_concentration(inp_art):
    """Return (conc_ng_ul, warning_str) for an input artifact.

    Follows the same pattern as bravo_csv.py: if the input came from a
    'Diluting Samples' step, read 'Final Concentration' instead of
    'Concentration'.  Returns (None, message) when no value is found.
    """
    try:
        try:
            if inp_art.parent_process.type.name == "Diluting Samples":
                conc = float(inp_art.udf["Final Concentration"])
                units = inp_art.udf["Conc. Units"] if "Conc. Units" in inp_art.udf else "ng/uL"
                return _to_ng_ul(conc, units), None
        except AttributeError:
            pass

        if "Concentration" not in inp_art.udf:
            return None, (
                f"No 'Concentration' UDF found for sample "
                f"'{inp_art.samples[0].name}' — no value to calculate dilution from."
            )

        conc = float(inp_art.udf["Concentration"])
        units = inp_art.udf["Conc. Units"] if "Conc. Units" in inp_art.udf else "ng/uL"
        return _to_ng_ul(conc, units), None

    except Exception as e:
        return None, (
            f"Could not read concentration for '{inp_art.samples[0].name}': {e}"
        )


def calculate_volumes(conc_ng_ul, target_conc, blanket_df=None, force_blanket=False):
    """Calculate sample and EB volumes for one dilution.

    When force_blanket is True the blanket_df is used for every sample,
    ignoring any measured concentration.  When False (default), the measured
    concentration is used and blanket_df is only a fallback for samples that
    have no concentration value.

    Returns (sample_vol, eb_vol) rounded to 2 d.p., or (None, None) when
    neither a measured concentration nor a blanket dilution factor is available.
    """
    if force_blanket and blanket_df is not None:
        df = max(float(blanket_df), MIN_DILUTION_FACTOR)
    elif not force_blanket and conc_ng_ul is not None and conc_ng_ul > 0:
        df = max(conc_ng_ul / target_conc, MIN_DILUTION_FACTOR)
    elif blanket_df is not None:
        df = max(float(blanket_df), MIN_DILUTION_FACTOR)
    else:
        return None, None

    sample_vol = TARGET_TOTAL_VOL / df
    if sample_vol < MIN_SAMPLE_VOL:
        # Increase total volume to honour minimum pipetting volume
        sample_vol = MIN_SAMPLE_VOL

    total_vol = sample_vol * df
    eb_vol = total_vol - sample_vol
    return round(sample_vol, 2), round(eb_vol, 2)


def main(lims, args):
    currentStep = Process(lims, id=args.pid)

    # --- Required: instrument name ---
    try:
        instrument_udf = currentStep.udf["Instrument Used"]
    except KeyError:
        sys.exit("ERROR: 'Instrument Used' UDF is not set on the step.")

    if "Fatboy" in instrument_udf:
        instrument = "fatboy"
    elif "Fergie" in instrument_udf:
        instrument = "fergie"
    else:
        sys.exit(
            f"ERROR: Could not determine instrument from 'Instrument Used' value: "
            f"'{instrument_udf}'. Expected a value containing 'Fatboy' or 'Fergie'."
        )

    # --- Optional: analysis type (dropdown) ---
    if "Analysis Type" in currentStep.udf:
        analysis_type = currentStep.udf["Analysis Type"]
    else:
        analysis_type = DEFAULT_ANALYSIS_TYPE

    if analysis_type not in ANALYSIS_TYPES:
        sys.exit(
            f"ERROR: Unknown Analysis Type '{analysis_type}'. "
            f"Valid values: {list(ANALYSIS_TYPES.keys())}"
        )
    target_conc = ANALYSIS_TYPES[analysis_type]

    # --- Optional: blanket dilution factor and force flag ---
    blanket_df = None
    if "Blanket Dilution Factor" in currentStep.udf:
        blanket_df = currentStep.udf["Blanket Dilution Factor"]

    force_blanket = False
    if "Force Blanket Dilution" in currentStep.udf:
        force_blanket = currentStep.udf["Force Blanket Dilution"] == "Yes"

    if force_blanket and blanket_df is None:
        sys.exit(
            "ERROR: 'Force Blanket Dilution' is checked but "
            "'Blanket Dilution Factor' is not set."
        )

    # --- Locate the Driver File output artifact ---
    driver_file_out = None
    for output in currentStep.all_outputs():
        if output.name == "Driver File":
            driver_file_out = output
            break

    if driver_file_out is None:
        sys.exit("ERROR: No 'Driver File' output artifact found in this step.")

    # --- Build dilution rows from input-output maps ---
    rows = []
    warnings = []

    for inp, out in currentStep.input_output_maps:
        if out["output-generation-type"] != "PerInput":
            continue
        if out["uri"].type != "Analyte":
            continue

        inp_art = inp["uri"]
        out_art = out["uri"]

        sample_name = inp_art.samples[0].name
        source_well = inp_art.location[1].replace(":", "")
        dest_well = out_art.location[1].replace(":", "")

        conc_ng_ul, warn = get_concentration(inp_art)
        if warn:
            warnings.append(warn)

        sample_vol, eb_vol = calculate_volumes(
            conc_ng_ul, target_conc, blanket_df, force_blanket
        )

        if sample_vol is None:
            warnings.append(
                f"No concentration or blanket dilution factor available for "
                f"'{sample_name}' — volumes left as N/A."
            )
            rows.append((sample_name, source_well, "N/A", "N/A", dest_well))
        else:
            rows.append((sample_name, source_well, sample_vol, eb_vol, dest_well))

    # --- Write CSV ---
    filename = f"frag_an_driver_{instrument}.csv"
    with open(filename, "w") as f:
        f.write("Sample Name,Source Well,Sample Volume (uL),EB Volume (uL),Destination Well\n")
        for row in rows:
            f.write(f"{row[0]},{row[1]},{row[2]},{row[3]},{row[4]}\n")

    lims.upload_new_file(driver_file_out, filename)

    if warnings:
        sys.stderr.write("\n".join(warnings) + "\n")
        sys.exit(1)


if __name__ == "__main__":
    parser = ArgumentParser(description=DESC)
    parser.add_argument("--pid", help="Lims id for current Process")
    args = parser.parse_args()

    lims = Lims(BASEURI, USERNAME, PASSWORD)
    lims.check_version()
    main(lims, args)
