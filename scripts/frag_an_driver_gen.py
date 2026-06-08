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
    "Smear Analysis": 1.0,  # Range 0.5 – 5 ng/uL, ideal 1 ng/uL
    "Single Peak": 0.05,  # Range 5 – 500 pg/uL (0.005 – 0.5 ng/uL), ideal ~50 pg/uL
}
DEFAULT_ANALYSIS_TYPE = "Smear Analysis"

MIN_DILUTION_FACTOR = 5.0  # Always at least 1:5 dilution
MIN_SAMPLE_VOL = 2.0  # Minimum sample volume (uL)
TARGET_TOTAL_VOL = 10.0  # Starting target for total volume (uL)


SUPPORTED_MASS_UNITS = {"ng/ul", "ng/uL", "ng/mL", "pg/uL", "pg/ul"}


def _to_ng_ul(conc, units):
    """Convert concentration value to ng/uL."""
    if units in ("ng/ul", "ng/uL"):
        return conc
    elif units == "ng/mL":
        return conc / 1000.0
    elif units in ("pg/uL", "pg/ul"):
        return conc / 1000.0
    return conc  # assume ng/uL if unit is unrecognised


def _get_conc_units(udf_dict):
    """Return concentration units from a UDF dict, defaulting to ng/uL."""
    return udf_dict["Conc. Units"] if "Conc. Units" in udf_dict else "ng/uL"


def _get_latest_ngul_measurement(inp_art):
    """Find the latest ng/uL concentration from a Quant-iT or Plate Reader
    measurement ResultFile for the same sample.

    Searches all ResultFile artifacts associated with the sample, filters those
    whose name contains 'Measurement' and whose 'Concentration' UDF is in a
    supported mass unit, then returns the value from the most recent process
    (sorted by the numeric part of the process ID).

    Returns (conc_ng_ul, None) on success, or (None, warning_str) on failure.
    """
    sample_id = inp_art.samples[0].id
    sample_name = inp_art.samples[0].name

    result_files = inp_art.lims.get_artifacts(samplelimsid=sample_id, type="ResultFile")

    candidates = []
    for rf in result_files:
        if "Measurement" not in rf.name:
            continue
        conc = rf.udf.get("Concentration")
        units = _get_conc_units(rf.udf)
        if conc is not None and units in SUPPORTED_MASS_UNITS:
            candidates.append(rf)

    if not candidates:
        return None, (
            f"No ng/uL Quant-iT or Plate Reader measurement found for '{sample_name}'."
        )

    # Most recent process has the highest numeric ID
    candidates.sort(
        key=lambda rf: int(rf.parent_process.id.split("-")[1]), reverse=True
    )
    latest = candidates[0]
    conc = float(latest.udf["Concentration"])
    units = _get_conc_units(latest.udf)
    return _to_ng_ul(conc, units), None


def get_concentration(inp_art):
    """Return (conc_ng_ul, warning_str) for an input artifact.

    Concentration is sourced in the following priority order:
      1. 'Final Concentration' UDF, when the parent step is 'Diluting Samples'.
      2. The most recent Quant-iT or Plate Reader measurement ResultFile (ng/uL).
      3. 'Concentration' UDF on the current artifact, if in a supported mass unit.

    Returns (None, message) when no usable value is found.
    """
    try:
        try:
            if inp_art.parent_process.type.name == "Diluting Samples":
                conc = float(inp_art.udf["Final Concentration"])
                units = _get_conc_units(inp_art.udf)
                return _to_ng_ul(conc, units), None
        except AttributeError:
            pass

        # Priority 1: latest ng/uL measurement ResultFile (Quant-iT / Plate Reader)
        conc_ng_ul, _ = _get_latest_ngul_measurement(inp_art)
        if conc_ng_ul is not None:
            return conc_ng_ul, None

        # Priority 2: concentration on the current artifact (mass unit only)
        if "Concentration" in inp_art.udf:
            conc = float(inp_art.udf["Concentration"])
            units = _get_conc_units(inp_art.udf)
            if units in SUPPORTED_MASS_UNITS:
                return _to_ng_ul(conc, units), None

        return None, (
            f"No ng/uL concentration found for sample '{inp_art.samples[0].name}'."
        )

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

    Returns (sample_vol, eb_vol, dilution_factor, used_default_df), where
    used_default_df is True only when neither concentration nor blanket dilution
    factor is available and the minimum default dilution factor (1:5) is applied.
    """
    used_default_df = False
    if force_blanket and blanket_df is not None:
        df = max(float(blanket_df), MIN_DILUTION_FACTOR)
    elif not force_blanket and conc_ng_ul is not None and conc_ng_ul > 0:
        df = max(conc_ng_ul / target_conc, MIN_DILUTION_FACTOR)
    elif blanket_df is not None:
        df = max(float(blanket_df), MIN_DILUTION_FACTOR)
    else:
        # No concentration and no blanket factor: use default 1:5 dilution
        df = MIN_DILUTION_FACTOR
        used_default_df = True

    sample_vol = TARGET_TOTAL_VOL / df
    if sample_vol < MIN_SAMPLE_VOL:
        # Increase total volume to honour minimum pipetting volume
        sample_vol = MIN_SAMPLE_VOL

    total_vol = sample_vol * df
    eb_vol = total_vol - sample_vol
    return round(sample_vol, 2), round(eb_vol, 2), round(df, 4), used_default_df


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
    samples_with_default_df = []
    concentration_warning_count = 0
    total_samples = 0
    samples_with_conc = 0

    # First pass: collect all samples and their concentration status
    all_samples = []
    for inp, out in currentStep.input_output_maps:
        if out["output-generation-type"] != "PerInput":
            continue
        if not out["uri"].location[1]:  # skip outputs without a well location
            continue

        inp_art = inp["uri"]
        out_art = out["uri"]

        sample_name = inp_art.samples[0].name
        source_well = inp_art.location[1].replace(":", "")
        dest_well = out_art.location[1].replace(":", "")
        total_samples += 1

        conc_ng_ul, warn = get_concentration(inp_art)
        if warn:
            concentration_warning_count += 1
        if conc_ng_ul is not None and conc_ng_ul > 0:
            samples_with_conc += 1

        all_samples.append(
            (inp_art.samples[0], sample_name, source_well, dest_well, conc_ng_ul)
        )

    # Check if no samples have concentration and no blanket factor is set
    if total_samples > 0 and samples_with_conc == 0 and blanket_df is None:
        sys.exit(
            "ERROR: No samples in this project have usable concentration and "
            "'Blanket Dilution Factor' is not set. "
            "Please set 'Blanket Dilution Factor' and check "
            "'Force Blanket Dilution' to 'Yes' in the LIMS step."
        )

    # Second pass: calculate volumes
    for sample, sample_name, source_well, dest_well, conc_ng_ul in all_samples:
        sample_vol, eb_vol, dilution_factor, used_default_df = calculate_volumes(
            conc_ng_ul, target_conc, blanket_df, force_blanket
        )

        if used_default_df:
            samples_with_default_df.append(sample_name)

        sample.udf["FA Dilution Fold"] = dilution_factor
        sample.put()

        rows.append((sample_name, source_well, sample_vol, eb_vol, dest_well))

    if concentration_warning_count:
        warnings.append(
            f"WARNING: {concentration_warning_count} sample(s) had missing or "
            "unusable concentration information."
        )

    if samples_with_default_df:
        warnings.append(
            f"WARNING: {len(samples_with_default_df)} sample(s) used standard "
            f"1:5 dilution due to missing concentration."
        )

    # --- Add ladder based on instrument ---
    if instrument == "fatboy":
        # Ladder in well 12 of the last column used
        if rows:
            columns_used = set(row[4][0] for row in rows)
            last_column = sorted(columns_used)[-1]
            ladder_well = f"{last_column}12"
            rows.append(("ladder", "", "N/A", "N/A", ladder_well))
    elif instrument == "fergie":
        # Ladder always in H12
        rows.append(("ladder", "", "N/A", "N/A", "H12"))

    # --- Write CSV ---
    filename = f"frag_an_driver_{instrument}.csv"
    with open(filename, "w") as f:
        f.write(
            "Sample Name,Source Well,Sample Volume (uL),EB Volume (uL),Destination Well\n"
        )
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
