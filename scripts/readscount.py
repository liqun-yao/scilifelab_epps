#!/usr/bin/env python
DESC = """EPP script to aggregate the number of reads from different demultiplexing runs,
based on the flag 'include reads' located at the same level as '# reads'

Denis Moreno, Science for Life Laboratory, Stockholm, Sweden
"""
import logging
import os
from argparse import ArgumentParser
from datetime import datetime as dt

from genologics.config import BASEURI, PASSWORD, USERNAME
from genologics.entities import Process, Sample
from genologics.lims import Lims

from scilifelab_epps.epp import attach_file
from scilifelab_epps.wrapper import epp_decorator

TIMESTAMP: str = dt.now().strftime("%y%m%d_%H%M%S")

# Step names to look for
DEMULTIPLEX_STEPS = [
    "Bcl Conversion & Demultiplexing (Illumina SBS) 4.0",
    "ONT Finish Sequencing v3",
    "ONT Finish Sequencing v2.1",
    "Bcl Conversion & Demultiplexing (AVITI) v1.0",
]
SEQUENCING_STEPS = [
    "Illumina Sequencing (Illumina SBS) 4.0",
    "MiSeq Run (MiSeq) 4.0",
    "Illumina Sequencing (HiSeq X) 1.0",
    "AUTOMATED - NovaSeq Run (NovaSeq 6000 v2.0)",
    "Illumina Sequencing (NextSeq) v1.0",
    "NovaSeqXPlus Run v1.0",
    "ONT Start Sequencing v3.0",
    "AVITI Run v1.0",
    "Illumina Sequencing (MiSeq i100) v1.0",
]


def sample_has_status(
    sample: Sample,
    status: str,
    udfs_to_check: list[str] = ["Status (manual)"],
) -> bool:
    """Check the given status of a given sample by looking through a list of status UDFs.
    Return true if the status is found in any of the UDFs, false otherwise.
    """
    assert status in ["In Progress", "Finished", "Aborted"]
    for sample_udf_name in udfs_to_check:
        if sample.udf.get(sample_udf_name) == status:
            return True
    return False


@epp_decorator(script_path=__file__, timestamp=TIMESTAMP)
def main(args):
    """This should be run at project summary level"""
    process = Process(lims, id=args.pid)

    sample_counter = 0
    error_counter = 0

    summary = {}  # { sample_name : { flowcell : { lane1, lane2, ... } } }  # dict -> dict -> set
    try:
        summary_artifact = [
            art
            for art in process.all_outputs()
            if art.type == "ResultFile" and art.name == "AggregationSummary"
        ][0]
    except IndexError:
        # Skip the summary artifact, made as a patch for open steps
        summary_artifact = None
        logging.info("Could not find the summary file slot artifact, skipping.")

    # Iterate across output analytes
    arts_out = sorted(
        [art for art in process.all_outputs() if art.type == "Analyte"],
        key=lambda art: art.name,
    )
    for art_out in arts_out:
        assert len(art_out.samples) == 1, (
            f"Found {len(art_out.samples)} samples for the output analyte {art_out.id}, that should not happen"
        )

        sample = art_out.samples[0]
        sample_counter += 1

        # Traceback and calculate the total number of reads
        total_reads = sum_reads(sample, summary)

        # Set total reads for sample and artifact UDFs
        sample.udf["Total Reads (M)"] = total_reads
        art_out.udf["Set Total Reads"] = total_reads
        logging.info(f"Total reads is {total_reads} (M) for sample '{sample.name}'")

        # Set min reads sample UDF from project UDF
        min_reads = sample.project.udf.get("Reads Min", 0) / 1e6
        logging.info(f"Updating '{sample.name}' UDF 'Reads Min' to {min_reads}")
        sample.udf["Reads Min"] = min_reads

        # Set sample UDFs for status and sequencing QC based on min reads and total reads
        if total_reads <= min_reads:
            logging.info(
                "Total reads is below minimum, setting status to 'In Progress' and 'Passed Sequencing QC' to 'False'"
            )
            sample.udf["Status (auto)"] = "In Progress"
            sample.udf["Passed Sequencing QC"] = "False"
        elif total_reads > min_reads:
            logging.info(
                "Total reads is above minimum, setting status to 'Finished' and 'Passed Sequencing QC' to 'True'"
            )
            sample.udf["Passed Sequencing QC"] = "True"
            sample.udf["Status (auto)"] = "Finished"

        # Commit changes to sample and sample artifact
        sample.put()
        art_out.put()

    # Write the csv file, separated by pipes, no cell delimiter
    with open("AggregationSummary.csv", "w") as f:
        f.write("sep=,\n")
        f.write(
            "sample name,number of flowcells,number of lanes,flowcell1:lane1|lane2;flowcell2:lane1|lane2|lane3 ...\n"
        )
        for sample in summary:
            view = []
            n_flowcells = len(summary[sample])
            n_lanes = 0
            for fc in summary[sample]:
                view.append("{}:{}".format(fc, "|".join(summary[sample][fc])))
                n_lanes += len(summary[sample][fc])
            f.write(
                "{},{},{},{}\n".format(sample, n_flowcells, n_lanes, ";".join(view))
            )

    # Upload the summary file
    if summary_artifact is not None:
        try:
            attach_file(
                os.path.join(os.getcwd(), "AggregationSummary.csv"), summary_artifact
            )
            logging.info(
                f"Updated {sample_counter} samples with {error_counter} errors"
            )
        except AttributeError:
            # Happens if the summary artifact does not exist, if the step has been started before the configuration changes
            logging.info("Could not upload the summary file")


def sum_reads(sample, summary):
    """For a given submitted sample and a summary object,
    calculate the total number of reads and append to the summary.
    """

    logging.info(f"Aggregating reads of sample '{sample.name}'...")

    # Append to summary
    if sample.name not in summary:
        summary[sample.name] = {}

    # Look for artifacts matching the sample name in the demultiplexing processes.
    # ONT demux artifacts are named "{sample} (FASTQ reads) {barcode}" so the exact-match
    # name filter cannot be used; query without name and filter by substring instead.
    ont_steps = [s for s in DEMULTIPLEX_STEPS if s.startswith("ONT")]
    non_ont_steps = [s for s in DEMULTIPLEX_STEPS if not s.startswith("ONT")]

    demux_arts = []
    if non_ont_steps:
        demux_arts += lims.get_artifacts(
            sample_name=sample.name,
            process_type=non_ont_steps,
            name=f"{sample.name} (FASTQ reads)",
        )
    if ont_steps:
        ont_arts = lims.get_artifacts(
            sample_name=sample.name,
            process_type=ont_steps,
        )
        demux_arts += [art for art in ont_arts if "(FASTQ reads)" in art.name]

    if not demux_arts:
        if sample_has_status(sample, "Aborted"):
            logging.info(
                f"Sample {sample.name} has no demux artifacts, but is set to 'Aborted'."
            )
        else:
            logging.warning(
                f"Could not find any demultiplexing artifacts for sample {sample.name}."
            )

    # Iterate across found demux artifacts to aggregate reads and collect flowcell information
    # Check for ongoing demux steps and deduplicate by flowcell+lane, keeping only the latest demux step
    ongoing_demux_arts = []
    demux_arts_deduplicated = {}  # key: (flowcell_id, lane), value: (demux_artifact, demux_art_parent)
    for demux_art in demux_arts:
        if demux_art.parent_process.date_run is None:
            ongoing_demux_arts.append(demux_art)
        # Get flowcell and lane identifiers early for deduplication
        demux_art_parents = [
            parent
            for parent in get_parent_inputs(demux_art)
            if sample.id in [parent_sample.id for parent_sample in parent.samples]
        ]
        if not demux_art_parents or len(demux_art_parents) != 1:
            continue
        demux_art_parent = demux_art_parents[0]

        # Determine flowcell key
        if "ONT flow cell ID" in demux_art_parent.udf:
            flowcell_key = (demux_art_parent.udf["ONT flow cell ID"], "ONT")
        else:
            flowcell_key = (
                demux_art_parent.location[0].name,
                demux_art_parent.location[1].split(":")[0],
            )

        # Keep only the latest by date_run, with artifact ID as tiebreaker
        if flowcell_key not in demux_arts_deduplicated:
            demux_arts_deduplicated[flowcell_key] = (demux_art, demux_art_parent)
        else:
            existing_art, _ = demux_arts_deduplicated[flowcell_key]
            # Compare date_run: newer (later date) wins
            existing_date = existing_art.parent_process.date_run
            new_date = demux_art.parent_process.date_run

            # Extract numeric ID for tiebreaker (e.g., "92-7405142" -> 7405142)
            existing_id_num = int(existing_art.id.split("-")[1])
            new_id_num = int(demux_art.id.split("-")[1])

            should_replace = False
            if new_date is not None and (
                existing_date is None or new_date > existing_date
            ):
                should_replace = True
            elif new_date == existing_date and new_id_num > existing_id_num:
                # Same date: use higher artifact ID (created more recently)
                should_replace = True

            if should_replace:
                logging.info(
                    f"Found newer demux artifact for flowcell {flowcell_key}: "
                    f"replacing {existing_art.id} ({existing_date}) with {demux_art.id} ({new_date})"
                )
                demux_arts_deduplicated[flowcell_key] = (demux_art, demux_art_parent)
            else:
                logging.info(
                    f"Skipping older demux artifact {demux_art.id} ({new_date}), "
                    f"keeping {existing_art.id} ({existing_date})"
                )

    if ongoing_demux_arts:
        ongoing_demux_steps_str = ", ".join(
            [
                f"'{art.parent_process.type.name}' ({art.parent_process.id})"
                for art in ongoing_demux_arts
            ]
        )
        logging.warning(
            f"Sample {sample.name} has demux artifacts in ongoing steps:"
            + f" {ongoing_demux_steps_str}. Finish them before proceeding."
        )

    logging.info(
        f"Deduplicated {len(demux_arts)} artifacts to {len(demux_arts_deduplicated)} "
        f"(keeping only latest per flowcell+lane)"
    )

    # Check for lane corrections (sample moved between lanes in different demux steps of same flowcell)
    flowcell_to_demux_steps = {}  # key: flowcell_id, value: {demux_step_id: set(lanes)}
    for (flowcell_id, lane), (demux_art, _) in demux_arts_deduplicated.items():
        if flowcell_id not in flowcell_to_demux_steps:
            flowcell_to_demux_steps[flowcell_id] = {}
        demux_step_id = demux_art.parent_process.id
        if demux_step_id not in flowcell_to_demux_steps[flowcell_id]:
            flowcell_to_demux_steps[flowcell_id][demux_step_id] = set()
        flowcell_to_demux_steps[flowcell_id][demux_step_id].add(lane)

    for flowcell_id, demux_steps in flowcell_to_demux_steps.items():
        if len(demux_steps) > 1 and "ONT" not in str(
            demux_steps
        ):  # Multiple demux steps
            all_lanes_across_steps = [lanes for lanes in demux_steps.values()]
            if len(all_lanes_across_steps) > 1:
                lanes_step1 = all_lanes_across_steps[0]
                lanes_step2 = all_lanes_across_steps[1]
                if lanes_step1 != lanes_step2:
                    raise ValueError(
                        f"ERROR: Sample '{sample.name}' appears in different lanes "
                        f"across demux steps of flowcell '{flowcell_id}': "
                        f"{dict(demux_steps)}. "
                        f"This indicates a lane correction in LIMS. "
                        f"Please verify the values and correct them manually if needed."
                    )

    tot_reads = 0
    flowcell_lane_list = []
    for demux_art, demux_art_parent in demux_arts_deduplicated.values():
        logging.info(
            f"Looking at '{demux_art.name}' ({demux_art.id}) of step"
            + f" '{demux_art.parent_process.type.name}'"
            + f" ({demux_art.parent_process.id})"
            + f"{' (ONGOING)...' if demux_art.parent_process.date_run is None else '...'}"
        )

        # Evaluate skip conditions
        if "# Reads" not in demux_art.udf:
            logging.warning("Missing or unpopulated UDF '# Reads', skipping.")
            continue

        if "Include reads" not in demux_art.udf:
            logging.warning(
                "Missing or unpopulated UDF 'Include_reads' filled, skipping."
            )
            continue

        if demux_art.udf["Include reads"] == "NO":
            logging.info("UDF 'Include reads' is set to 'NO', skipping.")
            continue

        assert demux_art.udf["Include reads"] == "YES"

        # Two cases to consider:
        if demux_art_parent.parent_process.type.name in SEQUENCING_STEPS:
            # Parent of demux artifact is a sequencing process output (ONT)
            seq_process = demux_art_parent.parent_process
        else:
            # Parent of demux artifact is a sequencing process input (Illumina, AVITI)
            seq_process = lims.get_processes(
                type=SEQUENCING_STEPS,
                inputartifactlimsid=demux_art_parent.id,
            )[0]

        # Check whether we are dealing with dual reads
        if (
            "Read 2 Cycles" in seq_process.udf
            and seq_process.udf["Read 2 Cycles"] is not None
        ):
            dual_reads = True
        else:
            dual_reads = False

        # Gather flowcell information
        if "ONT flow cell ID" in demux_art_parent.udf:
            # ONT
            ont_flowcell = demux_art_parent.udf["ONT flow cell ID"]
            if ont_flowcell not in flowcell_lane_list:
                flowcell_lane_list.append(ont_flowcell)
            if ont_flowcell not in summary[sample.name]:
                summary[sample.name][ont_flowcell] = set()
        else:
            # Illumina and AVITI
            flowcell_and_lane = "{}:{}".format(
                demux_art_parent.location[0].name,
                demux_art_parent.location[1].split(":")[0],
            )
            if flowcell_and_lane not in flowcell_lane_list:
                flowcell_lane_list.append(flowcell_and_lane)
            if demux_art_parent.location[0].name in summary[sample.name]:
                summary[sample.name][demux_art_parent.location[0].name].add(
                    demux_art_parent.location[1].split(":")[0]
                )
            else:
                summary[sample.name][demux_art_parent.location[0].name] = set(
                    demux_art_parent.location[1].split(":")[0]
                )

        # Aggregate reads to total
        if dual_reads:
            tallied_reads = float(demux_art.udf["# Reads"]) / 2
        else:
            tallied_reads = float(demux_art.udf["# Reads"])

        logging.info(
            f"Tallied {int(tallied_reads):,} {'dual' if dual_reads else 'single'} reads"
        )
        tot_reads += tallied_reads

    # Total is displayed as millions
    tot_reads_m = tot_reads / 1e6
    return tot_reads_m


def get_parent_inputs(art):
    input_arts = set()
    for input_output_tuple in art.parent_process.input_output_maps:
        if input_output_tuple[1]["uri"].id == art.id:
            input_arts.add(input_output_tuple[0]["uri"])

    return input_arts


if __name__ == "__main__":
    # Parse args
    parser = ArgumentParser(description=DESC)
    parser.add_argument("--pid", type=str, help="Lims ID for current Process")
    parser.add_argument("--log", type=str, help="Which log file slot to use")

    args = parser.parse_args()

    lims = Lims(BASEURI, USERNAME, PASSWORD)

    main(args)
