DESC = """EPP script to copy User Defined Fields (UDFs) across multiple levels in Clarity LIMS:
1. From Analyte level to Submitted Sample level (or Sample-Artifact if --aggregate is used).
2. From Process level to Project level.

The script is optimized for high-throughput projects (e.g., 1900+ samples) by using 
unified API sessions and optimized file I/O to prevent server timeouts. It filters 
for valid NGI samples using the pattern 'P[0-9]+_[0-9]+' to skip controls.

Execution:
Can be executed in the background or triggered by a user pressing a "blue button" 
during step transitions.

Log Files:
1. Status Changelog (-c): Detailed audit trail including technician, date, and 
   specific UDF value changes for every sample and project.
2. Execution Log (--log): General runtime information, performance metrics, and 
   error handling.

Error handling:
If a source UDF is blank or undefined for any input or process, the script logs 
the omission and continues with the remaining items without failing the entire 
automation.
Written by Liqun Yao, Science for Life Laboratory, Stockholm, Sweden
"""
import logging
import re
import sys
from argparse import ArgumentParser

from genologics.config import BASEURI, PASSWORD, USERNAME
from genologics.entities import Process
from genologics.lims import Lims

from scilifelab_epps.epp import EppLogger

NGISAMPLE_PAT = re.compile("P[0-9]+_[0-9]+")



def safe_copy(source_obj, dest_obj, source_udf, dest_udf, changelog_f):
    """Best-effort safe copy. Skips if missing or unchanged."""
    if source_udf not in source_obj.udf:
        logging.warning(
            f"Source UDF '{source_udf}' missing on {source_obj}"
        )
        return False

    value = source_obj.udf.get(source_udf)

    if dest_obj.udf.get(dest_udf) == value:
        return False  # no update needed

    dest_obj.udf[dest_udf] = value
    dest_obj.put()

    changelog_f.write(
        f"Updated {dest_udf} on {dest_obj} from {source_udf}\n"
    )
    return True


def main(lims, args, epp_logger):

    process = Process(lims, id=args.pid)
    artifacts, _ = process.analytes()

    if args.status_changelog:
        epp_logger.prepend_old_log(args.status_changelog)

    ############################################
    # Open changelog ONCE (performance critical)
    ############################################
    changelog_f = open(args.status_changelog, "a") if args.status_changelog else None

    ############################################
    # PART 1 — Artifact ➜ Sample
    ############################################

    art_updates = 0
    art_skipped = 0

    if args.art_source_udf:

        if not args.art_dest_udf:
            args.art_dest_udf = args.art_source_udf

        for source_udf, dest_udf in zip(args.art_source_udf, args.art_dest_udf):

            for artifact in artifacts:

                if not artifact.samples:
                    continue

                sample = artifact.samples[0]

                # Only NGI samples
                if not NGISAMPLE_PAT.search(sample.name):
                    continue

                dest_obj = sample.artifact if args.aggregate else sample

                updated = safe_copy(
                    artifact, dest_obj,
                    source_udf, dest_udf,
                    changelog_f
                )

                if updated:
                    art_updates += 1
                else:
                    art_skipped += 1

    ############################################
    # PART 2 — Process ➜ Project
    ############################################

    proc_updates = 0
    proc_skipped = 0

    if args.proc_source_udf:

        if not args.proc_dest_udf:
            args.proc_dest_udf = args.proc_source_udf

        # Collect unique projects efficiently
        projects = {
            s.project
            for a in artifacts
            for s in a.samples
            if s.project
        }

        for source_udf, dest_udf in zip(args.proc_source_udf, args.proc_dest_udf):

            for project in projects:

                updated = safe_copy(
                    process, project,
                    source_udf, dest_udf,
                    changelog_f
                )

                if updated:
                    proc_updates += 1
                else:
                    proc_skipped += 1

    ############################################
    # Close changelog
    ############################################
    if changelog_f:
        changelog_f.close()

    ############################################
    # SUMMARY
    ############################################

    summary = (
        f"Artifact updates: {art_updates} | "
        f"Artifact skipped: {art_skipped} | "
        f"Project updates: {proc_updates} | "
        f"Project skipped: {proc_skipped}"
    )

    print(summary, file=sys.stderr)


if __name__ == "__main__":

    parser = ArgumentParser(description="Production-grade combined UDF copier")

    parser.add_argument("--pid", required=True)
    parser.add_argument("--log")
    parser.add_argument("-c", "--status_changelog")
    parser.add_argument("--aggregate", action="store_true")

    # Artifact arguments
    parser.add_argument("--art_source_udf", nargs="*")
    parser.add_argument("--art_dest_udf", nargs="*")

    # Process arguments
    parser.add_argument("--proc_source_udf", nargs="*")
    parser.add_argument("--proc_dest_udf", nargs="*")

    args = parser.parse_args()

    lims = Lims(BASEURI, USERNAME, PASSWORD)
    lims.check_version()

    with EppLogger(log_file=args.log, lims=lims, prepend=True) as epp_logger:
        main(lims, args, epp_logger)