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
import atexit
import logging
import re
import sys
from argparse import ArgumentParser
from collections import OrderedDict
from datetime import datetime

from genologics.config import BASEURI, PASSWORD, USERNAME
from genologics.entities import Process
from genologics.lims import Lims

from scilifelab_epps.epp import EppLogger

NGISAMPLE_PAT = re.compile("P[0-9]+_[0-9]+")


class ProgressTracker:
    def __init__(self, path):
        self.path = path
        self.handle = open(path, "a")
        self.write(f"Run started: {datetime.now().isoformat()}")
        atexit.register(self.close)

    def write(self, message):
        self.handle.write(f"{datetime.now().isoformat()} {message}\n")
        self.handle.flush()

    def close(self):
        if getattr(self, "handle", None) and not self.handle.closed:
            self.write("Run finished")
            self.handle.close()


def build_artifact_targets(artifacts, aggregate):
    targets = []
    for artifact in artifacts:
        if not artifact.samples:
            continue

        sample = artifact.samples[0]
        if not NGISAMPLE_PAT.search(sample.name):
            continue

        dest_obj = sample.artifact if aggregate else sample
        targets.append((artifact, dest_obj))
    return targets


def slice_chunk(items, chunk_size, chunk_index):
    if chunk_size <= 0:
        return items
    start = chunk_index * chunk_size
    end = start + chunk_size
    return items[start:end]


def build_projects(artifacts):
    projects = OrderedDict()
    for artifact in artifacts:
        for sample in artifact.samples:
            if sample.project:
                projects[sample.project.id] = sample.project
    return list(projects.values())


def queue_copy(source_obj, dest_obj, source_udf, dest_udf, pending_updates):
    """Queue a copy and defer saving so each destination is written once."""
    if source_udf not in source_obj.udf:
        logging.warning(f"Source UDF '{source_udf}' missing on {source_obj}")
        return False

    value = source_obj.udf.get(source_udf)

    if dest_obj.udf.get(dest_udf) == value:
        return False

    dest_obj.udf[dest_udf] = value
    pending_update = pending_updates.setdefault(dest_obj.id, {"obj": dest_obj, "logs": []})
    pending_update["logs"].append(f"Updated {dest_udf} on {dest_obj} from {source_udf}\n")
    return True


def flush_updates(pending_updates, target_name, changelog_f, progress_tracker):
    saved = 0
    failed = 0
    for idx, pending_update in enumerate(pending_updates.values(), start=1):
        dest_obj = pending_update["obj"]
        try:
            dest_obj.put()
            saved += 1
            if changelog_f:
                changelog_f.writelines(pending_update["logs"])
        except Exception as e:
            failed += 1
            logging.error(
                f"Failed to update {target_name} '{getattr(dest_obj, 'id', dest_obj)}': {e}"
            )
            progress_tracker.write(
                f"Failed to update {target_name} '{getattr(dest_obj, 'id', dest_obj)}': {e}"
            )
        if idx % 250 == 0:
            logging.info(f"Saved {idx} {target_name} objects so far")
            progress_tracker.write(f"Saved {idx} {target_name} objects so far")
    return saved, failed


def main(lims, args, epp_logger):
    progress_tracker = ProgressTracker(args.progress_log)
    progress_message = f"Progress log: {progress_tracker.path}"
    progress_tracker.write(progress_message)
    print(progress_message, file=sys.stderr)
    process = Process(lims, id=args.pid)
    progress_tracker.write(f"Loaded process {args.pid}")
    artifacts, _ = process.analytes()
    progress_tracker.write(f"Fetched {len(artifacts)} analytes from process")
    artifact_targets = build_artifact_targets(artifacts, args.aggregate)
    all_artifact_targets = len(artifact_targets)
    artifact_targets = slice_chunk(
        artifact_targets, args.art_chunk_size, args.art_chunk_index
    )
    projects = build_projects(artifacts) if args.proc_source_udf else []
    progress_tracker.write(
        f"Artifact chunk selection total={all_artifact_targets} chunk_size={args.art_chunk_size} "
        f"chunk_index={args.art_chunk_index} selected={len(artifact_targets)}"
    )

    logging.info(
        "Artifact chunk selection: total=%s chunk_size=%s chunk_index=%s selected=%s",
        all_artifact_targets,
        args.art_chunk_size,
        args.art_chunk_index,
        len(artifact_targets),
    )

    if args.status_changelog and not args.skip_status_changelog_prepend:
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
    art_saved = 0
    art_failed = 0
    art_pending_updates = OrderedDict()

    if args.art_source_udf:

        if not args.art_dest_udf:
            args.art_dest_udf = args.art_source_udf
        elif len(args.art_dest_udf) != len(args.art_source_udf):
            logging.error("art_source_udf and art_dest_udf lists of arguments are uneven.")
            progress_tracker.write("Error: art_source_udf and art_dest_udf lists are uneven")
            sys.exit(-1)

        for source_udf, dest_udf in zip(args.art_source_udf, args.art_dest_udf):
            progress_tracker.write(
                f"Queueing artifact copy source_udf='{source_udf}' dest_udf='{dest_udf}'"
            )
            for artifact, dest_obj in artifact_targets:
                updated = queue_copy(
                    artifact, dest_obj,
                    source_udf, dest_udf,
                    art_pending_updates,
                )

                if updated:
                    art_updates += 1
                else:
                    art_skipped += 1

        art_saved, art_failed = flush_updates(
            art_pending_updates, "artifact destination", changelog_f, progress_tracker
        )
        progress_tracker.write(
            f"Artifact flush complete updates={art_updates} skipped={art_skipped} "
            f"saved={art_saved} failed={art_failed}"
        )

    ############################################
    # PART 2 — Process ➜ Project
    ############################################

    proc_updates = 0
    proc_skipped = 0
    proc_saved = 0
    proc_failed = 0
    proc_pending_updates = OrderedDict()

    run_project_copy = bool(args.proc_source_udf)
    if args.proc_only_first_chunk and args.art_chunk_size > 0 and args.art_chunk_index > 0:
        run_project_copy = False
        logging.info(
            "Skipping process-to-project copy because --proc_only_first_chunk is set "
            "and this is chunk %s",
            args.art_chunk_index,
        )
        progress_tracker.write(
            f"Skipping process-to-project copy for chunk {args.art_chunk_index}"
        )

    if run_project_copy:

        if not args.proc_dest_udf:
            args.proc_dest_udf = args.proc_source_udf
        elif len(args.proc_dest_udf) != len(args.proc_source_udf):
            logging.error("proc_source_udf and proc_dest_udf lists of arguments are uneven.")
            progress_tracker.write("Error: proc_source_udf and proc_dest_udf lists are uneven")
            sys.exit(-1)

        for source_udf, dest_udf in zip(args.proc_source_udf, args.proc_dest_udf):
            progress_tracker.write(
                f"Queueing project copy source_udf='{source_udf}' dest_udf='{dest_udf}'"
            )
            for project in projects:
                updated = queue_copy(
                    process, project,
                    source_udf, dest_udf,
                    proc_pending_updates,
                )

                if updated:
                    proc_updates += 1
                else:
                    proc_skipped += 1

        proc_saved, proc_failed = flush_updates(
            proc_pending_updates, "project", changelog_f, progress_tracker
        )
        progress_tracker.write(
            f"Project flush complete updates={proc_updates} skipped={proc_skipped} "
            f"saved={proc_saved} failed={proc_failed}"
        )

    ############################################
    # Close changelog
    ############################################
    if changelog_f:
        changelog_f.close()

    ############################################
    # SUMMARY
    ############################################

    summary = (
        f"Artifact chunk total: {all_artifact_targets} | "
        f"Artifact chunk size: {args.art_chunk_size} | "
        f"Artifact chunk index: {args.art_chunk_index} | "
        f"Artifact chunk selected: {len(artifact_targets)} | "
        f"Artifact updates: {art_updates} | "
        f"Artifact skipped: {art_skipped} | "
        f"Artifact objects saved: {art_saved} | "
        f"Artifact save failures: {art_failed} | "
        f"Project updates: {proc_updates} | "
        f"Project skipped: {proc_skipped} | "
        f"Project objects saved: {proc_saved} | "
        f"Project save failures: {proc_failed}"
    )

    progress_tracker.write(summary)
    print(summary, file=sys.stderr)
    progress_tracker.close()


if __name__ == "__main__":

    parser = ArgumentParser(description="Production-grade combined UDF copier")

    parser.add_argument("--pid", required=True)
    parser.add_argument("--log")
    parser.add_argument("-c", "--status_changelog")
    parser.add_argument("--aggregate", action="store_true")
    parser.add_argument("--art_chunk_size", type=int, default=0)
    parser.add_argument("--art_chunk_index", type=int, default=0)
    parser.add_argument("--skip_status_changelog_prepend", action="store_true")
    parser.add_argument("--proc_only_first_chunk", action="store_true")
    parser.add_argument(
        "--progress_log",
        default=f"copy_field_project_summary_progress_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log",
    )

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
