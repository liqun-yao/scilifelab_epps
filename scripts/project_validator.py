#!/usr/bin/env python

import re
import sys
from argparse import ArgumentParser

from genologics.config import BASEURI, PASSWORD, USERNAME
from genologics.entities import Project
from genologics.lims import Lims

DESC = """EPP used to validate a project
Author: Chuan Wang, Science for Life Laboratory, Stockholm, Sweden
"""

# Pre-compile regexes in global scope:
NGISAMPLE_PAT = re.compile("P[0-9]+_[0-9]+")


# Verify sample IDs
def verify_sample_ids(lims, project):
    """Validate project sample IDs for format, count, and sequence."""
    message = []

    # Get all samples in the project
    samples = lims.get_samples(projectname=project.name)

    if not samples:
        message.append(
            f"SAMPLE COUNT WARNING: No samples found for project {project.id}"
        )
        return message

    ngi_ids = []
    customer_names = []

    # Validate sample name format and collect data
    for sample in sorted(samples, key=lambda s: s.name):
        sample_id = sample.name
        customer_name = sample.udf.get("Customer Name")

        ngi_ids.append(sample_id)
        if customer_name:
            customer_names.append(customer_name)

        # Validate format and prefix match
        if not NGISAMPLE_PAT.search(sample_id):
            message.append(f"SAMPLE NAME WARNING: Bad sample ID format {sample_id}")
        elif sample_id.split("_")[0] != project.id:
            message.append(
                f"SAMPLE NAME WARNING: Sample ID {sample_id} does not match "
                f"project ID {project.id}"
            )

    # Check count consistency
    if len(ngi_ids) != len(customer_names):
        message.append(
            f"SAMPLE COUNT WARNING: Mismatch between NGI Sample IDs ({len(ngi_ids)}) "
            f"and customer sample names ({len(customer_names)})"
        )

    # Validate sample numbering: group by plate digit, check first sample and gaps
    if ngi_ids:
        try:
            # Group samples by plate digit (first digit of suffix: _1XXX, _2XXX, etc.)
            plates = {}
            for sample_id in ngi_ids:
                sample_suffix = sample_id.split("_")[1] if "_" in sample_id else None
                if sample_suffix:
                    plate_digit = sample_suffix[0]
                    suffix_num = int(sample_suffix)
                    if plate_digit not in plates:
                        plates[plate_digit] = []
                    plates[plate_digit].append(suffix_num)

            # Validate each plate: first sample must be X001 or X01, check gaps
            for plate_digit in sorted(plates.keys()):
                plate_suffixes = sorted(plates[plate_digit])

                # Check first sample starts with X001 or X01
                expected_values = (
                    int(f"{plate_digit}001"),
                    int(f"{plate_digit}01"),
                )
                if plate_suffixes[0] not in expected_values:
                    message.append(
                        f"SAMPLE SEQUENCE WARNING: Plate {plate_digit} first "
                        f"sample should be _{plate_digit}001 or _{plate_digit}01, "
                        f"but got _{plate_suffixes[0]}"
                    )

                # Check gaps within plate
                for curr, next_val in zip(plate_suffixes, plate_suffixes[1:]):
                    if next_val - curr != 1:
                        message.append(
                            f"SAMPLE SEQUENCE WARNING: Gap detected in plate "
                            f"{plate_digit} numbering between {curr} and "
                            f"{next_val}. Verify missing samples in the "
                            f"uploaded CSV file."
                        )
        except (ValueError, IndexError):
            pass  # Skip validation if suffix extraction fails

    return message


def main(lims, pid):
    """Validate a project and exit with appropriate status code."""
    message = []
    project = Project(lims, id=pid)

    # Validate sample IDs
    message += verify_sample_ids(lims, project)

    if not message:
        print(f"No issue detected for project {pid}")
    else:
        sys.stderr.write("; ".join(message))
        sys.exit(2)


if __name__ == "__main__":
    parser = ArgumentParser(description=DESC)
    parser.add_argument("--pid", help="Project ID for current Project")
    args = parser.parse_args()

    lims = Lims(BASEURI, USERNAME, PASSWORD)
    lims.check_version()
    main(lims, args.pid)
