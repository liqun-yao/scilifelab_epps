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
PLATE_PAT = re.compile(r"P(\d+)$")  # Extract plate number from project ID


# Verify sample IDs
def verify_sample_ids(lims, project):
    """Validate project sample IDs for format, count, and sequence."""
    message = []
    
    # Get all samples in the project
    samples = lims.get_samples(projectname=project.name)
    
    if not samples:
        message.append(f"SAMPLE COUNT WARNING: No samples found for project {project.id}")
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
                f"SAMPLE NAME WARNING: Sample ID {sample_id} does not match project ID {project.id}"
            )

    # Check count consistency
    if len(ngi_ids) != len(customer_names):
        message.append(
            f"SAMPLE COUNT WARNING: Mismatch between NGI Sample IDs ({len(ngi_ids)}) "
            f"and customer sample names ({len(customer_names)})"
        )

    # Check first sample numbering based on plate number
    plate_match = PLATE_PAT.search(project.id)
    if plate_match and ngi_ids:
        plate_num = plate_match.group(1)
        first_sample = ngi_ids[0]
        expected_suffixes = f"_{plate_num}(001|01)$"
        if not re.search(expected_suffixes, first_sample):
            message.append(
                f"SAMPLE SEQUENCE WARNING: First sample {first_sample} should end with _{plate_num}001 or _{plate_num}01"
            )

    # Check for gaps in sample numbering
    if ngi_ids:
        try:
            suffixes = sorted([int(s.split("_")[1]) for s in ngi_ids])
            for curr, next_val in zip(suffixes, suffixes[1:]):
                if next_val - curr != 1:
                    message.append(
                        f"SAMPLE SEQUENCE WARNING: Gap detected in sample numbering between {curr} and {next_val}. "
                        f"Verify missing samples in the uploaded CSV file."
                    )
        except ValueError:
            pass  # Skip gap detection if suffix extraction fails

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
