#!/usr/bin/env python

from argparse import ArgumentParser

from genologics.config import BASEURI, PASSWORD, USERNAME
from genologics.entities import Process
from genologics.lims import Lims

DESC = """EPP used to create csv files to drive the fragment analyzer"""


def main(lims, args):
    currentStep = Process(lims, id=args.pid)
    driver_file_out = None
    driver = []
    ar_driver = {}
    valid_cols = set()
    for output in currentStep.all_outputs():
        # 1. Catch the placeholder
        if output.name == "Driver File":
            driver_file_out = output
            continue
            
        # 2. Try to get the location
        try:
            location = output.location[1]
        except (AttributeError, IndexError):
            location = None

        # 3. Process if it's a sample/pool with a position
        if location:
            well = location.replace(":", "")
            row_letter = location.split(":")[0]
            
            valid_cols.add(row_letter)
            ar_driver[well] = output.name
            print(f"Successfully added {output.name} at {well}")

    col_idx = -1
    for column in sorted(list(valid_cols)):
        col_idx += 1
        for i in range(1, 13):
            location = f"{column}{i}"
            driver.append(
                (
                    col_idx * 12 + i,
                    location,
                    ar_driver.get(location, "ladder" if i == 12 else ""),
                )
            )

    with open("frag_an_driver.csv", "w") as f:
        for line in driver:
            f.write(f"{line[0]},{line[1]},{line[2]}\n")

    lims.upload_new_file(driver_file_out, "frag_an_driver.csv")


if __name__ == "__main__":
    parser = ArgumentParser(description=DESC)
    parser.add_argument("--pid", help="Lims id for current Process")
    args = parser.parse_args()

    lims = Lims(BASEURI, USERNAME, PASSWORD)
    lims.check_version()
    main(lims, args)
