#!/usr/bin/env python

from argparse import ArgumentParser
from genologics.config import BASEURI, PASSWORD, USERNAME
from genologics.entities import Process
from genologics.lims import Lims
import sys

DESC = """EPP used to create csv files to drive the fragment analyzer"""

def main(lims, args):
    currentStep = Process(lims, id=args.pid)
    driver_file_out = None
    driver = []
    ar_driver = {}
    valid_cols = set()

    # --- 1. DATA COLLECTION ---
    for output in currentStep.all_outputs():
        # Check for the file placeholder
        if output.name == "Driver File":
            driver_file_out = output
        
        # Check for Samples/Pools (Analytes)
        elif output.type == "Analyte" and output.location[1]:
            location = output.location[1]
            row_letter = location.split(":")[0]
            well_key = location.replace(":", "")
            
            valid_cols.add(row_letter)
            # FIX: Using output.name ensures Pools show their Pool Name
            ar_driver[well_key] = output.name
            # Debug print to help you see it in the terminal
            print(f"DEBUG: Found {output.name} at {well_key}")

    # --- 2. VALIDATION ---
    if not driver_file_out:
        print("ERROR: Could not find an output artifact named 'Driver File'. Check LIMS configuration.")
        sys.exit(1)
    
    if not valid_cols:
        print("ERROR: No samples with valid plate locations were found.")
        # We still want to upload something to avoid a zero-byte hang
        with open("frag_an_driver.csv", "w") as f:
            f.write("0,No_Samples_Found,Empty\n")
        lims.upload_new_file(driver_file_out, "frag_an_driver.csv")
        return

    # --- 3. LAYOUT GENERATION ---
    col_idx = -1
    for column in sorted(list(valid_cols)):
        col_idx += 1
        for i in range(1, 13):
            location = f"{column}{i}"
            # Check dictionary for sample, default to ladder for col 12, else empty
            sample_name = ar_driver.get(location, "ladder" if i == 12 else "")
            
            driver.append(
                (
                    col_idx * 12 + i,
                    location,
                    sample_name,
                )
            )

    # --- 4. FILE WRITING AND UPLOAD ---
    filename = "frag_an_driver.csv"
    with open(filename, "w") as f:
        for line in driver:
            f.write(f"{line[0]},{line[1]},{line[2]}\n")

    lims.upload_new_file(driver_file_out, filename)
    print(f"DONE: Successfully uploaded driver file with {len(driver)} lines.")

if __name__ == "__main__":
    parser = ArgumentParser(description=DESC)
    parser.add_argument("--pid", help="Lims id for current Process")
    args = parser.parse_args()

    lims = Lims(BASEURI, USERNAME, PASSWORD)
    lims.check_version()
    main(lims, args)