#!/usr/bin/env python

from argparse import ArgumentParser
from genologics.config import BASEURI, PASSWORD, USERNAME
from genologics.entities import Process
from genologics.lims import Lims
import sys
import traceback

def main(lims, args):
    currentStep = Process(lims, id=args.pid)
    driver_file_out = None
    log_file_out = None
    
    # Standard log file name
    log_filename = "epp_execution_log.txt"
    
    # Start capturing all print statements and errors to a file
    with open(log_filename, "w") as log_f:
        try:
            log_f.write(f"Script started for process {args.pid}\n")
            
            # --- 1. FIND FILE SLOTS ---
            for output in currentStep.all_outputs():
                if output.name == "Driver File":
                    driver_file_out = output
                elif output.name == "EPP Log": # Make sure this matches LIMS name!
                    log_file_out = output

            # --- 2. DATA COLLECTION ---
            # We check both Inputs and Outputs to be 100% safe
            artifacts = currentStep.all_inputs() + currentStep.all_outputs()
            ar_driver = {}
            valid_cols = set()

            for art in artifacts:
                if art == driver_file_out or art == log_file_out:
                    continue
                
                # Check for location. location[1] is the well (e.g., A:1)
                if art.location and art.location[1]:
                    loc = art.location[1]
                    well = loc.replace(":", "")
                    row = loc.split(":")[0]
                    
                    if well not in ar_driver:
                        valid_cols.add(row)
                        ar_driver[well] = art.name
                        log_f.write(f"INFO: Mapped {art.name} to {well}\n")

            # --- 3. LAYOUT GENERATION ---
            if not valid_cols:
                log_f.write("ERROR: No samples found with valid plate locations.\n")
            else:
                driver = []
                col_idx = -1
                for column in sorted(list(valid_cols)):
                    col_idx += 1
                    for i in range(1, 13):
                        location = f"{column}{i}"
                        name = ar_driver.get(location, "ladder" if i == 12 else "")
                        driver.append((col_idx * 12 + i, location, name))

                # --- 4. FILE WRITING ---
                csv_filename = "frag_an_driver.csv"
                with open(csv_filename, "w") as f:
                    for line in driver:
                        f.write(f"{line[0]},{line[1]},{line[2]}\n")
                
                # Upload the driver file
                if driver_file_out:
                    lims.upload_new_file(driver_file_out, csv_filename)
                    log_f.write("SUCCESS: Driver file uploaded.\n")

        except Exception as e:
            # If anything crashes, write the traceback to the log file
            log_f.write("\n--- SCRIPT CRASHED ---\n")
            log_f.write(traceback.format_exc())
            print(f"Script failed. Check {log_filename} for details.")

    # --- 5. FINAL LOG UPLOAD ---
    # This happens outside the try/except to ensure it always uploads
    if log_file_out:
        lims.upload_new_file(log_file_out, log_filename)