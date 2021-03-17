#!/usr/bin/env python3
import subprocess
import argparse
import pandas as pd
import numpy as np
import io
import os
import sys
import atexit
import signal
import time

def write_stdout(proc):
  for c in iter(lambda: proc.stdout.read(1), b''): 
    sys.stdout.buffer.write(c)
    sys.stdout.flush()

if __name__ == "__main__":
  
  parser = argparse.ArgumentParser(description="Edit a snakemake notebook for the first output file of a particular rule.")
  parser.add_argument("-r", "--rule", required=True, help="The name of the rule.")
  args = parser.parse_args()
  
  summary = subprocess.run("snakemake --filegraph | dot -Txdot_json", shell=True, capture_output=True)
  df = pd.read_csv(io.StringIO(summary.stdout.decode("utf-8")), sep='\t')

  rule_df = df.loc[df["rule"] == args.rule].reset_index(drop=True)

  if rule_df.shape[0] > 0:
    row = rule_df.iloc[0]
    output_file = row["output_file"]
    proc = subprocess.Popen(
      f"snakemake -j 1 --edit-notebook {output_file}",
      shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, preexec_fn=os.setsid
    )
    def cleanup():
      try:
        os.killpg(os.getpgid(proc.pid), signal.SIGINT)
        write_stdout(proc)
        time.sleep(2)
        os.killpg(os.getpgid(proc.pid), signal.SIGTERM)
      except ProcessLookupError:
        pass
    atexit.register(cleanup)
    write_stdout(proc)
  else:
    print(df)
    print("No output files found for the specified rule.")