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
from os.path import join

def try_open_notebook(stdout_str):
  if "Or copy and paste one of these URLs:" in stdout_str:
    url_i = stdout_str.index("Or copy and paste one of these URLs:")
    if(len(stdout_str) == url_i + 122):
      url = stdout_str[url_i+45:url_i+122]
      
      # Get the name of most recently created file, and open directly to that notebook in JupyterLab.
      notebook_files = [ f for f in os.listdir(join(".snakemake", "scripts")) if f.endswith(".py.ipynb") ]
      
      if len(notebook_files) > 0:
        notebook_files = sorted(notebook_files, key=lambda f: os.path.getmtime(join(".snakemake", "scripts", f)), reverse=True)
        most_recent_notebook_file = notebook_files[0]
        url_path = f"/lab/tree/{most_recent_notebook_file}?token="
        url = url.replace("/?token=", url_path)
      
      subprocess.run(f"open {url}", shell=True)

def write_stdout(proc):
  s = ""
  for c in iter(lambda: proc.stdout.read(1), b''): 
    sys.stdout.buffer.write(c)
    sys.stdout.flush()
    try:
      s += c.decode("utf-8")
      try_open_notebook(s)
    except UnicodeDecodeError:
      pass

if __name__ == "__main__":
  
  parser = argparse.ArgumentParser(description="Edit a snakemake notebook for the first output file of a particular rule.")
  parser.add_argument("-f", "--snakefile", required=True, help="The name of the snakefile.")
  parser.add_argument("-r", "--rule", required=False, help="The name of the rule.")
  parser.add_argument("-s", "--suffix", required=False, help="The suffix for an output file.")
  args = parser.parse_args()
  
  summary = subprocess.run(f"snakemake --snakefile {args.snakefile} --summary", shell=True, capture_output=True)
  df = pd.read_csv(io.StringIO(summary.stdout.decode("utf-8")), sep='\t')
  
  if args.rule is not None:
    rule_df = df.loc[df["rule"] == args.rule].reset_index(drop=True)
  elif args.suffix is not None:
    rule_df = df.loc[df["output_file"].str.endswith(args.suffix)]
  else:
    raise ValueError("Provide the --rule or --suffix argument.")

  if rule_df.shape[0] > 0:
    row = rule_df.iloc[0]
    output_file = row["output_file"]
    proc = subprocess.Popen(
      f"snakemake --snakefile {args.snakefile} -j 1 --edit-notebook {output_file}",
      shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, preexec_fn=os.setsid
    )
    def cleanup():
      try:
        os.killpg(os.getpgid(proc.pid), signal.SIGINT)
        write_stdout(proc)
      except ProcessLookupError:
        pass
    atexit.register(cleanup)
    write_stdout(proc)
  else:
    print("No output files were found for the specified rule or suffix.")