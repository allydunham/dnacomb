#!/usr/bin/env python3
"""
Utility functions for testing scripts
"""
import os
import time
import subprocess

def run_tool(f_file, r_file=None, lib_spec=None, output="test", mode="inframe",
             metric="exact", verbose=True, no_cache=False, sort=True,
             group=None, overwrite=True, library=None, threads=1,
             path=None, additional_args=None, rm_output=False):
    """
    Run the tool, returning stdout/err and the completed process and a time in seconds
    """
    args = [
        "target/release/dnacomb" if path is None else path,
        "--output", output,
        "--mode", mode,
        "--distance-metric", metric,
    ]

    if lib_spec is not None:
        args.append("--library-spec")
        args.append(lib_spec)

    if verbose:
        args.append("--verbose")

    if sort:
        args.append("--sort")

    if no_cache:
        args.append("--no-cache")

    if overwrite:
        args.append("--overwrite")

    if library is not None:
        args.append("--library")
        args.extend(library)

    if group is not None:
        args.append("--group")
        args.append(group)

    if threads > 1:
        args.append("--threads")
        args.append(str(threads))

    if additional_args is not None:
        args.extend(additional_args)

    # Add this to make sure F/R are properly identified
    args.append("--")

    args.append(f_file)
    if r_file is not None:
        args.append(r_file)

    start = time.monotonic()
    out = subprocess.run(args, capture_output=True)
    end = time.monotonic()

    if rm_output:
        if os.path.exists(f"{output}.counts.tsv"):
            os.remove(f"{output}.counts.tsv")
        if os.path.exists(f"{output}.library_counts.tsv"):
            os.remove(f"{output}.library_counts.tsv")
        if os.path.exists(f"{output}.summary.tsv"):
            os.remove(f"{output}.summary.tsv")

    return (out, end - start)
