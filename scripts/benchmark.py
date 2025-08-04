#!/usr/bin/env python3
"""
Run benchmark tests on a large seque using a variety of settings.

Prints a benchmark results table to stdout and logging to stderr.
"""
import argparse
import os
import re
from utils import run_tool
from itertools import product
from generate_test_data import generate_test_data, generate_library

REGION_RE = re.compile(
    "Extracted regions\: ([0-9]*) in ([0-9\.]*)([^ ]*) \| avg\. rate\: ([0-9\.]*)"
)

LIBRARY_RE = re.compile(
    "Compared combinations\: ([0-9]*)\/[0-9]* 100\% in ([0-9\.]*)([^ ]*) \| avg\. rate\: ([0-9\.]*)"
)

SUMMARY_RE = re.compile(
    "library matches\: ([0-9]*)\/[0-9]* 100\% in ([0-9\.]*)([^ ]*) \| avg\. rate\: ([0-9\.]*)"
)

def extract_time(time, unit):
    """
    Extract a time with s/ms/us/ns units
    """
    time = float(time)
    match unit:
        case "s":
            time = time
        case "ms":
            time = time/ 1e3
        case "ns":
            time = time / 1e9
        case _: # Mu has weird chars but is only other choice
            time = time / 1e6
    return time

BENCH_HEADERS = [
    "name", "fwd", "rev", "lib_spec", "mode", "metric",
    "no_cache", "sort", "group", "library_counts",
    "library_size", "read_length", "additional_args",
    "total_time", "reads", "region_time", "region_rate",
    "unique_regions", "library_time", "library_rate",
    "summary_size", "summary_time", "summary_rate"
]

def run_benchmark(name, outfile, f_file, r_file=None, lib_spec=None,
                  mode="inframe", metric="exact",
                  no_cache=False, sort=True, group=None,
                  threads=1,
                  library_counts=True, additional_args=None,
                  library_size="NA", read_length="NA",
                  outname="bench"):
    """
    Run benchmark suit for a given input
    """
    print(name, "... ", sep="", end="", flush=True)
    out, time = run_tool(f_file=f_file, r_file=r_file, lib_spec=lib_spec,
                         output=f"data/benchmark/{outname}", mode=mode, metric=metric,
                         verbose=True, no_cache=no_cache, sort=sort, threads=threads,
                         group=group, overwrite=True, library_counts=library_counts,
                         additional_args=additional_args, rm_output=True)
    if out.returncode != 0:
        print("failed:\n", out.stderr, sep="")
    else:
        print("done in ", round(time, 2), "s", sep="")
        with open(f"data/benchmark/{outname}.log", "w") as file:
            print(str(out.stderr), file=file)

    # Extract region count numbers
    reg = REGION_RE.search(str(out.stderr))
    reg_count = int(reg.group(1)) if reg is not None else "NA"
    reg_time = extract_time(reg.group(2), reg.group(3)) if reg is not None else "NA"
    reg_items = float(reg.group(4)) if reg is not None else "NA"

    # Extract library processing
    lib = LIBRARY_RE.search(str(out.stderr))
    lib_count = int(lib.group(1)) if lib is not None else "NA"
    lib_time = extract_time(lib.group(2), lib.group(3)) if lib is not None else "NA"
    lib_items = float(lib.group(4)) if lib is not None else "NA"

    # Extract library summarisation
    suma = LIBRARY_RE.search(str(out.stderr))
    suma_count = int(suma.group(1)) if suma is not None else "NA"
    suma_time = extract_time(suma.group(2), suma.group(3)) if suma is not None else "NA"
    suma_items = float(suma.group(4)) if suma is not None else "NA"

    print(name, f_file, r_file, lib_spec,
          mode, metric, no_cache, sort, group,
          library_counts, library_size, read_length,
          additional_args,
          time, reg_count, reg_time, reg_items,
          lib_count, lib_time, lib_items,
          suma_count, suma_time, suma_items,
          sep="\t", file=outfile)

def main():
    """Main"""
    args = parse_args()

    root = args.output

    os.makedirs("data/benchmark", exist_ok=True)

    print("\nGenerating data:")

    param_combs = product(
        ["grna", "grna_sensor", "pegrna"],
        [100, 1000, 10000]
    )

    for lib, size in param_combs:
        if not os.path.exists(f"data/benchmark/{root}_{lib}_{size}"):
            print(f"    {lib} {size} library... ", end="", flush=True)
            _ = generate_library(lib_spec=f"config/{lib}.json", n=size,
                                 path=f"data/benchmark/{root}_{lib}_{size}")
            print("done")

        if not os.path.exists(f"data/benchmark/{root}_{lib}_{size}.fq"):
            print(f"    {lib} {size} fastq... ", end="", flush=True)
            generate_test_data(lib_spec=f"data/benchmark/{root}_{lib}_{size}.json",
                               number=10000000, library_size=size,
                               output=f"data/benchmark/{root}_{lib}_{size}",
                               recombination_rate=0.01, contamination_rate=0.01,
                               mismatch_rate=0.01, truncation_rate=0.001,
                               sub_rate=0.001, indel_rate=0.0001)
            print("done")

    print("\nRunning Benchmarks:")

    with open(f"data/benchmark/{root}.tsv", "w") as file:
        print(*BENCH_HEADERS, sep="\t", file=file)

        param_combs = product(
            ["full-read", "inframe", "pattern", "align"],
            ["exact", "hamming", "bounded-levenshtein", "levenshtein"],
            [("grna", 156), ("grna_sensor", 244), ("pegrna", 313)],
            [100, 1000, 10000],
            [100000, 1000000, 10000000],
            [True, False],
            [True, False]
        )

        for mode, metric, (lib, read_length), lib_size, n_reads, nocache, paired in param_combs:
            if mode == "full-read" and not (metric == "exact" and lib_size == 100):
                continue

            if nocache and not mode == "align":
                continue

            name = f"mode:{mode}|metric:{metric}|lib:{lib}|lib_size:{lib_size}|reads:{n_reads}|threads:{args.threads}|cache:{not nocache}|paired:{paired}"

            if paired:
                run_benchmark(name, library_size=lib_size, read_length=read_length,
                              outname=f"{root}_{name}", outfile=file,
                              f_file=f"data/benchmark/{root}_{lib}_{lib_size}_forward.fq",
                              r_file=f"data/benchmark/{root}_{lib}_{lib_size}_reverse.fq",
                              lib_spec=f"data/benchmark/{root}_{lib}_{lib_size}.json",
                              mode=mode, metric=metric, library_counts=True, no_cache=nocache,
                              threads=args.threads, additional_args=["--max-reads", str(n_reads)])
            else:
                run_benchmark(name, library_size=lib_size, read_length=read_length,
                              outname=f"{root}_{name}", outfile=file,
                              f_file=f"data/benchmark/{root}_{lib}_{lib_size}.fq",
                              lib_spec=f"data/benchmark/{root}_{lib}_{lib_size}.json",
                              mode=mode, metric=metric, library_counts=True, no_cache=nocache,
                              threads=args.threads, additional_args=["--max-reads", str(n_reads)])

def parse_args(arg_list=None):
    """
    Parse and validate script arguments
    """
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("--output", "-o", default="bench",
                        help="Output name")

    parser.add_argument("--threads", "-t", default=1,
                        help="Number of threads to use")

    parser.add_argument("--global", "-g", action="store_true",
                        help="Use globally installed DNAComb rather than locallay compiled copy")

    return parser.parse_args(arg_list)

if __name__ == "__main__":
    main()
