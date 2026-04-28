#!/usr/bin/env python3
"""
Run benchmark tests on a large seque using a variety of settings.

Prints a benchmark results table to stdout and logging to stderr.
"""
import argparse
import sys
import os
import re
from utils import run_tool
from itertools import product
from generate_test_data import generate_test_data, generate_library

EXTRACT_RE = re.compile(
    "Extracted regions\: ([0-9]*) in ([0-9\.]*)([^ ]*) \| avg\. rate\: ([0-9\.]*)"
)

MATCH_RE = re.compile(
    "Matched regions\: ([0-9]*)\/[0-9]* 100\% in ([0-9\.]*)([^ ]*) \| avg\. rate\: ([0-9\.]*)"
)

COMBS_RE = re.compile(
    "Compared combinations\: ([0-9]*)\/[0-9]* 100\% in ([0-9\.]*)([^ ]*) \| avg\. rate\: ([0-9\.]*)"
)

SUMMARY_RE = re.compile(
    "Summarised library matches\: ([0-9]*)\/[0-9]* 100\% in ([0-9\.]*)([^ ]*) \| avg\. rate\: ([0-9\.]*)"
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
    "no_cache", "sort", "group", "library",
    "library_size", "read_length", "skip_variants", "additional_args",
    "total_time", "reads", "extraction_time", "extraction_rate",
    "unique_regions", "region_matching_time", "region_matching_rate",
    "unique_combinations", "combination_time", "combination_rate",
    "summary_size", "summary_time", "summary_rate"
]

def run_benchmark(name, outfile, f_file, r_file=None, lib_spec=None,
                  mode="inframe", metric="exact", skip_variants=False,
                  no_cache=False, sort=True, group=None,
                  threads=1, library=None, additional_args=None,
                  library_size="NA", read_length="NA",
                  outname="bench", path=None):
    """
    Run benchmark suit for a given input
    """
    print(name, "... ", sep="", end="", flush=True)
    out, time = run_tool(f_file=f_file, r_file=r_file, lib_spec=lib_spec,
                         output=f"data/benchmark/{outname}",
                         library=library, mode=mode, metric=metric,
                         verbose=True, no_cache=no_cache, sort=sort, threads=threads,
                         group=group, overwrite=True,
                         additional_args=additional_args, rm_output=True,
                         path=path)
    if out.returncode != 0:
        print("failed:\n", out.stderr, sep="")
    else:
        print("done in ", round(time, 2), "s", sep="")
        with open(f"data/benchmark/{outname}.log", "w") as file:
            print(str(out.stderr), file=file)

    # Extract region extraction processing (possibly multiple options as threaded)
    extraction = EXTRACT_RE.findall(str(out.stderr))
    extraction_count = 0
    extraction_time = 0
    extraction_items = 0
    for i in extraction:
        extraction_count += int(i[0])
        extraction_time = max(extract_time(i[1], i[2]), extraction_time)
        extraction_items += float(i[3])

    extraction_count = extraction_count if extraction_count > 0 else "NA"
    extraction_time = extraction_time if extraction_time > 0 else "NA"
    extraction_items = extraction_items if extraction_items > 0 else "NA"

    # Extract region matching processing
    match = MATCH_RE.search(str(out.stderr))
    match_count = int(match.group(1)) if match is not None else "NA"
    match_time = extract_time(match.group(2), match.group(3)) if match is not None else "NA"
    match_items = float(match.group(4)) if match is not None else "NA"

    # Extract combination processing
    lib = COMBS_RE.search(str(out.stderr))
    lib_count = int(lib.group(1)) if lib is not None else "NA"
    lib_time = extract_time(lib.group(2), lib.group(3)) if lib is not None else "NA"
    lib_items = float(lib.group(4)) if lib is not None else "NA"

    # Extract library summarisation
    suma = SUMMARY_RE.search(str(out.stderr))
    suma_count = int(suma.group(1)) if suma is not None else "NA"
    suma_time = extract_time(suma.group(2), suma.group(3)) if suma is not None else "NA"
    suma_items = float(suma.group(4)) if suma is not None else "NA"

    print(name, f_file, r_file, lib_spec,
          mode, metric, no_cache, sort, group,
          "".join(library), library_size, read_length,
          skip_variants, additional_args, time,
          extraction_count, extraction_time, extraction_items,
          match_count, match_time, match_items,
          lib_count, lib_time, lib_items,
          suma_count, suma_time, suma_items,
          sep="\t", file=outfile)

def main():
    """Main"""
    args = parse_args()

    os.makedirs(args.root, exist_ok=True)

    inroot = f"{args.root}/{args.input}"
    outroot = f"{args.root}/{args.output}"
    outname = args.output

    print("\nGenerating data:")

    param_combs = product(
        ["grna", "grna_sensor", "pegrna"],
        [100, 1000, 10000]
    )

    for lib, size in param_combs:
        if not os.path.exists(f"{inroot}_{lib}_{size}"):
            print(f"    {lib} {size} library... ", end="", flush=True)
            _ = generate_library(lib_spec=f"config/{lib}.json", n=size,
                                 path=f"{inroot}_{lib}_{size}")
            print("done")

        if not os.path.exists(f"{inroot}_{lib}_{size}.fq"):
            print(f"    {lib} {size} fastq... ", end="", flush=True)
            generate_test_data(lib_spec=f"config/{lib}.json",
                               libraries=[f"{inroot}_{lib}_{size}.tsv"],
                               number=1000000, library_size=size,
                               output=f"{inroot}_{lib}_{size}",
                               recombination_rate=0.001, contamination_rate=0.01,
                               mismatch_rate=0.01, truncation_rate=0.001,
                               sub_rate=0.0001, indel_rate=0.00001)
            print("done")

    if args.gen_only:
        sys.exit(0)

    print("\nRunning Benchmarks:")

    with open(f"{outroot}.tsv", "w") as file:
        print(*BENCH_HEADERS, sep="\t", file=file)

        param_combs = product(
            ["full-read", "inframe", "pattern", "align"],
            ["exact", "hamming", "bounded-levenshtein", "levenshtein"],
            [("grna", 156), ("grna_sensor", 244), ("pegrna", 313)],
            [100, 1000, 10000],
            [10000, 100000, 1000000],
            [True, False],
            [True, False],
            [True, False]
        )

        for mode, metric, (lib, read_length), lib_size, n_reads, nocache, paired, skip_variants in param_combs:
            if mode == "full-read" and not (metric == "exact" and lib_size == 100):
                continue

            if mode == "inframe" and lib == "pegrna":
                continue

            # Only do library comparison for the best estimates from align
            if metric != "exact" and not (mode == "align" and not nocache and n_reads < 10000000):
                continue

            name = f"mode:{mode}|metric:{metric}|lib:{lib}|lib_size:{lib_size}|reads:{n_reads}|threads:{args.threads}|cache:{not nocache}|paired:{paired}|variants:{not skip_variants}"

            if paired:
                run_benchmark(name, library_size=lib_size, read_length=read_length,
                              outname=f"{outname}_{name}", outfile=file,
                              f_file=f"{inroot}_{lib}_{lib_size}_forward.fq",
                              r_file=f"{inroot}_{lib}_{lib_size}_reverse.fq",
                              lib_spec=f"config/{lib}.json",
                              mode=mode, metric=metric, library=[f"{inroot}_{lib}_{lib_size}.tsv"], no_cache=nocache, skip_variants=skip_variants,
                              threads=args.threads, additional_args=["--max-reads", str(n_reads)],
                              path=args.path)
            else:
                run_benchmark(name, library_size=lib_size, read_length=read_length,
                              outname=f"{outname}_{name}", outfile=file,
                              f_file=f"{inroot}_{lib}_{lib_size}.fq",
                              lib_spec=f"config/{lib}.json",
                              mode=mode, metric=metric, library=[f"{inroot}_{lib}_{lib_size}.tsv"], no_cache=nocache, skip_variants=skip_variants,
                              threads=args.threads, additional_args=["--max-reads", str(n_reads)],
                              path=args.path)

def parse_args(arg_list=None):
    """
    Parse and validate script arguments
    """
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("--input", "-i", default="bench",
                        help="Root sequence input name. Benchmarks with the same name share input, so take care the sequences are generated first")

    parser.add_argument("--gen_only", "-g", action="store_true",
                        help="Generate input sequences without running the benchmark, for instance to then run multiple repeats on the input")

    parser.add_argument("--output", "-o", default="bench",
                        help="Output name")

    parser.add_argument("--root", "-r", default="data/benchmark",
                        help="Root folder to work in ")

    parser.add_argument("--threads", "-t", default=1, type=int,
                        help="Number of threads to use")

    parser.add_argument("--path", "-p",
                        help="Use DNAComb instance found at path instead of /target/release/dnacomb")

    return parser.parse_args(arg_list)

if __name__ == "__main__":
    main()
