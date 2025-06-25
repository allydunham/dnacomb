#!/usr/bin/env python3
"""
Run benchmark tests on a large seque using a variety of settings.

Prints a benchmark results table to stdout and logging to stderr.
"""
import argparse
import os
import re
from utils import run_tool
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
                  library_counts=True, additional_args=None,
                  library_size="NA", read_length="NA",
                  outname="bench"):
    """
    Run benchmark suit for a given input
    """
    print(name, "... ", sep="", end="", flush=True)
    out, time = run_tool(f_file=f_file, r_file=r_file, lib_spec=lib_spec,
                         output=f"data/benchmark/{outname}", mode=mode, metric=metric,
                         verbose=True, no_cache=no_cache, sort=sort,
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

    print("Generating test data:")
    for fmt in ("fa", "fq"):
        if not os.path.exists(f"data/benchmark/{root}_grna.fq"):
            print(f"    Basic gRNA library ({fmt})... ", end="", flush=True)
            generate_test_data(lib_spec="config/grna.json", number=1000000, library_size=100,
                               output=f"data/benchmark/{root}_grna", seqformat=fmt)
            print("done")

    for lib in ("grna_sensor", "pegrna"):
        if not os.path.exists(f"data/benchmark/{root}_{lib}.fq"):
            print(f"    {lib} library... ", end="", flush=True)
            generate_test_data(lib_spec=f"config/{lib}.json", number=1000000, library_size=100,
                               output=f"data/benchmark/{root}_{lib}")
            print("done")

    lib_sizes = (1000, 10000, 100000)
    for size in lib_sizes:
        if not os.path.exists(f"data/benchmark/{root}_lib_{size}.fq"):
            print(f"    {size} sized library... ", end="", flush=True)
            _ = generate_library(lib_spec=f"config/grna.json", n=size,
                                 path=f"data/benchmark/{root}_lib_{size}")
            generate_test_data(lib_spec=f"data/benchmark/{root}_lib_{size}.json",
                               number=1000000, library_size=size,
                               output=f"data/benchmark/{root}_lib_{size}", recombination_rate=0.05,
                               contamination_rate=0.05, mismatch_rate=0.05, sub_rate=0.01,
                               indel_rate=0.001, truncation_rate=0.01)
            print("done")

    read_counts = (1000, 10000, 100000, 1000000, 10000000)
    for n in read_counts:
        if not os.path.exists(f"data/benchmark/{root}_reads_{n}.fq"):
            print(f"    gRNA #{n} seqs... ", end="", flush=True)
            generate_test_data(lib_spec=f"data/benchmark/{root}_lib_10000.json",
                               number=n, library_size=10000,
                               output=f"data/benchmark/{root}_reads_{n}", recombination_rate=0.05,
                               contamination_rate=0.05, mismatch_rate=0.05, sub_rate=0.01,
                               indel_rate=0.001, truncation_rate=0.01)
            print("done")

    print("\nRunning Benchmarks:")

    with open(f"data/benchmark/{root}.tsv", "w") as file:
        print(*BENCH_HEADERS, sep="\t", file=file)

        # Fa vs Fq
        for fmt in ("fq", "fa"):
            run_benchmark(f"Format ({fmt})", library_size=0,
                            lib_spec=f"config/grna.json",
                            read_length=156, f_file=f"data/benchmark/{root}_grna.{fmt}",
                            mode="inframe", library_counts=False, outfile=file,
                            outname=f"{root}_format_{fmt}")

        # Region Extraction
        for mode in ("full-read", "inframe", "pattern", "align"):
            for path, length in (("grna", 156), ("grna_sensor", 244), ("pegrna", 313)):
                if path == "pegrna" and mode == "inframe":
                    # can't do pegRNA variable regions with inframe
                    continue

                run_benchmark(f"Extraction ({path}/single/{mode})", library_size=0,
                              lib_spec=f"config/{path}.json",
                              read_length=length, f_file=f"data/benchmark/{root}_{path}.fq",
                              mode=mode, library_counts=False, outfile=file,
                              outname=f"{root}_extract_{mode}_{path}")

        # Paired end region extraction
        for mode in ("full-read", "inframe", "pattern", "align"):
            run_benchmark(f"Extraction (gRNA/paired/{mode})", library_size=0, read_length=156,
                          f_file=f"data/benchmark/{root}_grna_forward.fq",
                          r_file=f"data/benchmark/{root}_grna_reverse.fq",
                          lib_spec="config/grna.json", mode=mode, library_counts=False,
                          outfile=file, outname=f"{root}_paired_{mode}")

        # Distance metrics
        for metric in ("exact", "hamming", "bounded-levenshtein", "levenshtein"):
            run_benchmark(f"Library comparison ({metric})",
                          library_size=10000, read_length=156,
                          f_file=f"data/benchmark/{root}_lib_10000.fq",
                          lib_spec=f"data/benchmark/{root}_lib_10000.json",
                          mode="inframe", metric=metric,
                          library_counts=True, outfile=file,
                          outname=f"{args.output}_library_{metric}")

        # Library size
        for size in lib_sizes:
            run_benchmark(f"Library size ({size})",
                          library_size=size, read_length=156,
                          lib_spec=f"data/benchmark/{root}_lib_{size}.json",
                          f_file=f"data/benchmark/{root}_lib_{size}.fq",
                          mode="inframe", metric="hamming",
                          library_counts=True, outfile=file,
                          outname=f"{root}_lib_size_{size}")

        # Number of reads
        for n in read_counts:
            run_benchmark(f"Sequence file size ({n})",
                           library_size=10000, read_length=156,
                           lib_spec=f"data/benchmark/{root}_lib_10000.json",
                           f_file=f"data/benchmark/{root}_reads_{n}.fq",
                           mode="align", metric="hamming",
                           library_counts=True, outfile=file,
                           outname=f"{root}_file_size_{n}")

def parse_args(arg_list=None):
    """
    Parse and validate script arguments
    """
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("--output", "-o", default="bench",
                        help="Output name")

    return parser.parse_args(arg_list)

if __name__ == "__main__":
    main()
