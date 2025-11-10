#!/usr/bin/env python3
"""
Run end-to-end tests using a variety of input sequence mutation levels
"""
import sys
import os
from itertools import product
from utils import run_tool
from generate_test_data import generate_test_data

def create_test_runner(results):
    def f(name, critical=False, expected_code=0, **kwargs):
        """
        Test a specific arg set, printing output and
        """
        print(name, "... ", sep="", end="", flush=True)
        out, time = run_tool(**kwargs)

        if out.returncode != expected_code:
            print("failed in ", round(time, 3), "s", sep="")
            print(out.stderr, file=sys.stderr)
            results.append(False)
            if critical:
                print("Critical failure, exiting")
                sys.exit(1)
        else:
            print("passed in ", round(time, 3), "s", sep="")
            results.append(True)

        return out.returncode != expected_code
    return f

def main():
    """Main"""
    results = []
    run_test = create_test_runner(results)

    root = "data/tests"
    os.makedirs(root, exist_ok=True)

    run_test("Basic functionality", f_file="", additional_args=["--help"])

    print("Generating test data")
    if not os.path.exists(f"{root}/perfect_grna.fq"):
        print("    Perfect gRNA library... ", end="")
        generate_test_data(lib_spec="config/grna.json", number=1000, library_size=100,
                           output=f"{root}/perfect_grna")
        print("done")

    if not os.path.exists(f"{root}/mutant_grna_sensor.fq"):
        print("    Mutant sensor gRNA library... ", end="")
        generate_test_data(lib_spec="config/grna_sensor.json", number=1000, library_size=100,
                           output=f"{root}/mutant_grna_sensor", contamination_rate=0,
                           recombination_rate=0, mismatch_rate=0,
                           sub_rate=0.01, indel_rate=0)
        print("done")

    if not os.path.exists(f"{root}/indel_grna_sensor.fq"):
        print("    Indel sensor gRNA library... ", end="")
        generate_test_data(lib_spec="config/grna_sensor.json", number=1000, library_size=100,
                           output=f"{root}/indel_grna_sensor", contamination_rate=0,
                           recombination_rate=0, mismatch_rate=0,
                           sub_rate=0, indel_rate=0.01)
        print("done")

    if not os.path.exists(f"{root}/recombined_grna_sensor.fq"):
        print("    Recombined sensor gRNA library... ", end="")
        generate_test_data(lib_spec="config/grna_sensor.json", number=1000, library_size=100,
                           output=f"{root}/recombined_grna_sensor", contamination_rate=0,
                           recombination_rate=0.05, mismatch_rate=0,
                           sub_rate=0, indel_rate=0)
        print("done")

    if not os.path.exists(f"{root}/mutant_pegrna.fq"):
        print("    Mutant pegRNA library... ", end="")
        generate_test_data(lib_spec="config/pegrna.json", number=1000, library_size=100,
                           output=f"{root}/mutant_pegrna", contamination_rate=0,
                           recombination_rate=0, mismatch_rate=0,
                           sub_rate=0.01, indel_rate=0)
        print("done")

    if not os.path.exists(f"{root}/messy_pegrna.fq"):
        print("    Fully mutated pegRNA library... ", end="")
        generate_test_data(lib_spec="config/pegrna.json", number=1000, library_size=100,
                           output=f"{root}/messy_pegrna", contamination_rate=0.05,
                           recombination_rate=0.05, mismatch_rate=0.05,
                           sub_rate=0.01, indel_rate=0.001)
        print("done")

    print("\nRunning main tests:")
    # Basic function
    run_test("Basic function", f_file=f"{root}/perfect_grna.fq", output="basic", mode="full-read",
             verbose=True, sort=True, library_counts=False, rm_output=True, critical=True)

    run_test("With IDs", f_file=f"{root}/perfect_grna.fq", lib_spec="config/grna.json",
             output="ids", mode="inframe", metric="hamming", critical=True,
             verbose=True, sort=True, rm_output=True)

    run_test("Without IDs", f_file=f"{root}/perfect_grna.fq", lib_spec="config/grna_no_id.json",
             output="no_ids", mode="inframe", metric="hamming", critical=True,
             verbose=True, sort=True, rm_output=True)

    # Full read counts
    run_test("Total counts", f_file=f"{root}/perfect_grna.fq",
             output=f"{root}/total_counts", mode="full-read",
             verbose=True, sort=True, library_counts=False, rm_output=False)
    run_test("Paired total counts", f_file=f"{root}/perfect_grna_forward.fq",
             r_file=f"{root}/perfect_grna_reverse.fq", output=f"{root}/total_pairs",
             mode="full-read", verbose=True, sort=True, library_counts=False,
             rm_output=False)

    # Combination of different params
    param_combs = product(
        ["inframe", "pattern", "align"],
        ["exact", "hamming", "bounded-levenshtein", "levenshtein"],
        ["perfect_grna", "mutant_grna_sensor", "indel_grna_sensor", "indel_grna_sensor",
         "recombined_grna_sensor", "mutant_pegrna", "messy_pegrna"]
    )
    lib_specs = {
        "perfect_grna": "grna",
        "mutant_grna_sensor": "grna_sensor",
        "indel_grna_sensor": "grna_sensor",
        "indel_grna_sensor": "grna_sensor",
         "recombined_grna_sensor": "grna_sensor",
         "mutant_pegrna": "pegrna",
         "messy_pegrna": "pegrna"
    }
    for (mode, dist, lib) in param_combs:
        expected_code = 1 if mode == "inframe" and lib in ["mutant_pegrna", "messy_pegrna"] else 0
        spec = f"config/{lib_specs[lib]}.json"

        run_test(f"{mode} {dist} {lib} single end", f_file=f"{root}/{lib}.fq",
                 output=f"{root}/params:{mode}:{dist}:{lib}:single:1",
                 mode=mode, metric=dist, verbose=True, sort=True,
                 lib_spec=spec,
                 expected_code=expected_code)

        run_test(f"{mode} {dist} {lib} paired end", f_file=f"{root}/{lib}_forward.fq",
                 r_file=f"{root}/{lib}_reverse.fq",
                 output=f"{root}/params:{mode}:{dist}:{lib}:paired:1",
                 mode=mode, metric=dist, verbose=True, sort=True,
                 lib_spec=spec,
                 expected_code=expected_code)

    for threads in [2, 4, 6]:
        run_test(f"{threads} Threads", f_file=f"{root}/mutant_pegrna.fq",
                 output=f"{root}/threads:align:bounded-levenshtein:mutant_pegrna:single:{threads}",
                 mode="align", metric="bounded-levenshtein", verbose=True, sort=True,
                 lib_spec="config/pegrna.json",
                 expected_code=0, additional_args=["--threads", str(threads)])

    # Filtering
    filter_args = ["--minimum-read-length", "80", "--maximum-read-length", "300",
                   "--alignment-tolerance", "0.9", "--mean-quality-threshold", "38"]
    run_test("Filtering", f_file=f"{root}/mutant_pegrna.fq",
             output=f"{root}/filtering", mode="align", metric="bounded-levenshtein",
             verbose=True, sort=True, lib_spec="config/pegrna.json",
             expected_code=0, additional_args=filter_args)

    print(f"\nTesting complete {sum(results)}/{len(results)} passed")
    if all(results):
        print(f"All tests passed")
        sys.exit(0)
    else:
        print(f"{sum(not i for i in results)} failures")
        sys.exit(1)

if __name__ == "__main__":
    main()