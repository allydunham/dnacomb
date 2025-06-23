#!/usr/bin/env python3
"""
Run end-to-end tests using a variety of input conditions
"""
import sys
import os
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

    if not os.path.isdir("tests"):
        os.mkdir("tests")

    run_test("Basic functionality", f_file="", additional_args=["--help"])

    print("Generating test data")
    if not os.path.exists("tests/grna.fq"):
        print("    Basic gRNA library... ", end="")
        generate_test_data(lib_spec="config/grna.json", number=10000, library_size=100,
                           output="tests/grna")
        print("done")

    if not os.path.exists("tests/grna.fa"):
        print("    Basic gRNA Fasta library... ", end="")
        generate_test_data(lib_spec="config/grna.json", number=10000, library_size=100,
                           output="tests/grna", seqformat="fa")
        print("done")

    if not os.path.exists("tests/grna_sensor.fq"):
        print("    Mutant sensor gRNA library... ", end="")
        generate_test_data(lib_spec="config/grna_sensor.json", number=10000, library_size=100,
                           output="tests/grna_sensor", contamination_rate=0.05,
                           recombination_rate=0.05, mismatch_rate=0.05,
                           sub_rate=0.01, indel_rate=0.001)
        print("done")

    if not os.path.exists("tests/pegrna.fq"):
        print("    pegRNA library... ", end="")
        generate_test_data(lib_spec="config/pegrna.json", number=10000, library_size=100,
                           output="tests/pegrna", contamination_rate=0.05,
                           recombination_rate=0.05, mismatch_rate=0.05,
                           sub_rate=0.01, indel_rate=0.001)
        print("done")

    print("\nRunning main tests:")
    # Reading plain fasta
    run_test("Fasta reading", f_file="tests/grna.fa", output="fasta", mode="full-read",
             verbose=True, sort=True, library_counts=False, rm_output=True)
    run_test("Fastq reading", f_file="tests/grna.fq", output="fastq", mode="full-read",
             verbose=True, sort=True, library_counts=False, rm_output=True)

    # Full read counts
    run_test("Total counts", f_file="tests/grna.fq", output="total_counts", mode="full-read",
             verbose=True, sort=True, library_counts=False, rm_output=True)
    run_test("Paired total counts", f_file="tests/grna_forward.fq", r_file="tests/grna_reverse.fq",
             output="total_pairs", mode="full-read", verbose=True, sort=True, library_counts=False,
             rm_output=True)

    # Exact matches using each mode in single/paired
    for i in ["inframe", "pattern", "align"]:
        run_test(f"Exact {i}", f_file="tests/grna_sensor.fq", output=f"tests/exact_sensor_{i}",
                 lib_spec="config/grna_sensor.json", mode=i, metric="exact", verbose=True,
                 sort=True)
        run_test(f"Exact paired {i}", f_file="tests/grna_sensor_forward.fq",
                r_file="tests/grna_sensor_reverse.fq", output=f"tests/exact_sensor_{i}_paired",
                lib_spec="config/grna_sensor.json", mode=i, metric="exact", verbose=True,
                sort=True)

    # Variable length regions
    for i in ["inframe", "pattern", "align"]:
        run_test(f"Variable exact {i}", f_file="tests/pegrna.fq",
                 output=f"tests/pegrna_exact_{i}",
                 mode=i, metric="exact", verbose=True, sort=True,
                 lib_spec="config/pegrna.json",
                 expected_code=1 if i == "inframe" else 0)
        run_test(f"Variable exact paired {i}", f_file="tests/pegrna_forward.fq",
                 r_file="tests/pegrna_reverse.fq", output=f"tests/pegrna_exact_{i}_paired",
                 mode=i, metric="exact", verbose=True, sort=True,
                 lib_spec="config/pegrna.json",
                 expected_code=1 if i == "inframe" else 0)

    # Distance metrics
    for i in ["hamming", "bounded-levenshtein", "levenshtein"]:
        run_test(f"{i}", f_file="tests/grna_sensor.fq", output=f"tests/{i}_sensor_inframe",
                 mode="inframe", metric=i, verbose=True, sort=True,
                 lib_spec="config/grna_sensor.json")

    print(f"\nTesting complete {sum(results)}/{len(results)} passed")
    if all(results):
        print(f"All tests passed")
        sys.exit(0)
    else:
        print(f"{sum(not i for i in results)} failures")
        sys.exit(1)

if __name__ == "__main__":
    main()