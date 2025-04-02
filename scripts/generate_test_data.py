#!/usr/bin/env python3
"""
Generate mock sequence data from an input LibrarySpec. Uses a simple approach that
doesn't capture all the possible error modes and doesn't account for multi-matching
when adding mutations (however, this should generally be rare and handled by other tests).
"""
import argparse
import sys
import json
import numpy as np
from collections import OrderedDict, Counter
from dataclasses import dataclass
from Bio.Seq import Seq

GENERATOR = np.random.Generator(np.random.PCG64())

MUTS = {
    "A": ["C", "G", "T"],
    "C": ["A", "G", "T"],
    "G": ["A", "C", "T"],
    "T": ["A", "C", "G"]
}

@dataclass(eq=True, frozen=True)
class GeneratedCombination:
    """
    Store a generated combination. Don't account for
    multi-matching and library variance - just store the
    expected sequence and the mutant, and similarly the
    intended library match.
    """
    group: str
    true_seqs: tuple
    observed_seqs: tuple
    combination_status: str
    combinations_in_library: str
    combination_index: str

def mutate_seq(s, sub_rate=0, indel_rate=0):
    """
    Add random indels and subs to a sequence
    """
    if sub_rate == 0 and indel_rate == 0:
        return s

    new = []
    i = 0
    while i < len(s):
        p = GENERATOR.random()
        if p < sub_rate:
            new.extend(GENERATOR.choice(MUTS[s[i]], size=1))
            i += 1

        # Otherwise indel chance
        elif p - sub_rate < indel_rate:
            # Decreasing chance of longer indel runs
            l = np.clip(GENERATOR.negative_binomial(1, 0.75), 1, 10)
            # Deletion
            if GENERATOR.random() < 0.5:
                i += l

            # Insertion
            else:
                new.extend(GENERATOR.choice(["A", "C", "G", "T"], size=l, replace=True))
                i += 1

        # No mutation
        else:
            new.append(s[i])
            i += 1

    return "".join(new)

def write_seq(name, seq, quality, filetype, file):
    """
    Write a fastq or fasta record
    """
    if filetype == "fa":
        print(">", name, "\n", seq, sep="", file=file)
    elif filetype == "fq":
        print("@", name, "\n", seq, "\n", "+", "\n", quality, sep="", file=file)
    else:
        raise ValueError("filetype must be fa or fq")

def sample_quality(n):
    """
    Sample quality scores using a rough Poisson methon
    """
    return 42 - np.clip(GENERATOR.poisson(5, size=n), 0, 42)

def generate_test_data(lib_spec, number=100, library_size=100, output="test_seqs", seqformat="fq",
                       ngroup=None, contamination_rate=0, recombination_rate=0, mismatch_rate=0,
                       sub_rate=0, indel_rate=0, truncation_rate=0):
    """
    Generate test data and write to files
    """
    # Load LibSpec
    with open(lib_spec, "r") as lib_spec_file:
        lib_spec = json.load(lib_spec_file)

    regions = {i["id"]: i for i in lib_spec["regions"]}
    variable_regions = [i["id"] for i in lib_spec["regions"] if i["seq_type"] in ("Library")]

    # Load library
    library = {}
    try:
        with open(lib_spec["library"], "r") as lib_file:
            l_regs = next(lib_file).strip().split("\t")

            for r in l_regs:
                library[r] = []

            for line in lib_file:
                line = line.strip().split("\t")
                for r, s in zip(l_regs, line):
                    library[r].append(s)

        lib_size = len(library[l_regs[0]])
    except:
        print("No library file, generating random library", file=sys.stderr)
        lib_size = library_size

    # Generate any missing variable regions
    for r in [i for i in variable_regions if not i in library.keys()]:
        library[r] = []
        for _ in range(lib_size):
            l = GENERATOR.integers(
                regions[r]["min_length"], regions[r]["max_length"], endpoint=True
            )
            library[r].append(''.join(GENERATOR.choice(["A", "C", "G", "T"], size=l, replace=True)))

    # Generate non-library mismatch elements
    non_library = {}
    for r in variable_regions:
        # Choose 10 random options for region mismatches
        non_library[r] = []
        attempts = 0
        while len(non_library[r]) < 10 and attempts < 100:
            attempts += 1
            l = GENERATOR.integers(
                regions[r]["min_length"], regions[r]["max_length"], endpoint=True
            )
            s = ''.join(GENERATOR.choice(["A", "C", "G", "T"], size=l, replace=True))
            if s not in library[r] and s not in non_library[r]:
                non_library[r].append(s)

    # Generate contaminants
    contaminants = []
    for _ in range(5):
        contaminants.append("".join(GENERATOR.choice(
            ["A", "C", "G", "T"], size=GENERATOR.integers(20, 120, size=1), replace=True
        )))

    # Generate sequences/counts
    counts = Counter()
    with (open(f"{output}.{seqformat}", "w") as mol_file,
          open(f"{output}_forward.{seqformat}", "w") as f_file,
          open(f"{output}_reverse.{seqformat}", "w") as r_file):
        for i in range(number):
            # Group
            group = f" group{GENERATOR.integers(0, ngroup)}" if ngroup else ""

            # Contamination
            if GENERATOR.random() < contamination_rate:
                name = f"seq{i}{group}"

                seq = contaminants[GENERATOR.integers(0, len(contaminants))]
                f_seq = seq[:lib_spec["forward_read_length"]]
                r_seq = str(Seq(seq).reverse_complement())[:lib_spec["reverse_read_length"]]

                # Quality
                quality = "".join(chr(i+33) for i in sample_quality(len(seq)))
                f_quality = "".join(chr(i+33) for i in sample_quality(len(f_seq)))
                r_quality = "".join(chr(i+33) for i in sample_quality(len(r_seq)))

                # Write seq
                write_seq(name, seq, quality, seqformat, mol_file)
                write_seq(name, f_seq, f_quality, seqformat, f_file)
                write_seq(name, r_seq, r_quality, seqformat, r_file)

                # Save match
                comb = GeneratedCombination(
                    group,
                    tuple("" for _ in range(len(variable_regions))),
                    tuple("" for _ in range(len(variable_regions))),
                    "nonmatch",
                    "0",
                    ""
                )

                counts[comb] += 1
                continue

            # Choose library member
            ind = GENERATOR.integers(0, lib_size)

            true_seqs = {r: library[r][ind] for r in variable_regions}
            match_type = "match"

            # Recombination
            if GENERATOR.random() < recombination_rate:
                recomb_ind = GENERATOR.integers(1, len(variable_regions))

                other_ind = ind
                while other_ind == ind:
                    other_ind = GENERATOR.integers(0, lib_size)

                for r in variable_regions[recomb_ind:]:
                    true_seqs[r] = library[r][other_ind]
                match_type = "recombination"

            observed_seqs = OrderedDict(
                (k, true_seqs[k] if k in true_seqs else v["seq"]) for k,v in regions.items()
            )

            if GENERATOR.random() < mismatch_rate:
                match_type = "mismatch"
                r = GENERATOR.choice(variable_regions)
                observed_seqs[r] = non_library[r][GENERATOR.integers(0, len(non_library[r]))]

            # Mutation
            if (sub_rate > 0 or indel_rate > 0):
                for k, v in observed_seqs.items():
                    observed_seqs[k] = mutate_seq(v, sub_rate, indel_rate)

            # Truncation
            if GENERATOR.random() < truncation_rate:
                trunc_ind = GENERATOR.integers(0, len(observed_seqs))
                keys = list(observed_seqs.keys())

                # Delete from random point in chosen region
                trunc_region = observed_seqs[keys[trunc_ind]]
                if len(trunc_region) > 2:
                    observed_seqs[keys[trunc_ind]] = trunc_region[:GENERATOR.integers(2, len(trunc_region))]

                for k in list(observed_seqs.keys())[trunc_ind + 1:]:
                    observed_seqs[k] = ""

                match_type = "nonmatch"

            name = f"seq{i}{group}"

            seq = "".join(observed_seqs[r] for r in regions.keys())
            f_seq = seq[:lib_spec["forward_read_length"]]
            r_seq = str(Seq(seq).reverse_complement())[:lib_spec["reverse_read_length"]]

            # Quality
            quality = "".join(chr(i+33) for i in sample_quality(len(seq)))
            f_quality = "".join(chr(i+33) for i in sample_quality(len(f_seq)))
            r_quality = "".join(chr(i+33) for i in sample_quality(len(r_seq)))

            # Write seq
            write_seq(name, seq, quality, seqformat, mol_file)
            write_seq(name, f_seq, f_quality, seqformat, f_file)
            write_seq(name, r_seq, r_quality, seqformat, r_file)

            # Save match
            comb = GeneratedCombination(
                    group,
                    tuple(true_seqs[r] for r in variable_regions),
                    tuple(observed_seqs[r] for r in variable_regions),
                    match_type,
                    1 if match_type == "match" else 0,
                    ind if match_type == "match" else ""
                )

            counts[comb] += 1

    # Write counts
    with open(f"{output}.true_counts.tsv", "w") as out_file:
        print("group", *variable_regions, *[f"{i}_nearest" for i in variable_regions],
              "combination_status", "combinations_in_library", "combination_indexes",
              "count", sep="\t", file=out_file)
        for comb, count in counts.items():
            print(comb.group, *comb.observed_seqs, *comb.true_seqs, comb.combination_status,
                comb.combinations_in_library, comb.combination_index,
                count, sep="\t", file=out_file)

def generate_library(lib_spec, n, path=None):
    """
    Generate a random test library, optionally writing it to a TSV and a new
    JSON pointing to the TSV
    """
    # Load LibSpec
    with open(lib_spec, "r") as lib_spec_file:
        lib_spec = json.load(lib_spec_file)

    regions = {i["id"]: i for i in lib_spec["regions"]}
    variable_regions = [i["id"] for i in lib_spec["regions"] if i["seq_type"] in ("Library")]

    library = {}
    for r in [i for i in variable_regions if not i in library.keys()]:
        library[r] = []
        for _ in range(n):
            l = GENERATOR.integers(
                regions[r]["min_length"], regions[r]["max_length"], endpoint=True
            )
            library[r].append(''.join(GENERATOR.choice(["A", "C", "G", "T"], size=l, replace=True)))

    if path is not None:
        with open(f"{path}.tsv", "w") as file:
            print(*variable_regions, sep="\t", file=file)
            for i in range(n):
                seqs = []
                for r in variable_regions:
                    seqs.append(library[r][i])
                print(*seqs, sep="\t", file=file)

        with open(f"{path}.json", "w") as file:
            lib_spec["library"] = f"{path}.tsv"
            json.dump(lib_spec, file)

    return library

def main():
    """
    Main
    """
    # from scripts.generate_test_data import *
    # args = parse_args(["benchmark/puffin.json"])
    args = parse_args()
    generate_test_data(**args)

def parse_args(arg_list=None):
    """
    Parse and validate script arguments
    """
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("lib_spec", metavar="J", help="LibSPec JSON file")

    parser.add_argument("--number", "-n", type=int, default=10,
                        help="Number of sequences to generate")

    parser.add_argument("--library", "-l", type=int, default=None,
                        help="Library size to generate if no library TSV")

    parser.add_argument("--output", "-o", type=str, default="test_seqs",
                        help="Root output path")

    parser.add_argument("--seqformat", "-f", type=str, default="fa", choices=("fa", "fq"),
                        help="Output sequence format")

    parser.add_argument("--ngroup", "-g", type=int, default=None,
                        help="Add random read groups formatted as >seqN groupG")

    parser.add_argument("--recombination_rate", "-r", type=float, default=0,
                        help="Recombination rate")

    parser.add_argument("--mismatch_rate", "-m", type=float, default=0,
                        help="Mismatch rate per region")

    parser.add_argument("--sub_rate", "-s", type=float, default=0,
                        help="Substitution rate")

    parser.add_argument("--indel_rate", "-i", type=float, default=0,
                        help="Indel rate")

    parser.add_argument("--truncation_rate", "-t", type=float, default=0,
                        help="Truncation rate")

    parser.add_argument("--contamination_rate", "-c", type=float, default=0,
                        help="Contamination rate")

    return parser.parse_args(arg_list)

if __name__ == "__main__":
    main()