#!/usr/bin/env python3
"""
Generate accuracy-profile datasets by running DNAComb on
simulated sequences with varying levels of perturbation.

Produces:
  - simulated FASTQ/true-count files from generate_test_data.py
  - DNAComb output files
  - a manifest TSV describing every run
"""
import argparse
import os
import sys
from dataclasses import dataclass
from itertools import product

from generate_test_data import generate_test_data
from utils import run_tool

@dataclass(frozen=True)
class PerturbationLevel:
    level: int
    sub_rate: float
    indel_rate: float
    recombination_rate: float
    contamination_rate: float
    mismatch_rate: float

@dataclass(frozen=True)
class Design:
    name: str
    lib_spec: str
    libraries: list[str]
    library_size: int
    allow_inframe: bool

PERTURBATION_LEVELS = [
    PerturbationLevel(0, 0.000, 0.000, 0.00, 0.00, 0.00),
    PerturbationLevel(1, 0.005, 0.000, 0.00, 0.00, 0.01),
    PerturbationLevel(2, 0.010, 0.001, 0.00, 0.01, 0.01),
    PerturbationLevel(3, 0.010, 0.005, 0.01, 0.01, 0.01),
    PerturbationLevel(4, 0.010, 0.010, 0.02, 0.02, 0.02),
    PerturbationLevel(5, 0.050, 0.010, 0.05, 0.05, 0.05),
    PerturbationLevel(6, 0.050, 0.020, 0.10, 0.10, 0.05),
]

DESIGNS = [
    Design(
        name="grna_sensor",
        lib_spec="config/grna_sensor.json",
        libraries=["config/grna_sensor.tsv"],
        library_size=100,
        allow_inframe=True,
    ),
    Design(
        name="pegrna",
        lib_spec="config/pegrna.json",
        libraries=["config/pegrna.tsv"],
        library_size=100,
        allow_inframe=False,
    ),
]

MODES = ["align", "pattern", "inframe"]
METRICS = ["exact", "hamming", "bounded-levenshtein"]
ENDS = ["single", "paired"]

def generate_dataset(
    design: Design,
    level: PerturbationLevel,
    replicate: int,
    n_reads: int,
    root: str,
) -> dict[str, str]:
    """
    Generate a dataset for a given sequence design and perturbation level
    """
    stem = os.path.join(
        root,
        "inputs",
        f"{design.name}.level{level.level}.rep{replicate}",
    )

    os.makedirs(os.path.dirname(stem), exist_ok=True)

    print(
        f"Generating {design.name} level {level.level} replicate {replicate}...",
        flush=True,
    )

    generate_test_data(
        lib_spec=design.lib_spec,
        libraries=design.libraries,
        number=n_reads,
        library_size=design.library_size,
        output=stem,
        contamination_rate=level.contamination_rate,
        recombination_rate=level.recombination_rate,
        mismatch_rate=level.mismatch_rate,
        sub_rate=level.sub_rate,
        indel_rate=level.indel_rate,
    )

    return {
        "stem": stem,
        "single_fq": f"{stem}.fq",
        "paired_fq_f": f"{stem}_forward.fq",
        "paired_fq_r": f"{stem}_reverse.fq",
        "true_counts": f"{stem}.true_counts.tsv",
    }


def write_manifest_header(handle):
    print(
        "design",
        "level",
        "replicate",
        "sub_rate",
        "indel_rate",
        "recombination_rate",
        "contamination_rate",
        "mismatch_rate",
        "mode",
        "metric",
        "end",
        "lib_spec",
        "libraries",
        "input_f",
        "input_r",
        "true_counts",
        "output_prefix",
        "counts_tsv",
        "library_counts_tsv",
        "summary_tsv",
        "filtered_tsv",
        "returncode",
        sep="\t",
        file=handle,
    )


def write_manifest_row(
    handle,
    design: Design,
    level: PerturbationLevel,
    replicate: int,
    mode: str,
    metric: str,
    end: str,
    input_f: str,
    input_r: str,
    true_counts: str,
    output_prefix: str,
    returncode: int,
):
    """
    Write a row to the manifest table describing a run
    """

    print(
        design.name,
        level.level,
        replicate,
        level.sub_rate,
        level.indel_rate,
        level.recombination_rate,
        level.contamination_rate,
        level.mismatch_rate,
        mode,
        metric,
        end,
        design.lib_spec,
        ",".join(design.libraries),
        input_f,
        input_r,
        true_counts,
        output_prefix,
        f"{output_prefix}.counts.tsv",
        f"{output_prefix}.library_counts.tsv",
        f"{output_prefix}.summary.tsv",
        f"{output_prefix}.filtered.tsv",
        returncode,
        sep="\t",
        file=handle,
        flush=True,
    )

def main():
    """
    Generate test data and run DNAComb on it over a range of designs
    and accuracy levels with replicates.
    """
    root = "data/accuracy_profile"
    reads = 10000
    replicates = 5

    os.makedirs(root, exist_ok=True)
    os.makedirs(os.path.join(root, "outputs"), exist_ok=True)

    manifest_path = os.path.join(root, "accuracy_profile_manifest.tsv")

    with open(manifest_path, "w") as manifest:
        write_manifest_header(manifest)

        for design, replicate, level in product(
            DESIGNS,
            range(1, replicates + 1),
            PERTURBATION_LEVELS,
        ):
            dataset = generate_dataset(
                design=design,
                level=level,
                replicate=replicate,
                n_reads=reads,
                root=root,
            )

            for mode, metric, end in product(MODES, METRICS, ENDS):
                if mode == "inframe" and not design.allow_inframe:
                    continue

                input_f = dataset["single_fq"]
                input_r = None

                if end == "paired":
                    input_f = dataset["paired_fq_f"]
                    input_r = dataset["paired_fq_r"]

                run_name = (
                    f"{design.name}"
                    f".level{level.level}"
                    f".rep{replicate}"
                    f".{end}"
                    f".{mode}"
                    f".{metric}"
                )

                output_prefix = os.path.join(root, "outputs", run_name)

                print(f"Running {run_name}... ", end="", flush=True)

                result, time = run_tool(
                    f_file=input_f,
                    r_file=input_r,
                    lib_spec=design.lib_spec,
                    output=output_prefix,
                    mode=mode,
                    metric=metric,
                    verbose=True,
                    sort=True,
                    overwrite=True,
                    library=design.libraries,
                )

                returncode = result.returncode

                if result.returncode != 0:
                    print("failed in ", round(time, 3), "s", sep="")
                    print(result.stderr, file=sys.stderr)
                    sys.exit(1)
                else:
                    print("completed in ", round(time, 3), "s", sep="")

                write_manifest_row(
                    handle=manifest,
                    design=design,
                    level=level,
                    replicate=replicate,
                    mode=mode,
                    metric=metric,
                    end=end,
                    input_f=input_f,
                    input_r=input_r or "",
                    true_counts=dataset["true_counts"],
                    output_prefix=output_prefix,
                    returncode=returncode,
                )

    print(f"\nManifest written to: {manifest_path}")

if __name__ == "__main__":
    main()