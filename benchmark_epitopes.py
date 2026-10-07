#!/usr/bin/env python3
"""
Epitope-level benchmark: documented neoantigens on their own restricting allele
===============================================================================

WHAT THIS IS FOR

This is the smallest honest test of the binding predictor. It takes the
peptides in external_validation.KNOWN_IMMUNOGENIC_EPITOPES, scores each one
against the HLA allele it is actually documented as restricted to, and reports
how many are recovered.

It exists because the previous benchmark in 06_validate_pipeline.py scored
every positive control against a fixed six-allele panel that did not contain
the restricting alleles for several of them. KRAS G12D is restricted to
HLA-A*11:01 and HLA-C*08:02; neither was in the panel, so the control was
scored against molecules that cannot present it and then recorded as a miss.

WHY PROVENANCE IS ENFORCED HERE

Only epitopes with a real published source count towards the headline number.
Two entries in the list are computational predictions with no experimental
measurement behind them. Scoring a predictor against another predictor's output
measures agreement, not accuracy, and those two happen to produce the strongest
affinities in the table. Counting them would inflate the result and prove
nothing. They are printed, clearly labelled, and excluded from the total.

REQUIREMENTS

MHCflurry must be installed and its models downloaded:

    pip install mhcflurry
    mhcflurry-downloads fetch models_class1_presentation

Without MHCflurry there is nothing meaningful to run. The built-in PSSM in
03_predict_binding.py defines matrices for four alleles only (A*02:01, A*01:01,
A*03:01, B*07:02) and cannot score HLA-A*11:01 or HLA-C*08:02 at all, so a
PSSM-only run of this benchmark would report misses that are artefacts of the
fallback rather than properties of the prediction.

Usage:
    python benchmark_epitopes.py
"""

import argparse
import sys

import pandas as pd

from external_validation import (
    KNOWN_IMMUNOGENIC_EPITOPES,
    EXPERIMENTAL_PROVENANCE,
)

# Standard interpretation thresholds for predicted affinity, in nM.
STRONG_BINDER_NM = 50
BINDER_NM = 500


def load_predictor():
    """Load the MHCflurry presentation predictor, or exit with instructions."""
    try:
        from mhcflurry import Class1PresentationPredictor
    except ImportError:
        sys.exit(
            "MHCflurry is not installed, and this benchmark is not meaningful "
            "without it.\n\n"
            "    pip install mhcflurry\n"
            "    mhcflurry-downloads fetch models_class1_presentation\n"
        )

    try:
        return Class1PresentationPredictor.load()
    except Exception as exc:  # models not fetched, or a version mismatch
        sys.exit(
            f"MHCflurry is installed but its models could not be loaded: {exc}\n\n"
            "    mhcflurry-downloads fetch models_class1_presentation\n"
        )


def score_epitopes(predictor, threshold_nm=BINDER_NM):
    """Score each documented epitope against its own restricting allele."""
    rows = []
    for entry in KNOWN_IMMUNOGENIC_EPITOPES:
        peptide = entry["peptide"]
        allele = entry["hla_restriction"]
        provenance = entry.get("provenance", "none")

        try:
            result = predictor.predict(
                peptides=[peptide], alleles=[allele], verbose=0
            ).iloc[0]
            affinity = float(result["affinity"])
            presentation_pct = float(result["presentation_percentile"])
            error = None
        except Exception as exc:
            affinity = None
            presentation_pct = None
            error = str(exc)[:60]

        rows.append({
            "gene": entry["gene"],
            "mutation": entry["mutation"],
            "peptide": peptide,
            "allele": allele,
            "provenance": provenance,
            "counts_towards_total": provenance in EXPERIMENTAL_PROVENANCE,
            "affinity_nM": round(affinity, 1) if affinity is not None else None,
            "presentation_pct": (
                round(presentation_pct, 3) if presentation_pct is not None else None
            ),
            "recovered": (
                affinity is not None and affinity < threshold_nm
            ),
            "strong": (
                affinity is not None and affinity < STRONG_BINDER_NM
            ),
            "error": error,
        })

    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--threshold", type=float, default=BINDER_NM,
        help=f"Affinity threshold in nM for calling a binder (default {BINDER_NM})",
    )
    parser.add_argument(
        "--out", default="output/validation/epitope_benchmark.csv",
        help="Where to write the full results table",
    )
    args = parser.parse_args()

    predictor = load_predictor()
    results = score_epitopes(predictor, threshold_nm=args.threshold)

    pd.set_option("display.width", 240)
    display_cols = [
        "gene", "mutation", "peptide", "allele", "provenance",
        "affinity_nM", "presentation_pct", "recovered",
    ]
    print()
    print(results[display_cols].to_string(index=False))
    print()

    counted = results[results["counts_towards_total"]]
    excluded = results[~results["counts_towards_total"]]

    recovered = int(counted["recovered"].sum())
    total = len(counted)

    print(f"Recovered {recovered}/{total} experimentally-sourced epitopes "
          f"at <{args.threshold:.0f} nM.")
    print()
    print("Excluded from that total, with reasons:")
    for _, row in excluded.iterrows():
        print(f"  {row['gene']} {row['mutation']} on {row['allele']}: "
              f"{row['provenance']} "
              f"({row['affinity_nM']} nM)")
    print()
    print("Read the misses before quoting the number. A miss can mean the")
    print("predictor is weak on that allele rather than that the epitope is")
    print("absent; HLA-C is substantially less well covered by affinity")
    print("predictors than HLA-A and HLA-B, because far less training data")
    print("exists for it. With four scoreable epitopes, this benchmark is a")
    print("smoke test, not an accuracy measurement.")

    import os
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    results.to_csv(args.out, index=False)
    print()
    print(f"Full table written to {args.out}")


if __name__ == "__main__":
    main()
