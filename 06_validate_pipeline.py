#!/usr/bin/env python3
"""
Pipeline Validation & End-to-End Simulation
=============================================
Runs an offline end-to-end simulation of the neoantigen vaccine pipeline
using well-characterized benchmark mutations from Multiple Myeloma.

This bypasses the GDC API fetch step and uses curated test data to:
1. Validate each pipeline step produces correct output
2. Benchmark predictions against known immunogenic mutations
3. Cross-reference results against COSMIC/IEDB databases
4. Generate a validation report

Usage:
    python 06_validate_pipeline.py [--config config.yaml]
"""

import os
import sys
import json
import time
from pathlib import Path
from io import StringIO

import pandas as pd
import numpy as np
import yaml

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from importlib import import_module

# Import pipeline modules
parse_mod = import_module("02_parse_mutations")
binding_mod = import_module("03_predict_binding")
design_mod = import_module("04_design_vaccine")


# =============================================================================
# Benchmark Mutation Dataset
# =============================================================================
# Well-characterized mutations from Multiple Myeloma and pan-cancer studies
# with known immunogenic properties. These are used to validate the pipeline.
# =============================================================================

BENCHMARK_MUTATIONS = [
    # --- KRAS mutations (most common oncogene in MM, ~23%) ---
    {
        "case_id": "BENCH-001", "submitter_id": "MMRF_BENCH_001",
        "gene_symbol": "KRAS", "gene_id": "ENSG00000133703",
        "chromosome": "chr12", "start_position": 25245350,
        "aa_change": "p.G12D", "consequence_type": "missense_variant",
        "ssm_id": "bench_kras_g12d",
        "reference_allele": "C", "tumor_allele": "T",
        "expected_immunogenic": True,
        "notes": "Most common KRAS mutation in MM, well-validated neoantigen",
    },
    {
        "case_id": "BENCH-001", "submitter_id": "MMRF_BENCH_001",
        "gene_symbol": "KRAS", "gene_id": "ENSG00000133703",
        "chromosome": "chr12", "start_position": 25245350,
        "aa_change": "p.G12V", "consequence_type": "missense_variant",
        "ssm_id": "bench_kras_g12v",
        "reference_allele": "C", "tumor_allele": "A",
        "expected_immunogenic": True,
        "notes": "Known HLA-A*02:01 restricted neoantigen",
    },
    # --- NRAS mutations (very common in MM, ~20%) ---
    {
        "case_id": "BENCH-002", "submitter_id": "MMRF_BENCH_002",
        "gene_symbol": "NRAS", "gene_id": "ENSG00000213281",
        "chromosome": "chr1", "start_position": 114713908,
        "aa_change": "p.Q61K", "consequence_type": "missense_variant",
        "ssm_id": "bench_nras_q61k",
        "reference_allele": "G", "tumor_allele": "T",
        "expected_immunogenic": True,
        "notes": "Common NRAS hotspot in MM, charge-changing mutation",
    },
    {
        "case_id": "BENCH-002", "submitter_id": "MMRF_BENCH_002",
        "gene_symbol": "NRAS", "gene_id": "ENSG00000213281",
        "chromosome": "chr1", "start_position": 114713908,
        "aa_change": "p.Q61R", "consequence_type": "missense_variant",
        "ssm_id": "bench_nras_q61r",
        "reference_allele": "T", "tumor_allele": "C",
        "expected_immunogenic": True,
        "notes": "Charge-changing mutation, highly immunogenic",
    },
    # --- BRAF V600E (found in ~4% of MM) ---
    {
        "case_id": "BENCH-003", "submitter_id": "MMRF_BENCH_003",
        "gene_symbol": "BRAF", "gene_id": "ENSG00000157764",
        "chromosome": "chr7", "start_position": 140753336,
        "aa_change": "p.V600E", "consequence_type": "missense_variant",
        "ssm_id": "bench_braf_v600e",
        "reference_allele": "A", "tumor_allele": "T",
        "expected_immunogenic": True,
        "notes": "Most studied cancer mutation, validated in melanoma vaccines",
    },
    # --- TP53 mutations (found in ~8% of MM, poor prognosis) ---
    {
        "case_id": "BENCH-004", "submitter_id": "MMRF_BENCH_004",
        "gene_symbol": "TP53", "gene_id": "ENSG00000141510",
        "chromosome": "chr17", "start_position": 7675088,
        "aa_change": "p.R175H", "consequence_type": "missense_variant",
        "ssm_id": "bench_tp53_r175h",
        "reference_allele": "C", "tumor_allele": "T",
        "expected_immunogenic": True,
        "notes": "Most common TP53 hotspot, charge-changing (R→H)",
    },
    {
        "case_id": "BENCH-004", "submitter_id": "MMRF_BENCH_004",
        "gene_symbol": "TP53", "gene_id": "ENSG00000141510",
        "chromosome": "chr17", "start_position": 7675229,
        "aa_change": "p.R248W", "consequence_type": "missense_variant",
        "ssm_id": "bench_tp53_r248w",
        "reference_allele": "G", "tumor_allele": "A",
        "expected_immunogenic": True,
        "notes": "TP53 hotspot, charge-changing (R→W)",
    },
    # --- DIS3 (MM-specific, ~11%) ---
    {
        "case_id": "BENCH-005", "submitter_id": "MMRF_BENCH_005",
        "gene_symbol": "DIS3", "gene_id": "ENSG00000083520",
        "chromosome": "chr13", "start_position": 73346251,
        "aa_change": "p.R780K", "consequence_type": "missense_variant",
        "ssm_id": "bench_dis3_r780k",
        "reference_allele": "C", "tumor_allele": "T",
        "expected_immunogenic": False,
        "notes": "Conservative mutation (R→K), both positive charge",
    },
    # --- FAM46C/TENT5C (MM-specific, ~11%) ---
    {
        "case_id": "BENCH-005", "submitter_id": "MMRF_BENCH_005",
        "gene_symbol": "FAM46C", "gene_id": "ENSG00000154277",
        "chromosome": "chr1", "start_position": 27105367,
        "aa_change": "p.D249N", "consequence_type": "missense_variant",
        "ssm_id": "bench_fam46c_d249n",
        "reference_allele": "C", "tumor_allele": "T",
        "expected_immunogenic": True,
        "notes": "Charge-changing (D→N), MM-specific driver",
    },
    # --- Passenger mutation (negative control) ---
    {
        "case_id": "BENCH-006", "submitter_id": "MMRF_BENCH_006",
        "gene_symbol": "TTN", "gene_id": "ENSG00000155657",
        "chromosome": "chr2", "start_position": 179390716,
        "aa_change": "p.A1234V", "consequence_type": "missense_variant",
        "ssm_id": "bench_ttn_a1234v",
        "reference_allele": "G", "tumor_allele": "A",
        "expected_immunogenic": False,
        "notes": "TTN is a very large gene with frequent passenger mutations",
    },
    {
        "case_id": "BENCH-006", "submitter_id": "MMRF_BENCH_006",
        "gene_symbol": "MUC16", "gene_id": "ENSG00000181143",
        "chromosome": "chr19", "start_position": 8959520,
        "aa_change": "p.S5000L", "consequence_type": "missense_variant",
        "ssm_id": "bench_muc16_s5000l",
        "reference_allele": "G", "tumor_allele": "A",
        "expected_immunogenic": False,
        "notes": "MUC16 (CA-125) frequent passenger, not a cancer driver",
    },
]


# =============================================================================
# Documented HLA restriction for the benchmark positives
# =============================================================================
#
# WHY THIS TABLE EXISTS
#
# Earlier revisions scored every benchmark positive against whichever allele in
# the default panel gave the best affinity. The default panel is HLA-A*02:01,
# A*01:01, A*03:01, A*24:02, B*07:02 and B*08:01. KRAS G12D is restricted to
# HLA-C*08:02 and HLA-A*11:01, neither of which is in that panel, so the
# pipeline was asked to find the epitope on molecules that cannot present it
# and then marked wrong for not finding it. Seven of eight positives failed
# that way, which made the validation run look like a catastrophic sensitivity
# problem when it was an allele bookkeeping problem.
#
# A positive control is only meaningful against its own restricting allele.
# Where no restriction is documented, the entry cannot be scored at all and is
# reported as INDETERMINATE rather than counted as a miss.
#
# Evidence grades mirror external_validation.EXPERIMENTAL_PROVENANCE:
#   published_pmid     citation with a PMID recorded (PMID not re-verified here)
#   published_no_pmid  author-and-year citation only
#   unsourced          no usable citation; not scoreable
#   computational      predicted only; cannot serve as a positive control for a
#                      prediction pipeline
#   none               no epitope record exists for this mutation
#
# Only published_pmid and published_no_pmid rows count towards sensitivity.

BENCHMARK_RESTRICTION = {
    "p.G12D": {
        "gene": "KRAS",
        "restricting_alleles": ["HLA-A*11:01", "HLA-C*08:02"],
        "evidence": "published_pmid",
        "source": "Wang et al. 2019 PMID:31537801; Tran et al. 2016 PMID:27959684",
    },
    "p.G12V": {
        "gene": "KRAS",
        "restricting_alleles": ["HLA-A*02:01", "HLA-A*11:01"],
        "evidence": "published_no_pmid",
        "source": "Cafri et al. 2019; A*11:01 record is unsourced",
    },
    "p.R175H": {
        "gene": "TP53",
        "restricting_alleles": ["HLA-A*02:01"],
        "evidence": "published_no_pmid",
        "source": "Lo et al. 2019",
    },
    "p.V600E": {
        "gene": "BRAF",
        "restricting_alleles": ["HLA-A*02:01"],
        "evidence": "unsourced",
        "source": "epitope record cited only 'Melanoma vaccine trials'",
    },
    "p.Q61K": {
        "gene": "NRAS",
        "restricting_alleles": ["HLA-A*01:01"],
        "evidence": "computational",
        "source": "predicted binding only, no experimental measurement",
    },
    "p.Q61R": {
        "gene": "NRAS",
        "restricting_alleles": ["HLA-A*01:01"],
        "evidence": "computational",
        "source": "predicted binding only, no experimental measurement",
    },
    "p.R248W": {
        "gene": "TP53",
        "restricting_alleles": [],
        "evidence": "none",
        "source": "no epitope record; restriction not established here",
    },
}

SCOREABLE_EVIDENCE = {"published_pmid", "published_no_pmid"}


def create_benchmark_dataframe():
    """Create a DataFrame from benchmark mutations."""
    # Strip extra keys not in the pipeline's expected columns
    records = []
    for m in BENCHMARK_MUTATIONS:
        record = {k: v for k, v in m.items()
                  if k not in ("expected_immunogenic", "notes")}
        records.append(record)
    return pd.DataFrame(records)


def run_validation(config_path="config.yaml"):
    """Run full pipeline validation with benchmark data."""
    config = yaml.safe_load(open(config_path))
    results = {
        "steps": {},
        "benchmark_results": [],
        "summary": {},
    }

    print("=" * 70)
    print("  PIPELINE VALIDATION & END-TO-END SIMULATION")
    print("  Using benchmark mutations from MM literature")
    print("=" * 70)

    # =========================================================================
    # STEP 1: Mutation Parsing (bypass GDC fetch)
    # =========================================================================
    print(f"\n{'─' * 70}")
    print("  STEP 1: Parsing benchmark mutations")
    print(f"{'─' * 70}")

    mutations_df = create_benchmark_dataframe()
    print(f"  Loaded {len(mutations_df)} benchmark mutations")
    print(f"  Patients: {mutations_df['case_id'].nunique()}")
    print(f"  Genes: {', '.join(mutations_df['gene_symbol'].unique())}")

    # Filter mutations
    filtered_df = parse_mod.filter_mutations(mutations_df, config)
    print(f"  After filtering: {len(filtered_df)} mutations")

    # Process into neoantigen candidates
    candidates_df = parse_mod.process_mutations(filtered_df, config)
    print(f"  Generated {len(candidates_df)} peptide candidates")

    results["steps"]["parsing"] = {
        "input_mutations": len(mutations_df),
        "filtered_mutations": len(filtered_df),
        "candidates_generated": len(candidates_df),
    }

    # =========================================================================
    # STEP 1b: UniProt protein sequence lookup (if available)
    # =========================================================================
    try:
        from uniprot_lookup import generate_real_peptides, KNOWN_MM_PROTEINS
        print(f"\n{'─' * 70}")
        print("  STEP 1b: UniProt protein sequence enrichment")
        print(f"{'─' * 70}")
        print(f"  Pre-cached proteins: {', '.join(KNOWN_MM_PROTEINS.keys())}")

        uniprot_count = 0
        for gene in mutations_df['gene_symbol'].unique():
            if gene in KNOWN_MM_PROTEINS:
                uniprot_count += 1
                seq = KNOWN_MM_PROTEINS[gene]
                print(f"  {gene}: {len(seq)} aa sequence available")

        print(f"  {uniprot_count}/{len(mutations_df['gene_symbol'].unique())} "
              f"genes have real protein sequences")
        results["steps"]["uniprot"] = {"genes_with_sequences": uniprot_count}
    except ImportError:
        print("\n  [Skip] uniprot_lookup module not available")

    # =========================================================================
    # STEP 2: MHC Binding Prediction
    # =========================================================================
    print(f"\n{'─' * 70}")
    print("  STEP 2: MHC binding prediction")
    print(f"{'─' * 70}")

    # Limit candidates for speed
    max_cands = min(len(candidates_df), 2000)
    if len(candidates_df) > max_cands:
        candidates_df = candidates_df.nlargest(max_cands, "immunogenicity_score")
        print(f"  Limited to top {max_cands} candidates by immunogenicity score")

    # The benchmark panel must contain every allele the positive controls are
    # restricted to, or those controls cannot be scored. The default config
    # panel omits HLA-A*11:01 and HLA-C*08:02, which between them carry the
    # KRAS G12D evidence, so they are added here for the validation run only.
    #
    # Caveat: the built-in PSSM in 03_predict_binding.py defines matrices for
    # four alleles only (A*02:01, A*01:01, A*03:01, B*07:02). Without MHCflurry
    # installed, the added alleles fall back to a generic approximation and the
    # resulting numbers are not interpretable. Install MHCflurry before reading
    # anything into this benchmark.
    required_alleles = sorted({
        allele
        for entry in BENCHMARK_RESTRICTION.values()
        for allele in entry.get("restricting_alleles", [])
    })
    panel = config["neoantigen"]["default_hla_alleles"]["class_i"]
    added = [a for a in required_alleles if a not in panel]
    if added:
        config["neoantigen"]["default_hla_alleles"]["class_i"] = panel + added
        print(f"  Added restricting alleles to benchmark panel: {', '.join(added)}")
    results["steps"]["benchmark_panel"] = {
        "configured_panel": panel,
        "added_for_benchmark": added,
        "pssm_supported_alleles": sorted(binding_mod.PSSM_MODELS.keys()),
    }

    binding_df = binding_mod.run_binding_predictions(candidates_df, config)
    print(f"  Completed predictions for {len(binding_df)} candidate-allele pairs")

    # Summarize
    binders = binding_df[binding_df["classification"].isin(["strong_binder", "weak_binder"])]
    strong = binding_df[binding_df["classification"] == "strong_binder"]
    print(f"  Binders: {len(binders)} ({100*len(binders)/max(len(binding_df),1):.1f}%)")
    print(f"  Strong binders: {len(strong)} ({100*len(strong)/max(len(binding_df),1):.1f}%)")

    results["steps"]["binding"] = {
        "predictions": len(binding_df),
        "binders": len(binders),
        "strong_binders": len(strong),
    }

    # Rank
    ranked_df = binding_mod.rank_neoantigens(binding_df, config)

    # =========================================================================
    # STEP 2b: External validation (if available)
    # =========================================================================
    try:
        from external_validation import annotate_candidates
        print(f"\n{'─' * 70}")
        print("  STEP 2b: COSMIC/IEDB cross-referencing")
        print(f"{'─' * 70}")

        annotated_df = annotate_candidates(ranked_df)
        hotspots = annotated_df[annotated_df.get("cosmic_status", pd.Series()) == "hotspot"]
        iedb_validated = annotated_df[annotated_df.get("iedb_validated", pd.Series()) == True]
        print(f"  COSMIC hotspot mutations: {len(hotspots)}")
        print(f"  IEDB-validated epitopes: {len(iedb_validated)}")
        ranked_df = annotated_df
        results["steps"]["external_validation"] = {
            "cosmic_hotspots": len(hotspots) if len(hotspots) > 0 else 0,
            "iedb_validated": len(iedb_validated) if len(iedb_validated) > 0 else 0,
        }
    except (ImportError, Exception) as e:
        print(f"\n  [Skip] external_validation: {e}")

    # =========================================================================
    # STEP 2c: Expression filtering (if available)
    # =========================================================================
    try:
        from expression_filter import enrich_with_expression
        print(f"\n{'─' * 70}")
        print("  STEP 2c: Gene expression enrichment")
        print(f"{'─' * 70}")

        enriched_df = enrich_with_expression(ranked_df)
        if "expression_category" in enriched_df.columns:
            expr_counts = enriched_df["expression_category"].value_counts()
            for cat, count in expr_counts.items():
                print(f"  {cat}: {count} candidates")
        ranked_df = enriched_df
        results["steps"]["expression"] = {"enriched": True}
    except (ImportError, Exception) as e:
        print(f"\n  [Skip] expression_filter: {e}")

    # =========================================================================
    # STEP 2d: Safety screening (if available)
    # =========================================================================
    try:
        from safety_screen import screen_candidates
        print(f"\n{'─' * 70}")
        print("  STEP 2d: Safety screening")
        print(f"{'─' * 70}")

        screened_df = screen_candidates(ranked_df)
        if "safety_flag" in screened_df.columns:
            safety_counts = screened_df["safety_flag"].value_counts()
            for flag, count in safety_counts.items():
                print(f"  {flag}: {count} candidates")
        ranked_df = screened_df
        results["steps"]["safety"] = {"screened": True}
    except (ImportError, Exception) as e:
        print(f"\n  [Skip] safety_screen: {e}")

    # =========================================================================
    # STEP 3: Vaccine Design
    # =========================================================================
    print(f"\n{'─' * 70}")
    print("  STEP 3: Vaccine construct design")
    print(f"{'─' * 70}")

    selected_df = design_mod.select_vaccine_epitopes(ranked_df, config)
    if len(selected_df) == 0:
        print("  WARNING: No epitopes selected. Using top candidates as fallback.")
        selected_df = ranked_df.head(10)

    epitopes = selected_df["mutant_peptide"].tolist()
    construct = design_mod.build_mrna_construct(epitopes, config)

    print(f"  Epitopes in construct: {construct['n_epitopes']}")
    print(f"  mRNA length: {construct['mrna_length_nt']} nt")
    print(f"  GC content: {construct['gc_content']:.1%}")
    print(f"  Estimated MW: {construct['estimated_molecular_weight_kda']} kDa")

    results["steps"]["vaccine_design"] = {
        "n_epitopes": construct["n_epitopes"],
        "mrna_length": construct["mrna_length_nt"],
        "gc_content": construct["gc_content"],
    }

    # =========================================================================
    # STEP 4: Benchmark Validation
    # =========================================================================
    print(f"\n{'─' * 70}")
    print("  STEP 4: Benchmark mutation validation")
    print(f"{'─' * 70}")

    # Check if known immunogenic mutations scored higher than passengers
    benchmark_lookup = {m["aa_change"]: m for m in BENCHMARK_MUTATIONS}

    allele_column = "hla_allele" if "hla_allele" in ranked_df.columns else None
    if allele_column is None:
        print("  WARNING: no hla_allele column in predictions. Restriction-aware")
        print("           scoring is unavailable and positives cannot be judged.")

    for aa_change, meta in benchmark_lookup.items():
        gene = meta["gene_symbol"]
        expected = meta["expected_immunogenic"]

        restriction = BENCHMARK_RESTRICTION.get(aa_change, {})
        restricting = restriction.get("restricting_alleles", [])
        evidence = restriction.get("evidence", "none" if expected else "n/a")

        matching = ranked_df[
            (ranked_df["gene_symbol"] == gene) &
            (ranked_df["aa_change"] == aa_change)
        ]

        if len(matching) == 0:
            print(f"  - {gene} {aa_change}: not found in predictions")
            results["benchmark_results"].append({
                "gene": gene, "mutation": aa_change,
                "expected_immunogenic": expected,
                "evidence": evidence,
                "status": "NOT_FOUND",
            })
            continue

        # A positive control is judged ONLY against its documented restricting
        # allele. Scoring it against the whole default panel asks the predictor
        # to present the epitope on molecules that cannot bind it.
        scored_on = None
        restricted = matching
        if expected and restricting and allele_column:
            restricted = matching[matching[allele_column].isin(restricting)]
            if len(restricted) == 0:
                # The restricting allele was never predicted against, so there
                # is nothing to judge. This is a panel gap, not a pipeline miss.
                missing = ", ".join(restricting)
                print(f"  ! {gene} {aa_change}: restricting allele not in panel "
                      f"({missing}); add it to score this control")
                results["benchmark_results"].append({
                    "gene": gene, "mutation": aa_change,
                    "expected_immunogenic": expected,
                    "evidence": evidence,
                    "restricting_alleles": restricting,
                    "status": "ALLELE_NOT_TESTED",
                })
                continue
            scored_on = restricted.iloc[0].get(allele_column)

        best = restricted.iloc[0]
        score = best.get("vaccine_priority_score", 0)
        ic50 = best.get("ic50_nM", 99999)
        classification = best.get("classification", "unknown")
        is_binder = classification in ("strong_binder", "weak_binder")

        # Positives with unsourced or computational-only evidence are not valid
        # controls. Report them, but never let them move the sensitivity figure.
        if expected and evidence not in SCOREABLE_EVIDENCE:
            status = "INDETERMINATE"
        elif expected:
            status = "PASS" if is_binder else "FAIL"
        else:
            status = "PASS" if not is_binder else "FAIL"

        results["benchmark_results"].append({
            "gene": gene,
            "mutation": aa_change,
            "expected_immunogenic": expected,
            "evidence": evidence,
            "restricting_alleles": restricting,
            "scored_on_allele": scored_on,
            "predicted_binder": is_binder,
            "best_ic50": round(ic50, 1),
            "priority_score": round(score, 1),
            "status": status,
        })

        marker = {"PASS": "✓", "FAIL": "✗", "INDETERMINATE": "~"}.get(status, "?")
        allele_note = f" on {scored_on}" if scored_on else ""
        print(f"  {marker} {gene} {aa_change}{allele_note}: IC50={ic50:.0f}nM, "
              f"score={score:.1f}, {classification} "
              f"(expected: {'immunogenic' if expected else 'passenger'}, "
              f"evidence: {evidence})")

    # Sensitivity and specificity are reported only over scoreable controls.
    scoreable = [
        r for r in results["benchmark_results"]
        if r["status"] in ("PASS", "FAIL")
    ]
    positives = [r for r in scoreable if r["expected_immunogenic"]]
    negatives = [r for r in scoreable if not r["expected_immunogenic"]]
    excluded = len(results["benchmark_results"]) - len(scoreable)

    results["benchmark_summary"] = {
        "scoreable_controls": len(scoreable),
        "excluded_controls": excluded,
        "positives_scored": len(positives),
        "positives_recovered": sum(1 for r in positives if r["status"] == "PASS"),
        "negatives_scored": len(negatives),
        "negatives_rejected": sum(1 for r in negatives if r["status"] == "PASS"),
        "note": (
            "Excluded controls are those with no documented HLA restriction, "
            "with unsourced or computational-only evidence, or whose restricting "
            "allele is absent from the configured panel. They are reported "
            "individually but do not contribute to these counts. A sensitivity "
            "figure over this few controls is indicative only."
        ),
    }

    print(f"\n  Scoreable controls: {len(scoreable)} "
          f"({excluded} excluded as unscoreable)")
    if positives:
        print(f"  Positives recovered: "
              f"{results['benchmark_summary']['positives_recovered']}/{len(positives)}")
    if negatives:
        print(f"  Negatives rejected:  "
              f"{results['benchmark_summary']['negatives_rejected']}/{len(negatives)}")

    # =========================================================================
    # STEP 5: Generate Validation Report
    # =========================================================================
    print(f"\n{'─' * 70}")
    print("  STEP 5: Generating validation outputs")
    print(f"{'─' * 70}")

    output_dir = Path("output/validation")
    output_dir.mkdir(parents=True, exist_ok=True)

    # Save benchmark results
    bench_df = pd.DataFrame(results["benchmark_results"])
    bench_path = output_dir / "benchmark_results.csv"
    bench_df.to_csv(bench_path, index=False)
    print(f"  Saved: {bench_path}")

    # Save binding predictions
    binding_path = output_dir / "validation_binding_predictions.csv"
    ranked_df.to_csv(binding_path, index=False)
    print(f"  Saved: {binding_path}")

    # Save selected epitopes
    epitopes_path = output_dir / "validation_selected_epitopes.csv"
    selected_df.to_csv(epitopes_path, index=False)
    print(f"  Saved: {epitopes_path}")

    # Save vaccine construct
    construct_json = {k: v for k, v in construct.items()
                      if k not in ("full_mrna_dna", "full_mrna_rna", "cassette_dna")}
    construct_path = output_dir / "validation_construct.json"
    with open(construct_path, "w") as f:
        json.dump(construct_json, f, indent=2)
    print(f"  Saved: {construct_path}")

    # Save mRNA FASTA
    fasta_path = output_dir / "validation_vaccine.fasta"
    with open(fasta_path, "w") as f:
        f.write(f">MM_benchmark_vaccine | {construct['n_epitopes']} epitopes | "
                f"{construct['mrna_length_nt']} nt\n")
        seq = construct["full_mrna_rna"]
        for i in range(0, len(seq), 80):
            f.write(seq[i:i+80] + "\n")
    print(f"  Saved: {fasta_path}")

    # Generate report
    report = design_mod.generate_vaccine_report(construct, selected_df, config)
    report_path = output_dir / "validation_vaccine_report.txt"
    with open(report_path, "w") as f:
        f.write(report)
    print(f"  Saved: {report_path}")

    # Save full results JSON
    results_path = output_dir / "validation_summary.json"
    with open(results_path, "w") as f:
        json.dump(results, f, indent=2, default=str)
    print(f"  Saved: {results_path}")

    # =========================================================================
    # Summary
    # =========================================================================
    passed = sum(1 for b in results["benchmark_results"] if b["status"] == "PASS")
    total = len(results["benchmark_results"])

    print(f"\n{'=' * 70}")
    print(f"  VALIDATION COMPLETE")
    print(f"{'=' * 70}")
    print(f"  Benchmark mutations: {passed}/{total} passed")
    print(f"  Pipeline steps completed: {len(results['steps'])}")
    print(f"  Vaccine construct: {construct['n_epitopes']} epitopes, "
          f"{construct['mrna_length_nt']} nt")
    print(f"\n  All outputs in: {output_dir}/")
    print(f"{'=' * 70}")

    return results


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Validate the neoantigen pipeline")
    parser.add_argument("--config", default="config.yaml", help="Config file")
    args = parser.parse_args()

    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    run_validation(args.config)
