# Multiple Myeloma Neoantigen Vaccine Pipeline

A computational pipeline for exploring neoantigen candidates in multiple myeloma (MM), using open-access genomic data from the [MMRF CoMMpass Study](https://themmrf.org/we-are-curing-multiple-myeloma/mmrf-commpass-study/) via the [NCI Genomic Data Commons (GDC)](https://portal.gdc.cancer.gov/).

> **This is a research prototype, not a vaccine design tool.** It produces
> ranked candidate lists and a representative construct layout. It has had no
> wet-lab validation of any kind, and the sections below set out exactly which
> parts are implemented and which are placeholders. If you are a patient or
> carer who has found this repository while looking for treatment, nothing here
> is a therapy or a route to one; please talk to your haematology team.

## Status: what runs, and what does not

Honest accounting, because the previous version of this README did not
distinguish these and that made the repository hard to evaluate.

### Implemented and working

- GDC API retrieval of MMRF CoMMpass somatic mutations and clinical data
- Mutation filtering to protein-altering variants, with mutant peptide generation
- MHC-I binding prediction via MHCflurry, including class II peptide lengths
- Candidate ranking, epitope selection and construct assembly
- An epitope-level benchmark (`benchmark_epitopes.py`) that scores documented
  neoantigens against their own restricting HLA allele
- A Streamlit dashboard for browsing results

### Stubbed, approximate, or not implemented

| Component | Actual state |
| --- | --- |
| **Patient HLA typing** | Not implemented. `config.yaml` names arcasHLA and OptiType, but the code always falls back to population-frequency priors for European ancestry. **Nothing this pipeline outputs is personalised to an individual**, because the allele set is a demographic assumption. Real typing from the CoMMpass RNA-seq is the single most important missing piece. |
| **Built-in PSSM predictor** | A teaching-grade approximation with matrices for four alleles only (A\*02:01, A\*01:01, A\*03:01, B\*07:02). Not comparable to a trained predictor. Install MHCflurry; without it, affinity numbers are not interpretable. |
| **Ensemble weighting** | The 0.7/0.3 split between MHCflurry and PSSM is a choice made for this repository, not a validated scheme. Read the per-predictor columns. |
| **COSMIC and IEDB annotation** | Offline dictionaries hand-built in `external_validation.py`. No API calls are made, so the counts reported in validation summaries describe that local file, not a database query. |
| **3' UTR sequence** | A 43-nucleotide placeholder. It is labelled in `config.yaml` as the published AES-mtRNR1 design; it is not that sequence. Do not treat any construct output as a real UTR architecture. |
| **Self-similarity safety screen** | Compares against a short hardcoded gene list, not the human proteome. It is not a tolerance or autoimmunity screen in any meaningful sense. |
| **Secondary validation cohort** | None. Previously configured as "TCGA-MM", which does not exist. |
| **Selection pressure (dN/dS)** | Disabled. Previously carried fifteen invented per-gene constants. |

### Known limits of the approach

Multiple myeloma has a low mutational burden. Across 664 CoMMpass patients the
mean was 63.9 missense mutations, 23.5 predicted neoantigens and **9.4
expressed neoantigens** per patient
([Miller et al., Blood Cancer Journal, 2017](https://www.nature.com/articles/bcj201794)).
The default config asks for 5 to 20 epitopes per construct, so for a typical
patient there is barely a shortlist to rank. A better binding predictor does not
fix this; the constraint is antigen supply. Myeloma-specific sources that this
pipeline does **not** currently model, and that would matter more than predictor
accuracy, include the clonal immunoglobulin idiotype, IGH translocation
breakpoint peptides, and selection by clonal persistence across the serial
CoMMpass timepoints.

## Pipeline Overview

```
┌─────────────────┐     ┌──────────────────┐     ┌──────────────────┐     ┌──────────────────┐
│  01: Fetch Data  │────▶│  02: Parse Muts   │────▶│  03: Predict     │────▶│  04: Design      │
│  (GDC API)       │     │  & Gen Candidates │     │  MHC Binding     │     │  mRNA Vaccine    │
└─────────────────┘     └──────────────────┘     └──────────────────┘     └──────────────────┘
```

| Step | Script | Description |
|------|--------|-------------|
| 1 | `01_fetch_mmrf_data.py` | Fetches somatic mutations, clinical data, and gene expression metadata from GDC |
| 2 | `02_parse_mutations.py` | Filters protein-altering mutations, generates mutant peptides, scores immunogenicity |
| 3 | `03_predict_binding.py` | Predicts MHC-I binding affinity (PSSM or MHCflurry), ranks by vaccine priority |
| 4 | `04_design_vaccine.py` | Selects epitopes, optimizes ordering, builds full mRNA construct with UTRs/signal peptide |
| — | `05_run_pipeline.py` | Orchestrator that runs all steps end-to-end |

## Quick Start

```bash
# Install dependencies
pip install -r requirements.txt

# Run the full pipeline
python 05_run_pipeline.py

# Or run individual steps
python 01_fetch_mmrf_data.py --max-cases 50
python 02_parse_mutations.py
python 03_predict_binding.py
python 04_design_vaccine.py
```

### Pipeline Options

```bash
# Start from a specific step (skips earlier steps)
python 05_run_pipeline.py --step 3

# Skip data fetch if you already have the data
python 05_run_pipeline.py --skip-fetch

# Dry run (show what would execute)
python 05_run_pipeline.py --dry-run

# Custom case/mutation limits
python 05_run_pipeline.py --max-cases 100 --max-mutations 20000
```

## Output Files

| File | Description |
|------|-------------|
| `data/mmrf_cases.csv` | Clinical data (demographics, diagnoses) |
| `data/mmrf_mutations.csv` | Raw somatic mutations from GDC |
| `data/neoantigen_candidates.csv` | Filtered candidates with immunogenicity scores |
| `data/binding_predictions.csv` | MHC binding predictions and rankings |
| `output/vaccine_report.txt` | Full vaccine design report |
| `output/vaccine_mrna.fasta` | mRNA sequence in FASTA format |
| `output/vaccine_construct.json` | Construct metadata (epitopes, properties) |
| `output/selected_epitopes.csv` | Final selected epitopes |

## Configuration

Edit `config.yaml` to customize:

- **GDC settings** — project, data types, case limits
- **Mutation filters** — VAF threshold, consequence types, driver genes
- **Neoantigen prediction** — peptide lengths, HLA alleles, binding thresholds
- **Vaccine design** — max epitopes, linker sequence, UTRs, signal peptide, poly-A length

## How It Works

### Data Source
The pipeline uses the **MMRF CoMMpass Study** (~995 newly diagnosed MM patients) with whole genome/exome sequencing, RNA-seq, and clinical data — all accessible through the GDC open-access API.

### Neoantigen Selection
Mutations are filtered to protein-altering variants (missense, frameshift, indels) and scored by:
- Physicochemical distance between wildtype and mutant amino acids
- Known MM driver gene status (KRAS, NRAS, BRAF, TP53, etc.), by membership of the `driver_genes` list in `config.yaml`
- MHC-I binding affinity (IC50) across the configured HLA panel, which is a population prior rather than a patient's genotype
- Agretopicity (mutant vs wildtype binding ratio)
- Foreignness score

The dN/dS selection-pressure term is disabled. It previously added a bonus
derived from fifteen hardcoded per-gene constants that had no source, so any
ranking produced with it was shaped by invented numbers. To restore it, run
`dndscv` against the MMRF MAF and use its real output, or take driver status
from OncoKB or IntOGen.

### mRNA Vaccine Design
The construct output follows the general *layout* of a published mRNA vaccine
architecture. It is a schematic, not a manufacturable sequence:

- **5' Cap** — m7G Cap1 structure (recorded as an annotation, not at sequence level)
- **5' UTR** — Human alpha-globin derived
- **Signal peptide** — Human tPA leader sequence
- **Epitope cassette** — Ordered epitopes joined by GSG linkers
- **3' UTR** — **placeholder, 43 nt.** `config.yaml` labels this as the
  AES-mtRNR1 combination used in published BioNTech constructs. It is not that
  sequence. Every construct this pipeline emits carries a stub here.
- **Poly-A tail** — 120 nucleotides
- **Modification** — N1-methylpseudouridine (m1Ψ) substitution is recorded as an
  annotation; no chemistry is modelled
- **Codon optimization** — Human codon usage bias, with optional GC balancing

## Dependencies

- Python 3.8+
- pandas, numpy, biopython
- requests (GDC API access)
- scikit-learn
- PyYAML, tqdm
- **[MHCflurry](https://github.com/openvax/mhcflurry) — effectively required.** It is
  listed as optional because the pipeline will run without it, but the built-in
  PSSM fallback covers four alleles and is a rough approximation, so affinity
  numbers produced without MHCflurry should not be interpreted:

  ```bash
  pip install mhcflurry
  mhcflurry-downloads fetch models_class1_presentation
  ```

## Validating a change

Run the epitope-level benchmark. It scores documented neoantigens against the
HLA allele each is actually restricted to, and takes seconds:

```bash
python benchmark_epitopes.py
```

At the time of writing it recovers 2 of 4 experimentally-sourced epitopes at the
500 nM threshold. Read the misses rather than the ratio. One is KRAS G12D on
HLA-C\*08:02, where the clinical evidence is strong but affinity predictors are
substantially weaker on HLA-C than on HLA-A and HLA-B for want of training data.
Four scoreable epitopes make this a smoke test, not an accuracy measurement.

The fuller `06_validate_pipeline.py` run is slower and exercises the whole chain.
Both now score each positive control against its own restricting allele; an
earlier version scored every control against a fixed six-allele panel that did
not contain several of them, which produced seven spurious failures.

## Live Dashboard

The dashboard is hosted at **https://mm-vaccine.streamlit.app** via Streamlit Community Cloud (free tier).

### Keeping it awake

Streamlit Cloud free tier apps sleep after ~7 days of inactivity and can take 30-60 seconds to wake on first visit. To keep the dashboard consistently responsive:

**Option A — UptimeRobot (free, recommended)**
1. Create a free account at [uptimerobot.com](https://uptimerobot.com)
2. Add a new HTTP(S) monitor pointing to: `https://mm-vaccine.streamlit.app`
3. Set check interval to 10 minutes
4. This pings the app regularly and prevents it sleeping

**Option B — Streamlit Community Cloud paid tier**
Upgrade at [share.streamlit.io](https://share.streamlit.io) for always-on hosting.

## Disclaimer

This pipeline is for **research and educational purposes only**. The vaccine designs are computationally generated and have not been experimentally validated. Clinical application requires extensive preclinical testing, regulatory approval, and patient-specific HLA typing. This tool is not a medical device and should not be used for clinical decision-making.

## References

- MMRF CoMMpass Study: https://themmrf.org/we-are-curing-multiple-myeloma/mmrf-commpass-study/
- GDC Data Portal: https://portal.gdc.cancer.gov/
- GDC API Documentation: https://docs.gdc.cancer.gov/API/Users_Guide/Getting_Started/
- MHCflurry: O'Donnell et al., Cell Systems 2020
- BioNTech individualized neoantigen vaccines: Sahin et al., Nature 2017
