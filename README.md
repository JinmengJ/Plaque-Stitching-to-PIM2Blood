# Plaque Stitching to PIM2Blood

Code accompanying the manuscript **“Reconstruction of a plaque progression cell atlas reveals monocyte state alerting near-term atherosclerotic plaque instability.”**

## Overview

Human atherosclerotic plaques are generally sampled only after clinically advanced disease has developed, whereas mouse models provide controlled temporal observations but do not fully reproduce advanced human plaque instability. This repository contains custom scripts used in two key components of the study:

1. **Cross-species diffusion stitching** — construction of an LDLR mouse temporal reference and projection of APOE mouse and human monocyte/macrophage states into the reference diffusion space.
2. **PIM2Blood panel development and evaluation** — discovery of a module-balanced blood marker panel from plaque instability monocyte (PIM)-associated genes, followed by panel-only evaluation using UK Biobank Olink proteomic data.

The repository provides analysis code only. Public input datasets should be downloaded from their original repositories. Controlled-access human data and UK Biobank participant-level data are not redistributed here.

## Repository contents

| File | Purpose |
| --- | --- |
| `cross_species_diffusion_stitching.R` | Builds the LDLR diffusion reference, projects APOE and human query cells, assigns reference-supported or reference-external stages, and exports mapping results. |
| `PIM2Blood_main.py` | Performs subject-level development splitting, NMF-based module discovery, candidate-gene filtering, module-balanced panel selection, and microarray evaluation. |
| `PIM2Blood_evaluation.py` | Evaluates a frozen panel-only score for incident AMI in UK Biobank Olink data without marker reselection or clinical-model fitting. |
| `DATA_INPUTS.md` | Required input files, columns, and object metadata. |
| `requirements.txt` | Python dependencies. |
| `R_REQUIREMENTS.md` | R dependencies and installation guidance. |

## Data availability

The public scRNA-seq datasets used in the study are available from GEO under accessions **GSE155512**, **GSE155513**, **GSE131778**, and **GSE184073**. The microarray datasets are available under accessions **GSE59867** and **GSE62646**.

Newly generated scRNA-seq and spatial transcriptomics raw sequence data from the study-specific human validation cohorts are available through the Genome Sequence Archive at the National Genomics Data Center, China National Center for Bioinformation/Beijing Institute of Genomics, Chinese Academy of Sciences, under BioProject **PRJCA068015** and GSA-Human accession **HRA019548**: <https://ngdc.cncb.ac.cn/gsa-human/browse/HRA019548>.

UK Biobank data are available to approved researchers through the UK Biobank Research Analysis Platform. This repository does not contain UK Biobank participant-level data, identifiers, derived participant-level scores, or Olink measurements.

## Installation

### Python

Python 3.10 or later is required.

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

### R

The cross-species workflow requires R and the Seurat, Matrix, destiny, and dplyr packages. See [`R_REQUIREMENTS.md`](R_REQUIREMENTS.md) for installation guidance.

## Input preparation

Detailed input requirements are provided in [`DATA_INPUTS.md`](DATA_INPUTS.md). In brief:

- `PIM2Blood_main.py` expects two gene-by-sample microarray tables represented as samples in rows and genes in columns, a PIM marker table, and optionally a UK Biobank/Olink header for cross-platform marker intersection.
- `PIM2Blood_evaluation.py` expects a frozen panel table and an approved local UK Biobank/Olink analysis table.
- `cross_species_diffusion_stitching.R` expects preprocessed Seurat objects for the LDLR reference and the APOE and human query datasets.

Raw CEL processing, initial scRNA-seq quality control, cell annotation, ortholog conversion, spatial transcriptomics processing, and construction of the approved UK Biobank analysis table are upstream procedures and are not performed by these three scripts.

## Usage

### 1. Cross-species diffusion stitching

Edit the `config` list near the beginning of `cross_species_diffusion_stitching.R` to specify the input directory, output directory, RDS filenames, assay, and metadata columns. Then run:

```bash
Rscript cross_species_diffusion_stitching.R
```

The workflow writes the frozen reference features, diffusion reference, mapping parameters, mapped Seurat objects, per-cell metadata, diffusion coordinates, and stage summaries to the configured output directory.

### 2. PIM2Blood panel construction

Edit the path constants near the beginning of `PIM2Blood_main.py`:

- `BASE`
- `DEV_RAW_FP`
- `EXT_RAW_FP`
- `PIMS_FP`
- `UKB_COHORT_FP`

Then run:

```bash
python PIM2Blood_main.py
```

The script saves the development split, NMF results, module membership, module and gene statistics, candidate panels, the frozen selected panel, panel scores, evaluation metrics, and a JSON method manifest beneath `SAVE`.

The random seed and analysis constants are defined at the top of the script. Any departure from the manuscript analysis should be documented before rerunning the workflow.

### 3. Frozen-panel UK Biobank evaluation

Run the evaluation only in an approved UK Biobank computing environment:

```bash
python PIM2Blood_evaluation.py \
  --panel-file /path/to/selected_panel_by_devtrain_external.csv \
  --ukb-file /path/to/approved_ukb_olink_table.csv \
  --out-dir /path/to/private/output_directory \
  --horizon-days 365
```

This script constructs an equal-weight mean z-score from the frozen panel and reports AUROC, AUPRC, and operating-point metrics at prespecified specificity targets. It does not reselect markers or fit a clinical prediction model.

**Do not commit any UK Biobank input or output file to this repository.** Several output tables contain participant-level identifiers or scores and must remain within the approved analysis environment.

## Reproducibility notes

- All longitudinal samples from a participant must remain in the same development split.
- Panel discovery and frozen-panel evaluation are separated into different scripts to reduce accidental marker reselection during outcome evaluation.
- The code contains dataset-specific file paths that must be replaced with local paths before execution.
- Large expression matrices, Seurat objects, controlled-access human data, and participant-level outputs are intentionally excluded from version control.
- Exact reproduction requires the same preprocessed inputs and metadata encodings described in `DATA_INPUTS.md`.

## Citation

If you use this repository, please cite the associated manuscript. Full citation information and a persistent software archive identifier will be added after publication.

## Contact

For questions about the code or study, please open a GitHub issue or contact the corresponding author listed in the manuscript.

## License

No reuse license has yet been specified. A license should be selected by the authors and their institution before public release, taking the associated pending patent application into account.
