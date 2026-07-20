# Input data requirements

This document describes the inputs consumed directly by the released scripts. It does not describe every upstream preprocessing step in the manuscript.

## `PIM2Blood_main.py`

### Development and external microarray tables

The files referenced by `DEV_RAW_FP` and `EXT_RAW_FP` must be comma-separated tables with one sample per row.

Required metadata columns:

| Column | Description |
| --- | --- |
| `SampleID` | Unique sample identifier. |
| `SubjectID` | Participant identifier used to keep longitudinal samples in the same split. |
| `Group` | Disease group. Values are normalized to `AMI` or `Stable`. |
| `Time` | Sampling time. Recognized aliases are normalized to `adm`, `dis`, `1m`, or `6m`. Stable samples are assigned `adm`. |

All remaining columns are interpreted as gene-level expression features. Gene symbols must be consistent across both microarray tables and the PIM marker table.

The released analysis assumes that raw CEL files have already been processed and collapsed to gene-level expression before these CSV files are created.

### PIM marker table

The CSV referenced by `PIMS_FP` must contain:

| Column | Description |
| --- | --- |
| `gene_in_data` | PIM-associated gene symbol represented in the processed microarray data. |

### Optional UK Biobank/Olink header

If `UKB_COHORT_FP` exists, its column names are used to restrict candidate genes to those with corresponding Olink measurements. Protein columns may be named with the gene symbol or with an `olink_` prefix; matching is case-insensitive after punctuation is removed.

The main construction script reads only the header at this stage. Do not move or publish the underlying controlled participant-level file.

## `PIM2Blood_evaluation.py`

### Frozen panel table

The file supplied through `--panel-file` must contain:

| Column | Required | Description |
| --- | --- | --- |
| `gene` | Yes | Frozen panel gene/protein symbol. |
| `panel_rank` | No | Panel order; generated if absent. |
| `module_use` | No | Module annotation retained in mapping output when present. |

The manuscript analysis uses a 20-marker panel. The script warns but continues if a different number of unique markers is supplied.

### Approved UK Biobank/Olink table

The file supplied through `--ukb-file` must contain:

| Column | Description |
| --- | --- |
| `EID_STD` | Participant identifier included in private score outputs. |
| `T0` | Baseline date/time field retained in private score outputs. |
| `followup_days` | Available follow-up duration after baseline. |
| `t_event_days` | Time from baseline to AMI; missing for participants without an observed event. |

It must also contain one Olink abundance column for every frozen panel marker. Optional columns used for subgroup summaries are `sex` and `age_recruit`.

All UK Biobank data and participant-level outputs must remain in the approved computing environment and must not be committed to GitHub.

## `cross_species_diffusion_stitching.R`

The input files are preprocessed Seurat RDS objects specified in the script-level `config` list.

### LDLR temporal reference

The reference object must contain:

- the configured RNA assay and normalized `data` slot;
- a stage metadata column specified by `ldlr_stage_col` (default: `TimePoints`);
- a reference cell-state/subcluster column specified by `ldlr_celltype_col` (default: `cellgroup`);
- the reductions required if the optional Seurat reference-mapping helper is used.

### APOE and human query objects

Each query object must contain:

- the configured RNA assay and normalized `data` slot;
- gene symbols harmonized with the LDLR reference;
- sufficient overlap with the locked reference feature set;
- for the final target-state summary, a query cell-state column specified by `query_celltype_col` (default: `predicted.celltype`).

Cross-species ortholog conversion, initial Seurat preprocessing, and construction of the query cell-state annotations must be completed upstream.

## Privacy and repository exclusions

Do not upload any of the following:

- UK Biobank participant-level inputs or outputs;
- controlled-access GSA-Human files;
- direct or indirect participant identifiers;
- local credentials, access tokens, or configuration files containing secrets;
- large RDS, h5ad, matrix, CEL, or intermediate result files unless redistribution is explicitly permitted.
