# Release checklist

Complete this checklist before making the repository public or citing it in the manuscript.

## Scientific consistency

- [ ] Confirm whether panel size was selected using the held-out Microarray-Dev set, as currently stated in the manuscript, or using Dev-train plus Microarray-External, as implemented in `PIM2Blood_main.py`.
- [ ] Confirm whether score-standardization parameters were fitted using Microarray-Dev training only or the combined Dev-train plus Microarray-External data.
- [ ] Confirm whether the external-cohort admission-to-discharge effect was part of the gene eligibility gate.
- [ ] Confirm whether within-module gene ranking used NMF loading alone or the composite `gene_score` implemented in the script.
- [ ] Confirm whether NMF directly identified four modules or whether four core modules were selected from a larger NMF solution.
- [ ] Confirm the final selected panel contains exactly 20 unique markers and add the frozen panel file if redistribution is permitted.
- [ ] Align the cross-species manuscript parameters and R script defaults, including feature selection, cells sampled per stage, diffusion dimensions, neighbor number, and fixed versus calibrated thresholds.
- [ ] Confirm whether Seurat label transfer is performed inside the released workflow or upstream in the supplied query RDS objects.

## Reproducibility

- [ ] Replace or parameterize local `/data2/...` paths.
- [ ] Add `sessionInfo()` from the original R environment and freeze the R environment with `renv.lock` if possible.
- [ ] Freeze exact Python package versions from the original environment.
- [ ] Test all scripts from a clean directory using shareable example or synthetic inputs.
- [ ] Add a small, non-sensitive example dataset or a schema-only example if permitted.
- [ ] Verify that every file named in the README is present.

## Privacy and permissions

- [ ] Confirm that no UK Biobank participant-level data, identifiers, or derived scores are tracked by Git.
- [ ] Confirm that no controlled GSA-Human data are tracked by Git.
- [ ] Scan the complete Git history for credentials, access tokens, participant identifiers, and restricted files.
- [ ] Confirm redistribution permissions for every included result or example file.

## Publication metadata

- [ ] Select a software license with institutional approval, considering the pending patent application.
- [ ] Add the final manuscript citation, DOI, and corresponding-author contact information.
- [ ] Create a versioned GitHub release matching the submitted manuscript code.
- [ ] Archive the release in Zenodo or another long-term repository and add the software DOI to the manuscript if desired.
