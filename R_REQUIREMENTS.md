# R environment

The cross-species diffusion stitching script requires:

- R 4.2 or later;
- Seurat 4.4.0 for closest alignment with the manuscript environment;
- Matrix;
- destiny;
- dplyr.

An example installation is:

```r
install.packages(c("Seurat", "Matrix", "dplyr"))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("destiny")
```

Package APIs can differ across Seurat major versions. For exact reproduction, record the output of `sessionInfo()` from the original analysis environment and use the same package versions. An `renv.lock` file should be generated from that environment before the archival release if it remains available.
