# LEO.SC <img src="man/figures/image-1.png" align="right" height="120" alt="leo.sc banner" />

<!-- badges: start -->
[![Dev status](https://img.shields.io/badge/dev%20status-experimental-orange.svg)](https://github.com/laleoarrow/leo.sc)
[![Version](https://img.shields.io/badge/version-0.1.0-blue.svg)](https://github.com/laleoarrow/leo.sc/releases)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![License: GPL-3](https://img.shields.io/badge/license-GPL--3-blue.svg)](LICENSE.md)
[![R](https://img.shields.io/badge/R-%3E%3D%204.2-blue.svg)](https://cran.r-project.org/)
<!-- badges: end -->

`leo.sc` is an exploratory toolkit for single-cell analysis in the LEO ecosystem. It provides practical helpers for QC, visualization, metadata processing, cell-type-level analysis, and in-silico perturbation experiments.

The package includes `silico_ko()`, a virtual knock-out / knock-in workflow that contrasts high-vs-low expressing cells of one target gene under matched cell-type composition, then runs differential expression and enrichment analysis.

<div align="center">
  <img src="man/figures/image-1.png" alt="leo.sc banner" width="900"/>
</div>

## Tutorial

SKO demo tutorial (Seurat built-in data): [articles/tutorial-sko.html](https://laleoarrow.github.io/leo.sc/articles/tutorial-sko.html)

## Installation

```r
options(repos = c(
  laleoarrow = "https://laleoarrow.r-universe.dev",
  CRAN = "https://cloud.r-project.org"
))
install.packages("leo.sc")
```

```r
remotes::install_github("laleoarrow/leo.sc")
```

## Links

- [Seurat](https://satijalab.org/seurat/)
- [SeuratObject `pbmc_small` dataset](https://rdrr.io/cran/SeuratObject/man/pbmc_small.html)
- [GitHub repository](https://github.com/laleoarrow/leo.sc)

## Citation

Lu A (2026). `leo.sc`: Layered Exploratory Omics (LEO for single-cell analysis). R package version 0.1.0, <https://github.com/laleoarrow/leo.sc>.

```bibtex
@Manual{leo.sc,
  title = {leo.sc: Layered Exploratory Omics (LEO for single-cell analysis)},
  author = {Ao Lu},
  year = {2026},
  note = {R package version 0.1.0},
  url = {https://github.com/laleoarrow/leo.sc}
}
```
