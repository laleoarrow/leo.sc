<div align="center">

# LEO.SC

*Single-cell analysis and visualization suite for the LEO universe.*

<p align="center">
  <a href="https://laleoarrow.r-universe.dev/leo.sc">
    <img src="https://img.shields.io/badge/dynamic/json?label=leo.sc&query=$.Version&url=https://laleoarrow.r-universe.dev/api/packages/leo.sc&style=for-the-badge&color=226e63&logo=r&logoColor=white" />
  </a>
  <a href="https://github.com/r-universe/laleoarrow/actions/workflows/build.yml">
    <img src="https://img.shields.io/github/actions/workflow/status/r-universe/laleoarrow/build.yml?branch=master&label=R-universe&style=for-the-badge&logo=r&logoColor=white" />
  </a>
  <a href="https://github.com/laleoarrow/leo.sc/actions/workflows/R-CMD-check.yaml">
    <img src="https://img.shields.io/github/actions/workflow/status/laleoarrow/leo.sc/R-CMD-check.yaml?branch=main&label=R-CMD-check&style=for-the-badge&logo=github" />
  </a>
  <a href="https://lifecycle.r-lib.org/articles/stages.html#experimental">
    <img src="https://img.shields.io/badge/Lifecycle-Experimental-339999?style=for-the-badge&logo=jekyll&logoColor=white" />
  </a>
</p>

<img src="man/figures/image-1.png" width="90%" />

*A focused exploratory toolkit for single-cell data analysis, visualization, and in-silico perturbation.*

---

[**Documentation**](https://laleoarrow.github.io/leo.sc/) | [**Tutorial**](https://laleoarrow.github.io/leo.sc/articles/tutorial-sko.html) | [**Author**](https://www.researchgate.net/profile/Ao-Lu-5)

</div>

### Quick Start


You can install `leo.sc` from [r-universe](https://laleoarrow.r-universe.dev) like so:

``` r
# Enable the r-universe repository
options(repos = c(
  laleoarrow = "https://laleoarrow.r-universe.dev",
  CRAN = "https://cloud.r-project.org"
))

# Install leo.sc
install.packages("leo.sc")
```

Alternatively, you can install the development version from [GitHub](https://github.com/laleoarrow/leo.sc) with:

``` r
# install.packages("devtools")
devtools::install_github("laleoarrow/leo.sc")
```

### How to Contribute
We encourage and welcome contributions to LEO from the community. If you are interested in contributing code or documentation, please contact me by submitting a issue.

### License
The `leo.sc` package is open-source software. Please refer to the [LICENSE.md](LICENSE.md) file for terms and restrictions (GPL-3.0 License).

### Contact Information
For more information or assistance, please contact us at [Ao Lu](mailto:luao@stu.cqmu.edu.cn).

### Citation
If you use `leo.sc` in your work, please cite:

```r
citation("leo.sc")
```

```bibtex
@Manual{leo.sc,
  title = {leo.sc: Layered Exploratory Omics (LEO for single-cell analysis)},
  author = {Ao Lu},
  year = {2026},
  note = {R package version 0.1.0},
  url = {https://github.com/laleoarrow/leo.sc}
}
```

**We look forward to collaborating with researchers and developers worldwide to advance innovation and progress in genomic data analysis!**

---

<p align="center">
  <a href="https://www.researchgate.net/profile/Ao-Lu-5">
    <img src="https://img.shields.io/badge/ResearchGate-Ao%20Lu-00CCBB?logo=researchgate&logoColor=white&style=for-the-badge" />
  </a>
  <a href="https://orcid.org/0009-0001-0927-4468">
    <img src="https://img.shields.io/badge/ORCID-0009--0001--0927--4468-A6CE39?logo=orcid&logoColor=white&style=for-the-badge" />
  </a>
  <a href="https://github.com/laleoarrow/leo.sc">
    <img src="https://img.shields.io/badge/GitHub-laleoarrow-181717?logo=github&logoColor=white&style=for-the-badge" />
  </a>
</p>
