# Codes-GHJ-morphology
Reproducible workflow for the analysis of 3D geometric morphometrics data presented in the manuscript "Looking for Functional Signals in the Glenohumeral Joint of Primates using three-dimensional (3D) Geometric Morphometrics".

R code accompanying the manuscript:

> Cluzeau et al. (2026). *Looking for Functional Signals in the Glenohumeral Joint of Primates using three-dimensional (3D) Geometric Morphometrics*. 

This repository contains the R analysis pipeline used to generate the geometric
morphometric results reported in the manuscript, including Procrustes superimposition, PCA/RCPCA, allometry correction,
disparity analysis, phylogenetic comparative methods (BM/OU model fitting with mvMORPH),
and the associated figures.

## Repository contents

| File | Description |
|---|---|
| `Cluzeau_et_al_2026_RSOS.R` | Main analysis script (data import → GPA → PCA → allometry → disparity → phylogenetic models → figures) |
| `deformation_illustration.R` | Script generating the shape deformation figure(s) in the manuscript (warped meshes along PC axes) |
| `Disparity.R` | Script generating intergroup disparity |

## Attribution

A portion of `Cluzeau_et_al_2026_RSOS.R` [the multiblock analysis section] is adapted from the multiblock analysis approach described
in Thomas et al. (2023) [Constructing a multiple-part morphospace using a multiblock method]. This section is indicated with an inline comment in the code (`# adapted from Thomas et al. 2023`).

## Requirements

- R ≥ 4.x
- Packages Main analysis script: `morphoBlocks`, `geomorph`, `FactoMineR`, `abind`, `ape`, `geiger`, `mvMORPH`,
  `phytools`, `RColorBrewer`, `dplyr`, `rgl`
- Packages Deformation illustration script: `Morpho`, `geomorph`, `rgl`, `ape`, `RANN`, `rvcg`
- Packages Disparity script: `geomorph`, `ggplot2`, `reshape2`

`morphoBlocks` is installed from GitHub:
```r
remotes::install_github("aharmer/morphoBlocks", build_vignettes = TRUE, dependencies = TRUE)
```

## Reproducing the analysis

```r
source("Cluzeau_et_al_2026_RSOS.R")
```

The deformation figure is generated separately:
```r
source("deformation_illustration.R")
```

## Citation

If you use this code, please cite:

> Cluzeau et al. (2026). *Looking for Functional Signals in the Glenohumeral Joint of Primates using three-dimensional (3D) Geometric Morphometrics*. Submitted to Royal Society Open Science.

A citable archived version of this exact code will be available on Zenodo: DOI: TBD
## License

This code is released under the MIT License.

## References

- Thomas, D. B., Harmer, A. M. T., Giovanardi, S., Holvast, E. J., McGoverin, C. M., & Tenenhaus, A. (2023). Constructing a multiple-part morphospace using a multiblock method. Methods in Ecology and Evolution, 14, 65–76. https://doi.org/10.1111/2041-210X.13781 —
  method/code adapted for the multiblock analysis (see `Attribution` above).

## Contact

Cluzeau Océane — oceane.cluzeau@kuleuven.be - ORCID: 0009-0009-8183-2805
