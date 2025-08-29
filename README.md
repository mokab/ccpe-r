# CCPE-R
Cell-Cycle Pseudotime Estimation in R

CCPE-R is an R re-implementation of CCPE (Cell-Cycle Pseudotime Estimation) by Liu et al. for single-cell RNA-seq analysis. CCPE maps high-dimensional expression profiles onto a 3-D helix, where the angular component captures cell-cycle phase (e.g., G1/S/G2M) and the axial progression provides a continuous pseudotime. This R port mirrors the original MATLAB algorithm’s alternating procedure—low-dimensional embedding coupled with helix fitting—so that cell-cycle structure can be quantified and visualized within typical R workflows.

## Original reference
Liu, J., Yang, M., Zhao, W., & Zhou, X. (2021). CCPE: cell cycle pseudotime estimation for single-cell RNA-seq data. Nucleic Acids Research.

## Original implementation (MATLAB)
[LiuJJ0327/CCPE](https://github.com/LiuJJ0327/CCPE)
