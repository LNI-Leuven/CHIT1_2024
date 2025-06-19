# CHIT1 AT DIAGNOSIS IN MS - CODE
[![DOI](https://zenodo.org/badge/663491793.svg)](https://zenodo.org/doi/10.5281/zenodo.11235175)

This public repository contains the scripts used to reproduce the main results and figures from our article **CHIT1 at diagnosis predicts faster disability progression and reflects early microglial activation in multiple sclerosis** published in *Nature Communications* ([Beliën & Swinnen et al.; 2024](https://www.nature.com/articles/s41467-024-49312-y)). In this study we provided a clear rationale for further validation of CHIT1 as a prognostic biomarker in multiple sclerosis (MS). To do so, we applied mixed-effects models, a machine learning approach and single-cell/nucleus transcriptomics. The code for these analyses is here available and structured as follows:

**Statistical models**
- Correlation, single-time point and multi-time point analysis: `Biomarker and CHIT1 analyses.Rmd`

This is an R Markdown document writtin in R.

**Machine learning**
- `code.ipynb` is a Jupyter notebook (Python) showcasing the calculation of permutation importances for this paper.
- `code_utils.py` introduces some utility functions that are used in `code.ipynb`.
- `matplotlibrc` contains the plot generation setup.
- `requirements.txt` lists required Python packages to run the code yourself.

**Single-cell/nucleus transcriptomics**
- Quality control: `1. Quality_Control_CSF.Rmd` and `2. Quality_Control_CNS.Rmd`
- Integration: `3. Integration_With_Harmony.Rmd`
- Myeloid subclustering: `4. Myeloid_Subclustering.Rmd`
- Cluster annotation: `5. Cluster_Annotation.Rmd`
- Trajectory analysis: `6. Trajectory_Analysis.Rmd`
- Distribution and DGE CHIT1+ cells: `7. Distribution_And_DGE_CHIT1_Cells.Rmd`
- Pathway analysis: `8. Pathway_Analysis.Rmd`
- Reviewers' comments: `Reviewers' comments.Rmd`

These are R Markdown documents written in R.

## Contributors
This code was created by members of the [Laboratory for Neuroimmunology](https://gbiomed.kuleuven.be/english/research/50000666/50000668/50525530/laboratory-for-neuroimmunology) as well as the Department of Public Health and Primary Care at the KU Leuven (Belgium).
