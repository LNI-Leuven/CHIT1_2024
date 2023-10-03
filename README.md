# CHIT1 AT DIAGNOSIS IN MS - CODE
This public repository contains the R code scripts to reproduce the main results and figures from our article **CHIT1 at diagnosis predicts faster disability progression and reflects early microglial activation in multiple sclerosis** published in *Science Translational Medicine* (Beliën et al.; 2023). In this study we provided a clear rationale for further validation of CHIT1 as a prognostic biomarker in multiple sclerosis (MS). To do so, we applied mixed-effects models, a machine learning approach and single-cell/nucleus transcriptomics. The R code for these analyses is here available and structured as follows:

**Statistical models**
- Correlation analysis:
- Single-time point analysis:
- Multi-time point analysis:

**Machine learning**
- `code.ipynb` is a Jupyter notebook (Python) showcasing the calculation of permutation importances for this paper.
- `code_utils.py` introduces some utility functions that are used in `code.ipynb`.
- `matplotlibrc` contains the plot generation setup.
- `requirements.txt` lists required Python packages to run the code yourself.

**Single-cell/nucleus transcriptomics**

These are R Markdown documents written in R.
- Quality control: `1. Quality_Control_CSF.Rmd` and `2. Quality_Control_CNS.Rmd`
- Integration: `3. Integration_With_Harmony.Rmd`
- Myeloid subclustering: `4. Myeloid_Subclustering.Rmd`
- Cluster annotation: `5. Cluster_Annotation.Rmd`
- Trajectory analysis: `6. Trajectory_Analysis.Rmd`
- Distribution and DGE CHIT1+ cells: `7. Distribution_And_DGE_CHIT1_Cells.Rmd`
- Pathway analysis: `8. Pathway_Analysis.Rmd`

## Contributors:
This code was created by members of the [Laboratory for Neuroimmunology](https://gbiomed.kuleuven.be/english/research/50000666/50000668/50525530/laboratory-for-neuroimmunology) as well as the Department of Public Health and Primary Care at the KU Leuven (Belgium).
