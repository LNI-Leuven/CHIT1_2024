# CHIT1 AT DIAGNOSIS IN MS - CODE
This public repository contains the R code scripts to reproduce the main results and figures from our article **CHIT1 at diagnosis reflects early microglial activation and predicts future disability progression in MS** published in *Science Translational Medicine* (Beliën et al.; 2023). In this study we provided a clear rationale for further validation of CHIT1 as a prognostic biomarker in multiple sclerosis (MS). To do so, we applied mixed-effects models, a machine learning approach and single-cell/nucleus transcriptomics. The R code for these analyses is here available and structured as follows:

**Statistical models**
- Correlation analysis:
- Single-time point analysis:
- Multi-time point analysis:

**Machine learning**
- `code.ipynb` is a Jupyter notebook (Python) showcasing the calculation of permutation importances for this paper.
- `code_utils.py` introduces some utility functions that are used in `code.ipynb`.
- `matplotlibrc` contains the plot generation setup.
- `requirements.txt` lists required Python packages to run the code yourself.
- `results` is a directory with figures, arrays and models resulting from `code.ipynb`.

**Single-cell/nucleus transcriptomics**
- Quality control: `1. Quality_Control_CSF.rmd` and `2. Quality_Control_CNS.rmd`
- Integration: `3. Integration_With_Harmony.rmd`
- Myeloid subclustering: `4. Myeloid_Subclustering.rmd`
- Cluster annotation: `5. Cluster_Annotation.rmd`
- Trajectory analysis: `6. Trajectory_Analysis.rmd`
- Distribution and DGE CHIT1+ cells: `7. Distribution_And_DGE_CHIT1_Cells.rmd`
- Pathway analysis: `8. Pathway_Analysis.rmd`

## Contributors:
This code was created by members of the [Laboratory for Neuroimmunology](https://gbiomed.kuleuven.be/english/research/50000666/50000668/50525530/laboratory-for-neuroimmunology) at the KU Leuven (Belgium).
