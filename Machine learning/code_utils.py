import copy
import numpy as np
import pandas as pd
from sklearn.model_selection import KFold
from scipy.interpolate import interp1d

biomarkernames = [
    "GPNMB_log","CCL18_log","CHIT1_log","NfL_log","sTREM2_log","YKL40_log"
]
paper_varnames = {
    "Disease_course": "Disease course",
    "GPNMB_log": "GPNMB",
    "CCL18_log": "CCL18",
    "CHIT1_log": "CHIT1",
    "NfL_log": "NfL",
    "sTREM2_log": "sTREM2",
    "YKL40_log": "CHI3L1",
    "Age_Onset": "Age at onset",
    "Age_LP": "Age at diagnosis",
    # "Diff_LP_onset": "Years between diagnosis and onset",
    "Diff_LP_onset": "Disease duration at LP",
    "Year_LP": "Year of diagnosis",
    "EDSS_t0": "EDSS at diagnosis",
    "WBC_count__ul_": "WBC count",
    "Albumin_CSF_mg_L_": "Albumin",
    "IgG_index": "IgG index",
    "OCB_count_labels": "OCB status",
    "OCB_count": "OCB count",
    "Gender": "Sex",
}

def load_data():
    """Load and construct all the relevant datasets.

    @return A 5-tuple:
        - X: The input dataset that can be used for modeling.
        - y: The targets per patient, i.e. average EDSS per bin.
        - biom: Table with patient characteristics.
        - edss: Table with EDSS measurements.
        - binedges: Edges of the bins on which y was based.
    """
    # =========LOAD EXISTING DATA=========
    # Load in the clean biomarker and EDSS tables
    data_folder = "/mnt/j/GBW-0089_Neuroimmunology/0006_AI/Data/"
    CLEAN_BIOM  = data_folder + "biom.csv"
    CLEAN_EDSS  = data_folder + "edss.csv"
    biom = pd.read_csv(CLEAN_BIOM)
    edss = pd.read_csv(CLEAN_EDSS)
    biom = biom[biom.Cohort == "Biomarker"].drop(columns="Cohort")
    edss = edss[edss.Cohort == "Biomarker"].drop(columns="Cohort")

    # Fix data types
    for var in ["Birth_date", "Date_Onset", "Datum_LP"]:
        biom[var] = pd.to_datetime(biom[var])
    edss = edss.astype({"EDSS_Measurement": pd.Int32Dtype()})
    for var in ["Birth_date", "Date_Onset", "Date_EDSS"]:
        edss[var] = pd.to_datetime(edss[var])
    biom.loc[biom.Disease_course.isna(),"Disease_course"] = pd.NA # for backwards compatibility

    # =========DROPPING PATIENTS=========
    # Drop patients with no EDSS measurements
    ids = edss.ID[edss.EDSS.isna()]
    assert all(pd.isna(edss.EDSS[edss.ID.isin(ids)])) # Are all visits effectively NA?
    edss = edss[~edss.ID.isin(ids)]
    biom = biom[~biom.ID.isin(ids)]

    # Drop patients with any of the following dates missing:
    biom = biom.dropna(subset=["Birth_date","Date_Onset","Datum_LP"])
    edss = edss.dropna(subset=["Birth_date","Date_Onset","Datum_LP"])

    # Reset the index after all this
    biom.reset_index(drop=True, inplace=True)
    edss.reset_index(drop=True, inplace=True)

    # =========CONSTRUCTING THE TARGET=========
    binedges = np.arange(11) # Locations of the bins
    targetnames  = [f"bin{i+1}" for i in range(len(binedges)-1)]
    targets = pd.DataFrame(dtype=float, index=biom.ID, columns=targetnames)
    for ID in biom.ID:
        x, y = get_time_series(edss, ID)
        binned_targets = build_binned_targets(x, y, binedges)
        for ibin in range(len(binned_targets)):
            targets.loc[ID, f"bin{ibin+1}"] = binned_targets[ibin]

    # =========CLEANING UP=========
    # Add targets and initial EDSS to biom
    biom.index = biom.ID
    biom[targetnames] = targets
    biom.reset_index(drop=True, inplace=True)

    # Reset index and extract X and y
    X = biom.drop(columns=["ID","Birth_date","Date_Onset","Datum_LP",*targetnames])
    y = biom[targetnames]
    X.reset_index(drop=True, inplace=True)
    y.reset_index(drop=True, inplace=True)
    return X, y, biom, edss, binedges

def get_time_series(edss, ID, artificial_timepoint=True):
    """Extract the visits for the given patient.
    
    @param edss: Table with EDSS measurements.
    @param ID: ID of the patient of interest.
    @param artificial_timepoint: Whether to include an artificial visit of
        EDSS=0 at disease onset.
    @return: Two arrays, the first one indicating the time after diagnosis in
        years and the second one giving the EDSS measurements at these times.
    """
    subset = edss.loc[edss.ID == ID]#.sort_values(by="Date_EDSS")
    x = subset["Years_Since_LP"].values
    y = subset["EDSS"].values

    if artificial_timepoint:
        # Add artificial time point: patient assumed to start at EDSS 0
        x = np.insert(x, 0, (subset["Years_Since_LP"]-subset["Years_Since_Onset"]).values[0])
        y = np.insert(y, 0, 0)
    return x, y

def build_binned_targets(x, y, binedges):
    """Build the binned targets, as described in doi.org/10.1007/978-3-031-34344-5_3
    
    @param x: Array with time after diagnosis in years.
    @param y: Array with EDSS scores corresponding to x.
    @param binedges: Array with the desired edges for the bins.
    @return: The #(binedges)-1 targets for this patient.
    """
    targets = []
    isort = np.argsort(x)
    x = x[isort]
    y = y[isort]
    for lbin, rbin in zip(binedges[:-1], binedges[1:]):
        # Gather all the data that can possibly be used for this bin
        iloc   = (lbin < x) & (x < rbin)
        before = (x <= lbin)
        after  = (x >= rbin)
        if any(before): # Add the last measurement before this bin (if any)
            iloc[np.argmin(before) - 1] = True
        if any(after): # Add the first measurement after this bin (if any)
            iloc[np.argmax(after)] = True
        xbin = x[iloc]
        ybin = y[iloc]

        # Continue if there is not enough data for this bin (e.g. only data before)
        if len(xbin) <= 1:
            continue
        
        # Calculate the weight for each section (pair of visits) within the bin
        clamp = lambda xxx: min(max(xxx, lbin), rbin)
        weights = [clamp(xbin[i+1]) - clamp(xbin[i]) for i in range(len(xbin)-1)]
        weights = [w / sum(weights) for w in weights]

        # Compute the section values per target
        slope = [(ybin[i+1]-ybin[i]) / (xbin[i+1]-xbin[i]) for i in range(len(xbin)-1)]
        means = [(ybin[i]+ybin[i+1])/2                     for i in range(len(xbin)-1)]

        # Correction for values outside the bin borders (only for some targets)
        if xbin[0] < lbin: # if any(before)
            means[0] += slope[0]*(lbin-xbin[0])/2 # Correct ybin[i] a bit
        if rbin < xbin[-1]: # if any(after)
            means[-1] -= slope[-1]*(xbin[-1]-rbin)/2 # Correct ybin[i] a bit

        # Take the average over all sections in this bin
        targets.append( np.dot(weights, means) )

    return targets


class CrossValidation:
    """Utility class providing cross-validation functionality."""

    def __init__(self, k=5, random_state=42):
        """Initializes this CrossValidation object.

        @param k: Number of folds in the cross-validation/
        @param random_state: Seed for RNG.
        """
        self.kfold = KFold(n_splits=k, shuffle=True, random_state=random_state)
        self.random_state = random_state

    def fit(self, X, y, model):
        """Fits the cross-validation object to the given data.

        @param X: Input dataframe for the model
        @param y: Target(s) that the model should predict.
        @param model: Object providing .fit(X,y) and .predict(X) functionality.
        @post: This object now has a y_pred attribute with cross-validated test
            set predictions.
        """
        y = pd.DataFrame(y)
        # self.X = X
        # self.y_true = y
        self.y_pred = y.copy(deep=True)
        self.y_pred[:] = pd.NA
        self.models = []
        for train_index, test_index in self.kfold.split(X):
            X_train = X.iloc[train_index].reset_index(drop=True)
            X_test  = X.iloc[test_index ].reset_index(drop=True)
            y_train = y.iloc[train_index].reset_index(drop=True)
            model.fit(X_train, y_train)
            self.models.append(copy.deepcopy(model))
            self.y_pred.iloc[test_index] = model.predict(X_test)
        return self
    
    def predict_altered_dataset(self, X):
        """Provide cross-validated test set predictions for an altered input dataset.

        @param X: Dataframe with the same patient on each row as the dataframe
            on which the model was trained.
        @return: Cross-validated test set predictions for X.
        """
        assert X.shape[0] == self.y_pred.shape[0]
        y_pred = self.y_pred.copy(deep=True)
        y_pred[:] = pd.NA
        for i, (train_index, test_index) in enumerate(self.kfold.split(X)):
            X_test = X.iloc[test_index ].reset_index(drop=True)
            y_pred.loc[test_index, :] = self.models[i].predict(X_test).values
        return y_pred


class RegressorChain:
    """Class of regression chains for predicting sequential targets."""

    def __init__(self, model):
        """Initializes this regressor chain.
        
        @param model: An initialized object providing .fit(X,y) and .predict(X)
            functionality.
        """
        self.model = model

    def fit(self, X, y):
        X = pd.DataFrame(X)
        y = pd.DataFrame(y)
        n_targets = y.shape[1]
        self.models = []
        self.y_columns = y.columns
        Xext = X.copy()
        for i in range(n_targets):
            self.model.fit(Xext, pd.DataFrame(y.iloc[:,i]))
            self.models.append(copy.deepcopy(self.model))
            yi_pred = self.model.predict(Xext)
            Xext[self.y_columns[i]] = yi_pred # Use predictions as new inputs
            # Xext[self.y_columns[i]] = y.iloc[:,i]
        return self

    def predict(self, X):
        assert len(self.models) == len(self.y_columns)
        X = pd.DataFrame(X)
        Xext = X.copy()
        for i, model in enumerate(self.models):
            yi_pred = model.predict(Xext)
            Xext[self.y_columns[i]] = yi_pred
        return Xext.iloc[:, -len(self.models):]


def compute_MAE(y_true, y_pred):
    """Compute the negative mean absolute error between true values and pred."""
    return -np.nanmean(np.abs(y_true - y_pred))

def compute_rho(y_true, y_pred):
    """Compute Pearson correlation between true values and pred."""
    i_nona = np.logical_and(~np.isnan(y_true), ~np.isnan(y_pred))
    return np.corrcoef(np.array(y_true)[i_nona].flatten(), 
                       np.array(y_pred)[i_nona].flatten())[0,1]

def pred_at_t(t, binedges, y_pred, interp_kind="linear"):
    """Interpolate the given binned predictions to any timepoint.

    @param t: Timepoint(s) of interest. Should be within the bin edges.
    @param binedges: Edges of the bins.
    @param y_pred: Array with bin predictions. Any missing values should only be 
        at the end of this array.
    @param interp_kind: How to interpolate the predictions. Possible options are
        a.o. 'linear' and 'quadratic'. See `scipy.interpolate.interp1d` for more 
        information.
    @return: Interpolations of y_pred to timepoints given in t.
    """
    bincenters = (binedges[:-1] + binedges[1:]) / 2

    nona = ~np.isnan(y_pred)
    y_pred     = y_pred[nona]
    bincenters = bincenters[nona]
    i_after_last = (t > (bincenters[-1] + (bincenters[0] - binedges[0])))

    f = interp1d(bincenters, y_pred, kind=interp_kind)
    t = t.copy()
    t[t < bincenters[ 0]] = bincenters[0] # or edss_t0
    t[t > bincenters[-1]] = bincenters[-1]
    interpolated = f(t)
    interpolated[i_after_last] = np.nan
    interpolated
    return interpolated
