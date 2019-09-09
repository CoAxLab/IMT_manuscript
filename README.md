Matlab code for IMT predict paper

# ===================================

# Scrape imaging data and generate samples

IMT_predict_paper_make_dataset.m
- constructs combined AHAB & PIP dataset for later analyses
- contains details on excluded subjects
- prints to screen sample sizes according to studies and tasks

IMT_paper_Table1_demographics.R
- takes output from IMT_predict_paper_make_dataset.m and pulls demographic information for all subjects from SPSS files
- produces descriptive stats of demographics (for Table 1)
- combines samples and writes demographic information to file

# Machine learning procedures

IMT_lassopcr.m
- central script to the entire paper
- further prepares dataset produced by IMT_predict_paper_make_dataset.m according to the task/condition of interest and other parameters 
- runs nested cross-validated LASSO-PCR
- calls lassopcr_cv.m within an outer cross-validation scheme 
	(5-fold for Faces; 2 folds for IAPS, according to study)
- specify mask for imaging data for running (amygdala) sufficient vs. necessary analyses

lassopcr_cv.m
- runs primary cross-validated LASSO-PCR routine given x = brain and y = outcome
- optimizes lambda via cross-validation, stratifying across y
- produces plots of 
	- lambda
	- # components retained
	- error metric (min-MSE, min-MedSE, max-r)

lassopcr_bootstrap.m
- takes the input (fmri_train) and output (stats_train) from a model trained in lassopcr_cv.m and draws bootstrap samples of observations, runs lassopcr_cv.m on these
- after bootstrapping, converts output to Z values and p-values
- also has option for nonparametric (percentile) p-values
- thresholds bootstrapped map 
- saves NIFTI images and cluster info

IMT_paper_figure.m
- generates brain figures for a given .nii image
- run this on :
	maineffects_fdr05_k50.nii
	weightmap_unthresholded.nii
	weightmap_bootstrap_p05_k50.nii

IMT_confounders_modifiers.m
- runs ancillary confounder and modifier analyses looking at age, sex, CMR

IMT_reliability.m
- calculates split-half Spearman-Brown corrected internal consistency
- multiple options:
	ROI average
	LASSO-PCR pattern
	generate voxelwise map
	
IMT_generalizability.m
- post-hoc tests of model generalizability across tasks
- can also apply PINES and SBP reactivity patterns to images

IMT_paper_Bayes_Factors.Rmd
- calculate BF for predicted-observed correlations

# ===================================

# Miscellaneous

cdtodrive.m
- changes directory to the ProjectDrive, flexible so that multiple users can use it and CD to the same place

lassopcr.m
- toy example of how the inner working of LASSO-PCR are handled in MATLAB
- not actually used much in analyses

lassopcr_cv_nested.m
- wrapper script that stratifies across y to call the lassopcr_cv within each outer fold
- produces multiple lambdas, # components retained, and fits for each fold
- bins holdout samples for each fold

amygdala_sufficient.m
- runs cross-validated LASSO-PCR using 2 predictor variables: mean left amygdala and mean right amygdala (different from masking and taking all amygdala voxels)
- this is different from use all voxels within an ROI mask

coverage.m
- calculate fMRI coverage from masks generated in SPM

surface.m
- tweaked from Tor Wager’s version to allow limbic left/right, and to show brainstem on the hi-res surface


