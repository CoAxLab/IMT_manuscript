
% test if the amygdala is sufficient to predict IMT

% this follows Luke's paper where he treated reach region as one value
% this uses 2 amygdala values (left/right) per subject in model building and evaluation

% requires loading respective fmri_train and fmri_test images for a given analysis

% this is DIFFERENT from looking at all voxels in an ROI

amyg_file = '/home/tek31/ProjectDrive/Users/Thomas/roi/VisceralControlCircuitROIs/Amygdala.nii';
amyg_region = region(fmri_data(amyg_file));

load('data.mat')


imagefiles = table2cell(dat(:, contains(dat.Properties.VariableNames, 'filepath')));
fmri_dat = fmri_data(imagefiles, out.mask);

amyg = extract_data(amyg_region, fmri_dat);
amyg = [amyg.dat];

%amyg_test = extract_data(amyg_region, fmri_test);
%amyg_test = [amyg_test.dat];

stats = lassopcr_cv(amyg, dat.imt, 'nested');

stats_train = lassopcr_cv(amyg_train, fmri_train.Y)
saveas(gcf, 'train_predicted_observed_amygdala_roi.png')

y_test_predicted = stats_train.constant + (amyg_test * stats_train.weight_obj.dat);
figure; scatter(y_test_predicted, y_test)
[test_corr, test_pval] = corr(y_test, y_test_predicted);
test_corr_ci = bootci(10000, @corr, y_test, y_test_predicted);
xlabel(sprintf('predicted %s in test sample', myvar_string)); 
ylabel(sprintf('observed %s in test sample', myvar_string)); lsline;
title(sprintf('r(predicted,observed) = %.2f [%.2f %.2f], p = %.3f', ...
    test_corr, test_corr_ci(1), test_corr_ci(2), test_pval));
saveas(gcf, 'test_predicted_observed_amygdala_roi.png')
