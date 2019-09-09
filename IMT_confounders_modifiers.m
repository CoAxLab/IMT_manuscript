function lm = IMT_confounders_modifiers(task)

%% Given a fitted dataset, run post-hoc analyses looking at CMR and other mediators/moderators
%       Need to first run out = IMT_lassopcr('task')
%       Go to folder containing the 'data.mat' output from IMT_lassopcr

if nargin < 1, task = 'IAPS'; end 
close all;

%% CD to ProjectDrive
mywd = cdtodrive; cd('AHAB_II/ML_projects/IMT_LassoPCR/');

%% Load dataset table
alldat = readtable('datasets/ID_IMT_demographics.csv', 'TreatAsEmpty', 'NA');

%% load task-specific data
switch task
    case 'IAPS'
        cd IAPS
    case 'OF'
        cd Faces/OF/AllFaces
end

load('data.mat'); 

%% Calculate CMR - mean of standardized variables (reverse-code HDL)
% calculate across all participants

%alldat.waist_circumference(alldat.waist_circumference > 700) = NaN;
%alldat.HDL(alldat.HDL > 700) = NaN;
%alldat.triglycerides(alldat.triglycerides > 700) = NaN;
%alldat.glucose(alldat.glucose > 700) = NaN;

alldat.CMR = nanmean([scale(alldat.waist) scale(-1*alldat.HDL) scale(alldat.triglycerides) scale(alldat.glucose) scale(alldat.SBP)], 2);

%% join yfit values with master dataset
yfit = table(out.ID, out.yfit);
yfit.Properties.VariableNames = {'ID', 'IMT_predicted'};
dat = innerjoin(alldat, yfit);
dat.IMT_observed = dat.IMT;

%% reduce to variables of interest
dat = dat(:, {'ID' 'age' 'sex' 'CMR' 'IMT_predicted' 'IMT_observed'});
writetable(dat, 'cmr/data_withCMR.csv')

%% make figures
% figure; corrplot([dat.CMR dat.waist_circumference dat.HDL dat.triglycerides dat.glucose dat.SBP], ...
%     'rows', 'pairwise', 'varNames', {'CMR' 'WC' 'HDL' 'TRIG' 'GLUC' 'SBP'})
% 
 figure; corrplot([dat.CMR dat.IMT_predicted dat.IMT_observed], 'rows', 'pairwise', ...
     'varNames', {'CMR' 'IMTpredicted' 'IMTobserved'})
 save('cmr/CMR_predicted_observed_correlations.png');

%% standardize variables if you want to report standardized betas via matlab

% these scatterplots are not longer 'real' though
datz = varfun(@zscore, dat);
datz.Properties.VariableNames = dat.Properties.VariableNames;

%% analyses of confounders

% Does the IAPS-related prediction relate to real IMT above and beyond CMR?
lm.CMR = fitlm(datz, 'IMT_observed ~ CMR')
writetable(lm.CMR.Coefficients, 'regression_IMT~CMR.csv', 'WriteRowNames', 1)

lm.CMR_IMTpredicted = fitlm(datz, 'IMT_observed ~ CMR + IMT_predicted')
writetable(lm.CMR_IMTpredicted.Coefficients, 'cmr/regression_IMT~CMR+brain.csv', 'WriteRowNames', 1)

lm.CMR_IMTpredicted_chgRsquared = lm.CMR_IMTpredicted.Rsquared.Ordinary - lm.CMR.Rsquared.Ordinary

% add age and sex

lm.CMR_age_sex = fitlm(datz, 'IMT_observed ~ CMR + age + sex')
writetable(lm.CMR_age_sex.Coefficients, 'cmr/regression_IMT~CMR+age+sex.csv', 'WriteRowNames', 1)

lm.CMR_age_sex_IMTpredicted = fitlm(datz, 'IMT_observed ~ CMR + age + sex + IMT_predicted')
writetable(lm.CMR_age_sex_IMTpredicted.Coefficients, 'cmr/regression_IMT~CMR+age+sex+brain.csv', 'WriteRowNames', 1)

lm.CMR_age_sex_IMTpredicted_chgRsquared = lm.CMR_age_sex_IMTpredicted.Rsquared.Ordinary - lm.CMR_age_sex.Rsquared.Ordinary

save('cmr/regression_CMR.mat', 'lm')

%% make scatterplots that show modifiers

%% center continuous predictors
%dat.IMT_predicted = scale(dat.IMT_predicted, 1);
%dat.age = scale(dat.age, 1);
%dat.CMR = scale(dat.CMR, 1);

%% fix sex
sex = {'Male' 'Female'}; % check this
dat.sex = sex(dat.sex)';

%% sex
lm_sex = fitlm(datz, 'IMT_observed ~ IMT_predicted*sex')
writetable(lm_sex.Coefficients, 'moderator_regression_sex.csv');
[r1, p1] = corr(table2array(dat(strcmp(dat.sex, 'Male'), 'IMT_predicted')), table2array(dat(strcmp(dat.sex, 'Male'), 'IMT_observed')))
[r2, p2] = corr(table2array(dat(strcmp(dat.sex, 'Female'), 'IMT_predicted')), table2array(dat(strcmp(dat.sex, 'Female'), 'IMT_observed')))

figure; gscatter(dat.IMT_predicted, dat.IMT_observed, dat.sex, 'br', 'oo'); lsline
xlabel('predicted'); ylabel('observed');
title(sprintf('sex x brain prediction interaction \n beta(se) = %.2f(%.2f) p = %.3f', ...
    table2array(lm_sex.Coefficients(end, 1)), table2array(lm_sex.Coefficients(end, 2)), table2array(lm_sex.Coefficients(end, end))));
g = gcf; 
g.Children(1).String{1} = sprintf('%s    r = %.2f p = %.3f', sex{1}, r1, p1);
g.Children(1).String{2} = sprintf('%s    r = %.2f p = %.3f', sex{2}, r2, p2);
g.Children(1).String(3:4) = [];
saveas(gcf, 'moderators/moderator_plot_sex.png');


%% age (median split)
lm_age = fitlm(datz, 'IMT_observed ~ IMT_predicted*age')
writetable(lm_age.Coefficients, 'moderator_regression_age.csv');
dat.age_medsplit = dat.age > median(dat.age);
ages = {'under median age' 'over median age'};
dat.age_medsplit = ages(dat.age_medsplit+1)';
[r1, p1] = corr(table2array(dat(dat.age < median(dat.age), 'IMT_predicted')), table2array(dat(dat.age < median(dat.age), 'IMT_observed')))
[r2, p2] = corr(table2array(dat(dat.age > median(dat.age), 'IMT_predicted')), table2array(dat(dat.age > median(dat.age), 'IMT_observed')))

figure; gscatter(dat.IMT_predicted, dat.IMT_observed, dat.age_medsplit, 'br', 'oo'); lsline
xlabel('predicted'); ylabel('observed');
title(sprintf('age x brain prediction interaction \n beta(se) = %.2f(%.2f) p = %.3f', ...
    table2array(lm_age.Coefficients(end, 1)), table2array(lm_age.Coefficients(end, 2)), table2array(lm_age.Coefficients(end, end))));
g = gcf; 
g.Children(1).String{1} = sprintf('%s    r = %.2f p = %.3f', ages{1}, r1, p1);
g.Children(1).String{2} = sprintf('%s    r = %.2f p = %.3f', ages{2}, r2, p2);
g.Children(1).String(3:4) = [];
saveas(gcf, 'moderators/moderator_plot_age.png');

%% CMR (interaction  - median split) 
lm_cmr = fitlm(datz, 'IMT_observed ~ IMT_predicted*CMR')
writetable(lm_cmr.Coefficients, 'moderator_regression_cmr.csv');
dat.CMR_medsplit = dat.CMR > median(dat.CMR);
cmr = {'under median CMR' 'over median CMR'};
dat.CMR_medsplit = cmr(dat.CMR_medsplit+1)';
[r1, p1] = corr(table2array(dat(dat.CMR < median(dat.CMR), 'IMT_predicted')), table2array(dat(dat.CMR < median(dat.CMR), 'IMT_observed')))
[r2, p2] = corr(table2array(dat(dat.CMR > median(dat.CMR), 'IMT_predicted')), table2array(dat(dat.CMR > median(dat.CMR), 'IMT_observed')))

figure; gscatter(dat.IMT_predicted, dat.IMT_observed, dat.CMR_medsplit, 'br', 'oo'); lsline
xlabel('predicted'); ylabel('observed');
title(sprintf('CMR x brain prediction interaction \n beta(se) = %.2f(%.2f) p = %.3f', ...
    table2array(lm_cmr.Coefficients(end, 1)), table2array(lm_cmr.Coefficients(end, 2)), table2array(lm_cmr.Coefficients(end, end))));
g = gcf; 
g.Children(1).String{1} = sprintf('%s    r = %.2f p = %.3f', cmr{1}, r1, p1);
g.Children(1).String{2} = sprintf('%s    r = %.2f p = %.3f', cmr{2}, r2, p2);
g.Children(1).String(3:4) = [];
saveas(gcf, 'moderators/moderator_plot_cmr.png');

