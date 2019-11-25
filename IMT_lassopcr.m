function out = IMT_lassopcr(task, varargin)

% IMT prediction using Lasso-PCR algorithm
% can input IAPS, OldFaces, or NewFaces imaging data
% can also specify affects for Faces tasks

% dataset generated using IMT_predict_paper_make_dataset.m

if nargin < 1, task = 'IAPS'; end % or OF or NF

%% Defaults
do_affect = 0;
mask = 'gm';
stratify_by_study = 0;
bootstrap_weights = 0;
iterations = 10;
switch task
    case 'IAPS'
        savedir = task;
    otherwise
        savedir = ['Faces/' task '/AllFaces'];
end   
imagefiles = sprintf('filepath_%s', task);
stratify_javi = 0;

fprintf('\nTask: %s\n', task);

%% parse varargin
for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
                
            case {'Affect' 'affect'}
                do_affect = 1;
                affect = varargin{i + 1};
                savedir = ['Faces/' task '/' affect];
                fprintf('\n Affect: %s \n', affect');
                %varargin{i + 1} = [];
                
            case {'bootstrap'}
                bootstrap_weights = 1;
                
            case {'iterations'}
                iterations = varargin{i + 1};
                
            case {'amygdala'}
                mask = 'amygdala';
                
            case {'lesion amygdala' 'amygdala lesion'}
                mask = 'lesion amygdala';
                
            case {'outerloops' 'outerfolds'}
                % specify how many outer folds
                % set to 1 to essentially ignore outer folds
                outerfolds = varargin{i + 1};
                
            case 'stratify_by_study' 
                % for IAPS, stratify according to study
                stratify_by_study = 1;
                
            case 'stratify_javi'
                stratify_javi = 1;
                                               
        end
    end
end

%% CD to ProjectDrive
mywd = cdtodrive; cd('AHAB_II/ML_projects/IMT_LassoPCR/');

%% Load dataset table
%dat = readtable('datasets/all_data.csv');
dat = readtable('datasets/ID_IMT_filepaths.csv');

%% Set variables, task and condition of interest
%task = 'OF'; % OF, NF
%condition = 'AllFaces'; % AllFaces, Anger, Fear, Happy, Neutral
myvar = 'IMT'; % IMT
myvar_string = 'IMT';

%% reduce dataset to participants with data

dat = rmmissing(dat(:, {'ID' 'study' myvar imagefiles}));
fprintf('\nDataset has %d participants.\n', size(dat, 1));

%% IMT is y
y = table2array(dat(:, myvar));

%% fix beginning of image file paths
dat.(imagefiles) = strrep(dat.(imagefiles), '/Volumes/ProjectDrive', mywd);

%% If doing Faces AFFECTS, edit image file locations
if do_affect
    switch affect
        case {'Anger' 'anger'}
            dat.(imagefiles) = strrep(dat.(imagefiles), 'Faces/con_0001', 'Faces_Affects/con_0001');
        case {'Fear' 'fear'}
            dat.(imagefiles) = strrep(dat.(imagefiles), 'Faces/con_0001', 'Faces_Affects/con_0002');
        case {'Happy' 'happy'}
            dat.(imagefiles) = strrep(dat.(imagefiles), 'Faces/con_0001', 'Faces_Affects/con_0003');
        case {'Neutral' 'neutral'}
            dat.(imagefiles) = strrep(dat.(imagefiles), 'Faces/con_0001', 'Faces_Affects/con_0004');
        case {'AngerFear'}
            dat.(imagefiles) = strrep(dat.(imagefiles), 'Faces/con_0001', 'Faces_Affects/con_0015');
    end
end

%% define mask
switch mask
    case 'amygdala'
        mask = [mywd '/AHAB_II/ML_projects/IMT_LassoPCR/masks/rAmygdala.nii'];  
        savedir = [savedir '/amygdala'];
    case 'lesion amygdala'
        mask = [mywd '/AHAB_II/ML_projects/IMT_LassoPCR/masks/amygdala_lesion.nii']; 
        savedir = [savedir '/amygdala_lesion'];
    otherwise
        mask = [mywd '/AHAB_II/ML_projects/IMT_LassoPCR/masks/resliced_grey25grey25.nii'];
end

%% load all imaging data
dat_fmri = fmri_data(dat.(imagefiles), mask);
dat_fmri.Y = y;
x = dat_fmri.dat';
    
%% LASSO-PCR on all data and save overall weightmap
if stratify_javi
    stats_train_all = lassopcr_cv(dat_fmri, [], 'stratify_javi');
else
    stats_train_all = lassopcr_cv(dat_fmri);
end
saveas(gcf, sprintf('%s/train_overall_predicted_observed.png', savedir));
weightmap = dat_fmri; weightmap.dat = stats_train_all.weight_obj.dat;
write(weightmap, 'fname', sprintf('%s/weightmap_unthresholded.nii', savedir));

%% Split the sample into train and test sets
if stratify_by_study
    % this is for IAPS since 2 studies have data
    % for the paper we will use train=AHAB but want to generate results for the reverse
    clear tmp; tmp.NumTestSets = 2;
    tmp.trIdx{1} = strcmp(dat.study, 'AHAB');
    tmp.teIdx{1} = strcmp(dat.study, 'PIP');
    tmp.trIdx{2} = strcmp(dat.study, 'PIP');
    tmp.teIdx{2} = strcmp(dat.study, 'AHAB');
elseif stratify_javi
    disp('Stratifying using discretized values')
    tmp = stratified_holdout_set(discretize(table2array(dat(:, myvar)), prctile(table2array(dat(:, myvar)), 0:20:100)), 5);
else
    tmp = stratified_holdout_set(table2array(dat(:, myvar)), 5);
    %dat.istest = tmp.teIdx;
    %dat.istrain = tmp.trIdx;
end

%% start outer loop

for i = 1:tmp.NumTestSets
    
    fprintf('\n\n Starting outer fold %d of %d\n\n', i, tmp.NumTestSets);
    
    %% assign x & y to train & test
        
    xtrain = x(tmp.trIdx{i}, :);
    ytrain = y(tmp.trIdx{i});
    xtest = x(tmp.teIdx{i}, :);
    ytest = y(tmp.teIdx{i});
    
    %% look at distributions
%     figure; subplot(3,1,1);
%     histogram(y);
%     title(sprintf('%s (all subjects) N = %d', myvar_string, length(y)));
%     g1 = gca;
%     subplot(3,1,2); histogram(ytrain);
%     title(sprintf('%s (training set) N = %d', myvar_string, length(ytrain)));
%     set(gca, 'XLim', g1.XLim);
%     subplot(3,1,3); histogram(ytest);
%     title(sprintf('%s (test set) N = %d', myvar_string, length(ytest)));
%     set(gca, 'XLim', g1.XLim);
%     saveas(gcf, sprintf('%s/IMT_distributions.png', savedir));    


    %% Train using in-house code
    if stratify_javi
        stats_train(i) = lassopcr_cv(xtrain, ytrain, 'stratify_javi');
    else
        stats_train(i) = lassopcr_cv(xtrain, ytrain);
    end
    %saveas(gcf, sprintf('%s/train_predicted_observed_%d.png', savedir, i));

    %% save unthresholded weightmap
    %weightmap = dat_fmri;
    %weightmap.dat = stats_train(i).weight_obj.dat;
    %write(weightmap, 'fname', sprintf('%s/weightmap_unthresholded_%d.nii', savedir, i))
    
    %% predict holdout data using each method
    % multiply multivariate weight map with participant's own response
    % add constant derived from model building

    stats_test(i).id = dat.ID(tmp.teIdx{i});
    stats_test(i).Y = ytest;
    stats_test(i).yfit = stats_train(i).constant + (xtest * stats_train(i).weight_obj.dat);
    
    [stats_test(i).test_corr, stats_test(i).test_pval] = corr(stats_test(i).Y, stats_test(i).yfit);
    stats_test(i).test_corr_ci = bootci(10000, @corr, stats_test(i).Y, stats_test(i).yfit)';    
    
    %% plot prediction in test
   
%     figure; scatter(stats_test(i).yfit, stats_test(i).Y);
%     xlabel('predicted Y in test'); 
%     ylabel('observed Y in test'); lsline;
%     title(sprintf('r(predicted,observed) = %.2f [%.2f %.2f], p = %.3f', ...
%     stats_test(i).test_corr, stats_test(i).test_corr_ci(1), stats_test(i).test_corr_ci(2), stats_test(i).test_pval));
    
end

%% make a single figure of outer fold data

%close all; 
figure;
fsize = 8;
for i = 1:tmp.NumTestSets
    % histograms
    subplot(tmp.NumTestSets, 4, ((i-1)*4) + 1)
    histogram(stats_train(i).Y); hold on;
    histogram(stats_test(i).Y); xlabel('IMT (mm)', 'FontSize', fsize); ylabel('Count', 'FontSize', fsize); %legend('Train', 'Test');
    % lasso stats
    subplot(tmp.NumTestSets, 4, ((i-1)*4) + 2)
    yyaxis left; plot(stats_train(i).lambdavals, stats_train(i).mse); ylabel('MSE', 'FontSize', fsize);
    yyaxis right; plot(stats_train(i).lambdavals, stats_train(i).num_components); ylabel('# Components', 'FontSize', fsize);
    xlabel('Lambda', 'FontSize', fsize); line([stats_train(i).best_lambda stats_train(i).best_lambda], [0 max(stats_train(i).num_components)], 'Color', 'r');
    title(sprintf('optimal lambda = %.4f\n%d components', stats_train(i).best_lambda, stats_train(i).num_components_from_best_lambda), ...
        'FontSize', fsize)
    % scatter train
    subplot(tmp.NumTestSets, 4, ((i-1)*4) + 3)
    scatter(stats_train(i).yfit, stats_train(i).Y, 6, 'filled', 'black'); 
    l = lsline; l.Color = 'blue';
    xlabel('predicted', 'FontSize', fsize); ylabel('observed', 'FontSize', fsize);
    title(sprintf('r = %.2f', stats_train(i).pred_outcome_r));
    text(100, 100, 'corr');
    % scatter test
    subplot(tmp.NumTestSets, 4, ((i-1)*4) + 4)
    scatter(stats_test(i).yfit, stats_test(i).Y, 6, 'filled', 'black'); 
    l = lsline; l.Color = 'blue';
    [stats_test(i).test_corr, stats_test(i).test_pval] = corr(stats_test(i).Y, stats_test(i).yfit);
    stats_test(i).test_corr_ci = bootci(10000, @corr, stats_test(i).Y, stats_test(i).yfit)';
    xlabel('predicted', 'FontSize', fsize); ylabel('observed', 'FontSize', fsize); 
    title(sprintf('r = %.2f', stats_test(i).test_corr));
end
set(gcf, 'Position', [10 10 1300 900])
saveas(gcf, sprintf('%s/nested_xval_output.png', savedir));
saveas(gcf, sprintf('%s/nested_xval_output.fig', savedir));
       

%% combine outputs

out.mask = mask;
out.cv_folds = tmp;
out.model_all_data = stats_train_all;
out.inner_training_folds = [stats_train];
out.inner_test_folds = [stats_test];
%out.best_lambda = [stats_train.best_lambda];
%out.num_components_from_best_lambda = [stats_train.num_components_from_best_lambda];
out.ID = vertcat(stats_test.id);
out.Y = vertcat(stats_test.Y); out.yfit = vertcat(stats_test.yfit);
[out.pred_outcome_r, out.pred_outcome_pval] = corr(out.Y, out.yfit);
out.pred_outcome_r_ci = bootci(10000, @corr, out.yfit, out.Y);
out.mean_weightmap = mean(struct2array([out.inner_training_folds.weight_obj]), 2)

%% write weightmaps for inner folds
w = dat_fmri; w.dat = struct2array([out.inner_training_folds.weight_obj]);
write(w, 'fname', sprintf('%s/weightmap_innerfolds_unthresholded.nii', savedir));

%% scatterplot for overall holdout performance
figure; scatter(out.yfit, out.Y, 6, 'filled', 'black'); l = lsline; l.Color = 'blue';
title(sprintf('r = %.2f [%.2f %.2f] p = %.3f', ...
    out.pred_outcome_r, out.pred_outcome_r_ci(1), out.pred_outcome_r_ci(2), out.pred_outcome_pval));
xlabel('predicted IMT (mm)'); ylabel('observed IMT (mm)');
saveas(gcf, sprintf('%s/nested_xval_overall_predicted_observed.png', savedir));
saveas(gcf, sprintf('%s/nested_xval_overall_predicted_observed.fig', savedir));

%% save files for later
save(sprintf('%s/data.mat', savedir), ...
    'stats_train', 'stats_test', 'dat', 'out');
%%
return;

stats_test.Y = y_test;
stats_test.yfit = stats_train.constant + (fmri_test.dat' * stats_train.weight_obj.dat);

[stats_test.test_corr, stats_test.test_pval] = corr(stats_test.Y, stats_test.yfit);
stats_test.test_corr_ci = bootci(10000, @corr, stats_test.Y, stats_test.yfit)';
test_err = stats_test.yfit - stats_test.Y;
stats_test.mae = nanmean(abs(test_err));
stats_test.mse = mean(test_err.^2);
figure; scatter(stats_test.yfit, stats_test.Y);
xlabel('predicted Y in test'); 
ylabel('observed Y in test'); lsline;
title(sprintf('r(predicted,observed) = %.2f [%.2f %.2f], p = %.3f', ...
    stats_test.test_corr, stats_test.test_corr_ci(1), stats_test.test_corr_ci(2), stats_test.test_pval));
saveas(gcf, sprintf('%s/test_predicted_observed.png', savedir));

