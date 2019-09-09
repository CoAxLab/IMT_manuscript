function out = lassopcr_cv_nested(x, y, varargin)

if isa(x, 'fmri_data')
    y = x.Y;
    x = x.dat'; % transposed from default loading of fmri_data
end

%% defaults
nlambdas = 1000;
nfolds = 5;
method = 1;
do_ols = 0;
savebootweights = 0;
savebootweights2 = 0;
makeplots = 1;
eval_metric = 'mse';
components_retained = 1;
center_y = 0;
do_nested_xval = 0;

%% parse varargin
% NOTE - see lassopcr_cv for other parameters
% this code will pass those parameters into lassopcr_cv
% the only parameter used in the outer loop is 'nfolds'

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
                
            case {'folds', 'nfolds'}
                nfolds = varargin{i + 1};
                           
        end
    end
end
    
%% add 'noplots' to varargin to suppress plotting from inner folds
%varargin = [varargin 'noplots'];
    
%% start outer loop
tmp = stratified_holdout_set(y, nfolds);

for i = 1:nfolds
    
    fprintf('\n\n Starting outer fold %d of %d\n\n', i, nfolds);
    
    xtrain = x(tmp.trIdx{i}, :);
    ytrain = y(tmp.trIdx{i});
    xtest = x(tmp.teIdx{i}, :);
    ytest = y(tmp.teIdx{i});
    
    stats_train(i) = lassopcr_cv(xtrain, ytrain, varargin{:});
    
    stats_test(i).Y = ytest;
    stats_test(i).yfit = stats_train(i).constant + (xtest * stats_train(i).weight_obj.dat);
    
    figure; scatter(stats_test(i).yfit, stats_test(i).Y);
    [stats_test(i).test_corr, stats_test(i).test_pval] = corr(stats_test(i).Y, stats_test(i).yfit);
    stats_test(i).test_corr_ci = bootci(10000, @corr, stats_test(i).Y, stats_test(i).yfit)';
    xlabel('predicted Y in test'); 
    ylabel('observed Y in test'); lsline;
    title(sprintf('r(predicted,observed) = %.2f [%.2f %.2f], p = %.3f', ...
    stats_test(i).test_corr, stats_test(i).test_corr_ci(1), stats_test(i).test_corr_ci(2), stats_test(i).test_pval));
    
end

%% combine outputs

out.inner_training_folds = [stats_train];
out.inner_test_folds = [stats_test];
out.best_labmda = [stats_train.best_lambda];
out.num_components_from_best_lambda = [stats_train.num_components_from_best_lambda];
out.Y = vertcat(stats_test.Y); out.yfit = vertcat(stats_test.yfit);
[out.pred_outcome_r, out.pred_outcome_pval] = corr(out.Y, out.yfit);
out.pred_outcome_r_ci = bootci(10000, @corr, out.yfit, out.Y);
figure; scatter(out.yfit, out.Y); lsline;
title(sprintf('r(predicted,observed) = %.2f [%.2f %.2f] p = %.3f', ...
    out.pred_outcome_r, out.pred_outcome_r_ci(1), out.pred_outcome_r_ci(2), out.pred_outcome_pval));



