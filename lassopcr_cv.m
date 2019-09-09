function out = lassopcr_cv(x, y, varargin)
% LASSO-PCR with cross-validation to optimize lambda
% This code has been heavily adapted from Tor Wager's predict.m function
% https://github.com/canlab/CanlabCore
% See Wager et al. (2011) J Neurosci for a description of the procedure
% And cite Tor's papers if you use this code!
% 
% Method 1: 
%   To optimize lambda in training set, do entire LASSO-PCR scheme on each fold 
%   (i.e., PCA data-to-components -> LASSO -> project components-to-data -> fit pattern on holdout)
%   This requires us to jerry-rig the cross-validation scheme
%   For each fold:
%       1. Run SVD on xtrain, 
%       2. Try multiple lambdas in the LASSO, 
%       3. Back-project coefficients from each fit to x (voxel) space and fit on ytrain
%       4. Calculate error for each observation once for each lambda
%       5. Calculate overall (binned) mse per each lambda
% Method 2: 
%   To optimize lambda in training set, do PCA on entire training set, 
%   and optimize lambda within this component space
%   This can be completed in one line within the lassoglm machinery
%   (but it doesn't save out the cross-validated predictions, only the error)
% Need to add options for covariates (Powell et al 2018 Network Neuroscience)

%% Instead of entering x and y separately, you can enter them together as a fmri_data object, where:
%   x is a CANlab fmri_data structure;
%   x.Y is the outcome
%   x.dat is the imaging data
% NOTE if you want to use varargin flags with this option, put a [] placeholder in the place of y (field 2)
% NOTE also that we transpose the fmridat.dat matrix around a bit within this script

if isa(x, 'fmri_data')
    y = x.Y;
    x = x.dat'; % transposed from default loading of fmri_data
end

%% defaults
nlambdas = 1000;
nfolds = 5;
nouterfolds = 5;
method = 1;
do_ols = 0;
savebootweights = 0;
makeplots = 1;
eval_metric = 'mse';
components_retained = 1;
center_y = 0;
do_nested_xval = 0;
stratify_javi = 0;

%% parse varargin
for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
                
            case 'Method2'
                method = 2;
                
            case 'noplots'
                makeplots = 0;
                
            case {'nlambdas', 'lambdas', 'num_lambdas'}
                nlambdas = varargin{i + 1};
                
            case {'folds', 'nfolds'}
                nfolds = varargin{i + 1};
                           
            case {'savebootweights'}
                savebootweights = 1;
                makeplots = 0;
                
            case {'do_ols'}
                do_ols = 1;
                
            case {'eval', 'eval_metric'} 
                % metric on which to evaluate holdout data to optimize lambda
                eval_metric = varargin{i + 1};
                
            case {'components_retained'} 
                % only keep a proportion [0:1] of components following PCA (reduce overfitting)
                components_retained = varargin{i + 1}; 
                
            case {'center_y'}
                center_y = 1;       
                
            case {'nested' 'do_nested_xval'}
                % stratifies x and y into k folds
                % recursively calls lassopcr_cv within each fold
                % saves holdout test data
                % concatenates test data at end
                do_nested_xval = 1;
                varargin{i} = [];
            
            case {'outerloops' 'outerfolds'}
                % specify how many outer folds
                % set to 1 to essentially ignore outer folds
                outerfolds = varargin{i + 1};
                
            case {'stratify_javi'}
                stratify_javi = 1;
                

        end
    end
end

%clear out cv
out.nfolds = nfolds; 
out.nlambdas = nlambdas;
out.method = method;

if center_y; y = scale(y, 1); end


%% nested cross validation - built into this code as an optional method
% this calls the actual LASSO-PCR code below in a loop and combines output
if do_nested_xval
    fprintf(1, '\nStarting nested cross-validation procedure\n');
    if stratify_javi
        disp('stratifying using discretized values')
        tmp = stratified_holdout_set(discretize(y, prctile(y, 0:20:100)), nouterfolds);
    else
        tmp = stratified_holdout_set(y, nouterfolds);
    end
    for i = 1:tmp.NumTestSets
        fprintf('\n\nStarting outer fold %d of %d\n\n', i, tmp.NumTestSets); 
        %% assign x & y to train & tesst      
        xtrain = x(tmp.trIdx{i}, :);
        ytrain = y(tmp.trIdx{i});
        xtest = x(tmp.teIdx{i}, :);
        ytest = y(tmp.teIdx{i});
        
        stats_train(i) = lassopcr_cv(xtrain, ytrain, varargin{:});
        stats_test(i).Y = ytest;
        stats_test(i).yfit = stats_train(i).constant + (xtest * stats_train(i).weight_obj.dat);

        [stats_test(i).test_corr, stats_test(i).test_pval] = corr(stats_test(i).Y, stats_test(i).yfit);
        stats_test(i).test_corr_ci = bootci(10000, @corr, stats_test(i).Y, stats_test(i).yfit)';   
    end
    %% combine outputs
    out.inner_training_folds = [stats_train];
    out.inner_test_folds = [stats_test];
    out.best_labmda = [stats_train.best_lambda];
    out.num_components_from_best_lambda = [stats_train.num_components_from_best_lambda];
    out.Y = vertcat(stats_test.Y); out.yfit = vertcat(stats_test.yfit);
    [out.pred_outcome_r, out.pred_outcome_pval] = corr(out.Y, out.yfit);
    out.pred_outcome_r_ci = bootci(10000, @corr, out.yfit, out.Y);
    out.mean_weightmap = mean(struct2array([out.inner_training_folds.weight_obj]), 2)
    return
end

%% ACTUAL LASSO-PCR CODE %%

%% first generate a range of lambdas by fitting the entire training set
% this gives a 'max' lambda for which all predictors are zero
% this depends on how many predictors were left by the SVD, 
% also depends on range of Y,
% so it is okay to run on the entire training set
% we can squeeze more lambda values between 0 and the max lambda

try
    [pc,~,~] = svd(scale(x,1)', 'econ');
catch exception % sometimes the SVD does not converge (especially during bootstrapping
    sprintf(['Error! ' exception.message '\nSkipping lassoPCR'])
    if savebootweights
        out = zeros(1, length(x));
    end
    return;
end

% option to retain a proportion of components following PCA (reduce overfitting)
%num_components = round(components_retained*size(pc, 2));
pc(:, end) = [];               

sc = x * pc; 

% further trim components if x (sc) is rank deficient (e.g., when bootstrapping)
if rank(sc) == size(sc,2)
    numcomps = rank(sc); 
elseif rank(sc) < size(sc,2)
    numcomps = rank(sc)-1;
end

[B, stats] = lassoglm(sc(:, 1:numcomps), y, 'normal');  

lambdavals = linspace(0, max(stats.Lambda), nlambdas)';
out.lambdavals = lambdavals;
                                   
[B, xstats] = lassoglm(sc(:, 1:numcomps), y, 'normal', 'Lambda', lambdavals); 
out.num_components = xstats.DF;

%% k-fold cross-validation to find optimal lambda
        
switch method
    case 1
        
        if stratify_javi
            disp('stratifying using discretized values')
            out.cv_folds = stratified_holdout_set(discretize(y, prctile(y, 0:20:100)), nfolds);
        else
            out.cv_folds = stratified_holdout_set(y, nfolds);
        end
                  
        for i = 1:out.cv_folds.NumTestSets
            xtrain = x(out.cv_folds.trIdx{i}, :);
            ytrain = y(out.cv_folds.trIdx{i});
            xtest = x(out.cv_folds.teIdx{i}, :);
            ytest = y(out.cv_folds.teIdx{i});
            
            try
                [pc,~,~] = svd(scale(xtrain,1)', 'econ');
            catch exception % sometimes the SVD does not converge (especially during bootstrapping
                sprintf(['Error! ' exception.message '\nSkipping lassoPCR'])
                if savebootweights
                    out = zeros(1, length(x));
                end
                return;
            end
            
            pc(:, end) = [];                   
            sc = xtrain * pc;
            
            % further trim components if x (sc) is rank deficient (e.g., when bootstrapping)
            if rank(sc) == size(sc,2)
                numcomps = rank(sc); 
            elseif rank(sc) < size(sc,2)
                numcomps = rank(sc)-1;
            end

            [B,stats] = lassoglm(sc(:, 1:numcomps), ytrain, 'normal', 'Lambda', lambdavals);   
           
%            if do_ols % transformation is different here from below since there are multiple sets of B
%                ols_betas = zeros(size(B));
%                 for b = 1:length(lambdavals)
%                     s = find(B(:, b)~=0);
%                     r = regress(ytrain, sc(:, s));
%                     for ii = 1:length(s)
%                         ols_betas(s(ii), b) = r(ii);
%                     end
%                 end
%                 B = ols_betas;
%            end
           
          % project to x (voxel) space
          vox_weights = pc(:, 1:numcomps) * B;  

           % fit to holdout data (in train) to calculate predicted y
           cv(i).ytest = ytest;
           cv(i).yfit = stats.Intercept + (xtest * vox_weights);
           %cv(i).yfit = (xtest * vox_weights);
           if ~ savebootweights
               fprintf(1, 'Done with fold %d of %d \n', i, out.cv_folds.NumTestSets)
           end
        end

    case 2
        
        out.cv_folds = stratified_holdout_set(y, nfolds);
        
        for i = 1:out.cv_folds.NumTestSets
            xtrain = sc(out.cv_folds.trIdx{i}, :);
            ytrain = y(out.cv_folds.trIdx{i});
            xtest = sc(out.cv_folds.teIdx{i}, :);
            ytest = y(out.cv_folds.teIdx{i});     
            
            [B,stats] = lassoglm(xtrain, ytrain, 'normal', 'Lambda', lambdavals);   
            
           cv(i).ytest = ytest;
           cv(i).yfit = stats.Intercept + (xtest * B);
           
        end
        
        % Alternative way of doing case 2 in fewer lines
        % drawback is that we cannot get cross-validated yfit values for each fold
        
        %[xB, xstats] = lassoglm(sc, y, 'normal', 'Lambda', lambdavals, 'CV', nfolds); 
        %out.best_lambda = xstats.LambdaMinDeviance;
        %out.mse = xstats.Deviance;
        
        
        
end


%% combine yfit and ytest(observed) from each fold

% calculate mse for each lambda
yfit = vertcat(cv.yfit);
ytest = vertcat(cv.ytest);
error = yfit - ytest;


%% find best lambda according to your 'eval' metric
% options are:
%   mse - mean squared error (default; find min)
%   medse - median square error (find min)
%   predicted-observed r (find max)

switch eval_metric % metric on which to evaluate holdout data to optimize lambda
    case 'mse' % mean square error
        out.evaluation_metric = 'minimum of mean squared error';
        out.mse = mean(error.^2)'; % also consider median squared error
        out.best_lambda_index = find(out.mse == min(out.mse));
        out.metric_best_lambda = out.mse(out.best_lambda_index);
        yvals = out.mse; % for plotting
    case 'medse' % median square error
        out.evaluation_metric = 'minimum of median squared error';
        out.medse = median(error.^2)';
        out.best_lambda_index = find(out.medse == min(out.medse));
        out.metric_best_lambda = out.medse(out.best_lambda_index);
        yvals = out.medse; % for plotting
    case 'mae' % mean absolute error
        out.evaluation_metric = 'minimum of mean absolute error';
        out.mae = mean(abs(error))'; % also consider median squared error
        out.best_lambda_index = find(out.mae == min(out.mae));
        out.metric_best_lambda = out.mae(out.best_lambda_index);
        yvals = out.mae; % for plotting  
    case 'r' % predicted-observed correlation
        out.evaluation_metric = 'maximum of predicted-observed correlation';
        for rr = 1:nlambdas
            out.r(rr) = corr(yfit(:, rr), ytest);
        end
        out.best_lambda_index = find(out.r == max(out.r));
        out.metric_best_lambda = out.r(out.best_lambda_index);
        yvals = out.r; % for plotting
end

%% very rarely the 'best lambda' was tested more than once and so length(out.best_lambda_index) > 1
if length(out.best_lambda_index) > 1
    out.best_lambda_index = out.best_lambda_index(1);
end

%% identify the best lambda
out.best_lambda = lambdavals(out.best_lambda_index);

% extract fitted (cross-validated) y for that best lambda
out.Y = ytest; 
out.yfit = yfit(:, out.best_lambda_index);

% report cross-validated predicted-observed correlation using that lambda
[out.pred_outcome_r, out.pred_outcome_pval] = corr(out.yfit, out.Y);
out.pred_outcome_r_ci = bootci(10000, @corr, out.yfit, out.Y);


%% refit on entire training set using the selected lambda
try
    [pc,~,~] = svd(scale(x,1)', 'econ');
catch exception % sometimes the SVD does not converge (especially during bootstrapping
    sprintf(['Error! ' exception.message '\nSkipping lassoPCR'])
    if savebootweights
        out = zeros(1, length(x));
    end
    return;
end

% option to retain a proportion of components following PCA (reduce overfitting)
num_components = round(components_retained*size(pc, 2));
pc(:, num_components:end) = [];    

sc = x * pc;

% further trim components if x (sc) is rank deficient (e.g., when bootstrapping)
if rank(sc) == size(sc,2)
    numcomps = rank(sc); 
elseif rank(sc) < size(sc,2)
    numcomps = rank(sc)-1;
end

[B,stats] = lassoglm(sc(:, 1:numcomps), y, 'normal', 'Lambda', out.best_lambda);
out.constant = stats.Intercept;
out.component_betas = B;
out.num_components_from_best_lambda = stats.DF;

%% OLS
if do_ols
    i = find(B~=0);
    r = regress(y, sc(:, i));
    betas = zeros(size(B));
    for ii = 1:length(i)
        betas(i(ii)) = r(ii);
    end
    B = betas;
end

%% make weight object from final model
out.weight_obj.dat = pc(:, 1:numcomps) * B;  

%% make plots
% Plot relationships of lambda, number of components, error
if makeplots
    figure; 
    subplot(1, 2, 1)
    yyaxis left; plot(lambdavals, yvals); ylabel(eval_metric);
    yyaxis right; plot(lambdavals, xstats.DF); ylabel('Number of Components');
    xlabel('Lambda'); line([out.best_lambda out.best_lambda], [0 max(xstats.DF)], 'Color', 'r');
    title(sprintf('optimal lambda = %.4f, \n %d components', out.best_lambda, out.num_components_from_best_lambda))
    subplot(1, 2, 2); scatter(out.yfit, out.Y); lsline;
    xlabel('predicted Y in train'); ylabel('observed Y in train');
    title(sprintf('r = %.2f [%.2f %.2f]', ... 
        out.pred_outcome_r, out.pred_outcome_r_ci(1), out.pred_outcome_r_ci(2)));
end

   
%% option to only save out weightmap (for bootstrapping)
if savebootweights
    out = out.weight_obj.dat';
    % Sometimes the weightmap is not the correct dimensions (unclear why this happens)
    % If this happens, make it the correct size and all ZEROS so it does not
    % disrupt bootstrapping process - can ignore later
    if length(out) ~= length(x)
        out = zeros(1, length(x));
        fprintf('\nBootstrapped map has incorrect dimensions...\n')
    end
end
