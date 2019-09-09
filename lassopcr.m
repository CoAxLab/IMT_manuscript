function out = lassopcr(x, y, varargin)
% LASSO-PCR given multivariate x and univariate y
% define lambda, or leave blank to run across a range of lambdas
% This code has been heavily adapted from Tor Wager's predict.m function
% https://github.com/canlab/CanlabCore
% See Wager et al. (2011) J Neurosci for a description of the procedure
% And cite Tor's papers if you use this code!

% defaults
definelambda = 0; % when this is 1, enter one or more lambdas to test
num_components = length(y);
showfigs = 0;
savebootweights = 0;
pcr = 1;

% need to add covariate options!

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            case {'lambda'}
                % follow this by one or more lambdas to test
                definelambda = 1;
                lambdavals = varargin{i+1};
                varargin{i} = [];
                varargin{i+1} = []; 
            case {'num_components'} 
                % in case you want to trim off components 
                % that explain little variance in x
                num_components = varargin{i+1};
                varargin{i} = [];
                varargin{i+1} = []; 
            case {'showfigures', 'figures'}
                showfigs = 1;
                varargin{i} = [];
            case {'savebootweights'}
                savebootweights = 1;
                showfigs = 0;
                varargin{i} = [];
            case 'nopcr'
                pcr = 0;
                varargin{i} = [];
        end
    end
end

if isa(x, 'fmri_data')
    y = x.Y;
    x = x.dat'; % transposed from default loading of fmri_data
end

%% Dimensionality Reduction
if pcr
    [pc, ~, ~] = svd(scale(x, 1)', 'econ'); % replace princomp with SVD on transpose to reduce running time. 
                                     % specifically, SVD on the p x n matrix of X 
                                     % (not n x p, so might need to transpose)
    pc(:, num_components:end) = [];                % remove the last component, which is close to zero
                                   % edit:replaced 'pc(:,size(xtrain,1)) = [];' with
                                   % end to accomodate predictor matrices with
                                   % fewer features (voxels) than trials. SG     
    sc = x * pc;


    % further trim components if there are redundant components (e.g., when bootstrapping)
    if rank(sc) == size(sc,2)
        numcomps = rank(sc); 
    elseif rank(sc) < size(sc,2)
        numcomps = rank(sc)-1;
    end
else
    sc = x;
end

%% Fit LASSO
if definelambda
    [B, stats] = lassoglm(sc(:, 1:numcomps), y, 'normal', 'Lambda', lambdavals); 
else
    [B, stats] = lassoglm(sc(:, 1:numcomps), y, 'normal');  
end

%% Refit using OLS (???)

%% Project to x space
if pcr
    out = pc(:, 1:numcomps) * B; 
else
    out = B;
end

if savebootweights
    out = out';
end


