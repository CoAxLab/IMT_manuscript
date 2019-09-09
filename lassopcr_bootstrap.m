function out = lassopcr_bootstrap(iterations, threshold_method)

% bootstrap LASSO-PCR procedure on a model
% run from folder containing 'data.mat' output structure from IMT_lassopcr.m 

if nargin < 2, threshold_method = 'parametric'; end
if nargin < 1, iterations = 100; end
do_bootstrap = 0;
do_permutations = 0;

t = tic;

%% open parallel pool
myPool = parpool;
opt = statset('UseParallel', true);

%% load data and images

load('data.mat')
imagefiles = table2cell(dat(:, contains(dat.Properties.VariableNames, 'filepath')));
fmri_train = fmri_data(imagefiles, out.mask);
fmri_train.Y = out.Y;

%% run LASSO-PCR once on the entire dataset (necessary if nested cross-validation was done prior)
stats_train = lassopcr_cv(fmri_train, [], 'noplots');

%% BOOTSTRAPPING and repeat the entire LASSO-PCR procedure
% run bootstrapping in small batches since it sometimes fails

t = tic; 
boot_all = [];
batchsize = 50; % how many bootstraps to run at once?

fprintf('\n\nBOOTSTRAPPING VOXELS in batches of %d\n\n', batchsize)
bootfun = @(x, y) lassopcr_cv(x, y, 'savebootweights');

while size(boot_all, 1) < iterations
    fprintf(1, '\n%d bootstraps completed\n\n', size(boot_all, 1));
    rng 'shuffle' % shuffle random number generator
    try
        clear boot;
        boot =  bootstrp(batchsize, bootfun, fmri_train.dat', fmri_train.Y, 'Options', opt);
        %boot = boot(any(boot, 2), :); % remove rows with all zeros
        boot_all = [boot_all; boot];
    catch exception
        fprintf(1, '/nError during bootstrap routine - starting new batch /n');
    end
end

out.boot_iterations = iterations;
out.seconds_per_iteration = toc(t)/iterations;

% 17sec per iteration on linux machine (without parallel pool)
% 6 sec with parallel pool

%% combine output
boot_all = boot_all(1:iterations, :);

%% optionally save EVERY bootstrap draw (this creates large files)
out.boot_all = boot_all; 

%% calculate voxelwise statistics from boot_all
out.threshold_method = threshold_method;
out.original_map = stats_train.weight_obj.dat;
out.boot_mean = nanmean(boot_all)';
out.boot_std = nanstd(boot_all)';
out.boot_z = (nanmean(boot_all)./nanstd(boot_all))'; 

%% construct p-values (parametric and nonparametric versions)
switch threshold_method
    case 'parametric'       
        out.pval = 2 * (1 - normcdf(abs(out.boot_z)));
    case 'nonparametric'
        for i = 1:length(out.boot_mean)
            if out.boot_mean(i) > 0
                out.pval(i) = (1 - sum(out.boot_all(:, i) > 0)/out.boot_iterations)*2;
            elseif out.boot_mean(i) < 0
                out.pval(i) = (1 - sum(out.boot_all(:, i) < 0)/out.boot_iterations)*2;
            else
                out.pval(i) = NaN;
            end
        end
        out.pval(out.pval > 1) = .999;
        out.pval = out.pval';
end

%% make figures, correlate original weightmap with mean / Z of bootstrapped/permuted

mkdir(sprintf('bootstrap/%diterations', iterations));

figure; subplot(2,1,1); hist(out.boot_z); xlabel('voxelwise Z values');
subplot(2,1,2); hist(out.pval); xlabel('voxelwise p values');
saveas(gcf, sprintf('bootstrap/%diterations/histogram_z_p.png', iterations));

for i = 1:size(boot_all, 1)
    out.corr_orig_with_each_boot(i) = corr(out.original_map, boot_all(i, :)', 'rows', 'pairwise');
end
figure; hist(out.corr_orig_with_each_boot); title('Distribution of observed-bootstrapped voxelwise correlations');
saveas(gcf, sprintf('bootstrap/%diterations/histogram_original_bootstrap_comparison.png', iterations));

figure; scatter(out.original_map, out.boot_mean);
xlabel('Observed Pattern Map'); ylabel('Bootstrap Mean');
title(sprintf('r = %.2f', corr(out.original_map, out.boot_mean, 'rows', 'pairwise')));
saveas(gcf, sprintf('bootstrap/%diterations/scatter_weightmap_original_bootstrap-mean_comparison.png', iterations));

figure; scatter(out.original_map, out.boot_z);
xlabel('Observed Pattern Map'); ylabel('Bootstrap Z');
title(sprintf('r = %.2f', corr(out.original_map, out.boot_z, 'rows', 'pairwise')));
saveas(gcf, sprintf('bootstrap/%diterations/scatter_weightmap_original_bootstrap-z_comparison.png', iterations));

%% mask weightmap with bootstrap output
% generate table on significant clusters

out.bootmap = statistic_image('dat', stats_train.weight_obj.dat, ...
    'volInfo', fmri_train.volInfo, ...
    'p', out.pval, ...
    'removed_voxels', fmri_train.removed_voxels);

write(out.bootmap, 'fname', sprintf('bootstrap/%diterations/weightmap_bootstrap_unthresholded.nii', iterations));

%% save thresholded maps
alphas = [.001 .01 .05]; k = [20 50]; % try many alphas and cluster sizes

switch threshold_method
    case 'parametric'
        thresh_string = '';
    case 'nonparametric'
        thresh_string = '_percentile';
end

for aa = 1:length(alphas)
    for kk = 1:length(k)
        alpha_string = num2str(alphas(aa));  
        bootmap_thresh{aa, kk} = threshold(out.bootmap, alphas(aa), 'unc', 'k', k(kk));
        r = region(bootmap_thresh{aa, kk});
        if max([r.numVox]) > 0 % proceed only if clusters survived threshold
            write(bootmap_thresh{aa, kk}, 'thresh', 'fname', sprintf('bootstrap/%diterations/weightmap_bootstrap_%s_p%s_k%d.nii', iterations, thresh_string, alpha_string(3:end), k(kk)));
            [pos, neg, results_table] = table(r);
            writetable(results_table, sprintf('bootstrap/%diterations/weightmap_bootstrap%s_p%s_k%d_clusters.csv', iterations, thresh_string, alpha_string(3:end), k(kk)));
        end
    end
end

%% save entire data file
save(sprintf('bootstrap/%diterations/data.mat', iterations), 'out');

%% display p < .05 k > 50 map
orthviews(bootmap_thresh{3, 2})

delete(myPool)

    
