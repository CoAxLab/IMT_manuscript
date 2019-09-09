function out = IMT_generalizability(imgs, pattern)

% run post-hoc tests of model generalizability across tasks
% e.g., to see if the IAPS pattern can predict IMT using OldFaces images:
%       imgs = 'OF'
%       pattern = 'IAPS'
% can also use Luke Chang's PINES or Pete's stressor-evoked SBP maps as patterns

%%
cdtodrive; cd AHAB_II/ML_projects/IMT_LassoPCR

%%

switch imgs
    case 'IAPS'
        imgset = load('IAPS/data.mat');
    case 'OF'
        imgset = load('Faces/OF/AllFaces/data.mat');  
    case 'NF'
        imgset = load('Faces/NF/AllFaces/data.mat');
end

switch pattern
    case 'IAPS'
        patternset = load('IAPS/data.mat');
    case 'OF'
        patternset = load('Faces/OF/AllFaces/data.mat');  
    case 'NF'
        patternset = load('Faces/NF/AllFaces/data.mat');
    case 'PINES'
        pines = fmri_data(which('Rating_Weights_LOSO_2.nii'), imgset.out.mask); 
        patternset.out.model_all_data.constant = 0;
        patternset.out.model_all_data.weight_obj.dat = pines.dat;
    case 'stressSBP'
        sbp = fmri_data('masks/avgChgSBP_onAvgTaskHE_weight.nii', imgset.out.mask); 
        patternset.out.model_all_data.constant = 0;
        patternset.out.model_all_data.weight_obj.dat = sbp.dat; 
end

%% apply pattern to contrasts and compare to observed IMT

imagefiles = table2cell(imgset.dat(:, contains(imgset.dat.Properties.VariableNames, 'filepath')));
images = fmri_data(imagefiles, imgset.out.mask);
out.yfit = patternset.out.model_all_data.constant + (images.dat' * patternset.out.model_all_data.weight_obj.dat);
[r, pval] = corr(out.yfit, imgset.out.Y);
figure; scatter(out.yfit, imgset.out.Y);
xlabel(sprintf('Predicted IMT using %s contrast maps and %s pattern', imgs, pattern));
ylabel('Observed IMT'); title(sprintf('r = %.2f, p = %.3f', r, pval)); lsline;
saveas(gcf, sprintf('Generalizability/predicted_observed_pattern-%s_contrasts-%s.png', pattern, imgs));
% write data to csv
out.Y = imgset.out.Y;
writetable(struct2table(out), sprintf('Generalizability/predicted_observed_pattern-%s_contrasts-%s.csv', pattern, imgs));

