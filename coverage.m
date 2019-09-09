wd = cdtodrive;

greymask = 'AHAB_II/ML_projects/IMT_LassoPCR/masks/resliced_grey25grey25.nii';
amygmask = 'AHAB_II/ML_projects/IMT_LassoPCR/masks/Amygdala.nii'


tasks = {'ER' 'OldFaces' 'Faces'}

for t = 1:3
    
    % list files
    ahab = filenames(sprintf('AHAB_II/SPM12/First_Level/*/%s/mask.nii', tasks{t}), 'absolute')
    pip = filenames(sprintf('PIP/SPM12/First_Level/*/%s/mask.nii', tasks{t}), 'absolute')

    all_masks = fmri_data([ahab; pip])

    % get IDs
    clear id; 
    for i = 1:size(all_masks.fullpath, 1)
        tmp = strsplit(all_masks.fullpath(i,:), '/');
        id(i) = str2double(tmp{end-2});
    end
    id = id';
    
    % get coverage in grey matter mask and amygdala mask
    tmp = extract_roi_averages(all_masks, greymask)
    gm_coverage = sum(tmp.all_data ~= 0, 2)/size(tmp.all_data, 2);
    
    tmp = extract_roi_averages(all_masks, amygmask)
    amyg_coverage = sum(tmp.all_data ~= 0, 2)/size(tmp.all_data, 2);
  
    tt = table(id, gm_coverage, amyg_coverage)
      
    writetable(tt, sprintf('AHAB_II/ML_projects/IMT_LassoPCR/datasets/coverage_%s.csv', tasks{t}));
    
end
%%
%     tmp = 
%         r = region(fmri_data(all_masks, greymask));
%         numVox(i) = r.numVox;
%         percCoverage(i) = round(numVox(i)/grey.numVox*100, 0);
%        % o2 = addblobs(o2, r, 'trans', 'onecolor')
%        % hText = text(0,200, sprintf('Subject %s  \n%.d Percent Coverage', ...
%        %     id{i}, percCoverage(i)));
%        % saveas(gcf, sprintf('PIP_spm12/Model_Second_Level/LASSO-pcr/ER_LNegLNeu_IMT/coverage/%s_coverage_%d_percent.png', ...
%        %     id{i}, percCoverage(i)));
%        % delete(hText);
%        % o2 = removeblobs(o2);
%         disp(i) 
%     end
% end
% 
%     
% 
% % find the subjects with the lowest coverage and plot their coverage
% [a, b] = sort(percCoverage)
