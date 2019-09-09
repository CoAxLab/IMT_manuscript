function out = IMT_reliability(task, roi)
%% compute Spearman Brown split-half internal consistency for the IMT paper
% roi can be 3 options
%   voxelwise - computes voxelwise map of split-half internal consistency
%           restricts analysis to voxels for which half of participants have coverage
%   pattern - computes split-half internal consistency for the multivariate pattern generated from LASSO-PCR
%   an actual ROI - computes split-half for the average response in an ROI (e.g., amygdala, dACC, anterior insula)

wd = cdtodrive; cd AHAB_II/ML_projects/IMT_LassoPCR

maskgm = char(filenames('masks/resliced_grey25grey25.nii', 'absolute'));

switch roi
    case 'Amygdala'
        mask = filenames('masks/Amygdala.nii', 'absolute');
    case 'dACC'
        mask = filenames('masks/dACC.nii', 'absolute');
    case 'insula'
        mask = filenames('masks/anteriorInsula.nii', 'absolute');
end

switch task
    case 'OF'
        cd Faces/OF/AllFaces/
    case 'NF'
        cd Faces/NF/AllFaces/
    case 'IAPS'
        cd IAPS/
end

%% load task-specific data set
load('data.mat');
if strcmp(task, 'OF') % remove 33071 when looking at OF
    out.cv_folds.fold_indicator(dat.ID == 33071) = [];
    dat(dat.ID == 33071, :) = [];
end

%% access image files
imagefiles = table2cell(dat(:, contains(dat.Properties.VariableNames, 'filepath')));
tmp = strsplit(imagefiles{1}, '/');
taskname = tmp{end-1};

%% load images from both split-half models
blockA = fmri_data(strrep(imagefiles, taskname, ['SplitHalf/' taskname '_1stHalf']), maskgm);
blockB = fmri_data(strrep(imagefiles, taskname, ['SplitHalf/' taskname '_2ndHalf']), maskgm);

%%
switch roi
    case 'voxelwise' %% if looking at voxelwise reliability
        clear dat
        dat(:,:,1) = blockA.dat; dat(:,:,2) = blockB.dat;

        %% convert zeros to NaN so they are ignored in the correlations
        dat(dat == 0 ) = NaN;

        %% mask according to how many people have all data (not NaN) in a voxel
        for i = 1:size(dat, 1)
            %voxelmask(i) = sum(all(squeeze(dat(i, : , :)),2));
            voxelmask(i) = sum(any(~isnan(squeeze(dat(i,:,:))), 2));
        end

        %% number of people who must have a value in a voxel for that voxel to be included
        % go with 50% for now
        thresh = round(size(dat,2)/2); 

        %% calculate correlation for each voxel
        parfor i = 1:size(dat, 1)
            corrAB(i) = corr(dat(i, :, 1)', dat(i, :, 2)', 'rows', 'pairwise');
        end

        %% calculate SBreliability
        SBreliability = (2*corrAB)./(1+corrAB);

        %% apply mask to ignore bad SB values
        % thresh based on number/proportion of people who must have data in a voxel
        % go with 50% for now
        thresh = round(size(dat,2)/2); 
        SBreliability(voxelmask< thresh) = 0;
        
        %% generate map
        out = blockA; out.dat = SBreliability';

        %% view and save
       orthviews(out)
       write(out, 'fname', 'Spearman_Brown_reliability_map.nii')
        
    case 'pattern' % get the reliability of the weightmap generated from LASSO-PCR
%         for i = 1:out.cv_folds.NumTestSets
%             patternA(i).dat = blockA.dat(:, out.cv_folds.fold_indicator == i)'*out.inner_training_folds(i).weight_obj.dat;
%             patternB(i).dat = blockB.dat(:, out.cv_folds.fold_indicator == i)'*out.inner_training_folds(i).weight_obj.dat;
%         end
%         patternA = vertcat(patternA.dat);
%         patternB = vertcat(patternB.dat);
            
        %quicker way, but might be overfitting
        patternA = blockA.dat'*out.model_all_data.weight_obj.dat;
        patternB = blockB.dat'*out.model_all_data.weight_obj.dat;
        
        corrAB = corr(patternA, patternB)

        out = (2*corrAB)./(1+corrAB);
        
        figure; scatter(patternA, patternB); lsline;
        xlabel('First Half'); ylabel('Second Half');
        title(sprintf('Internal Consistency for %s in the %s', task, roi));

    otherwise % individual ROI
        %% extract ROI data from each half
        roiA = extract_roi_averages(blockA, char(mask));
        out.roiA = roiA.dat;
        roiB = extract_roi_averages(blockB, char(mask));
        out.roiB = roiB.dat;

        %% calculate correlation and Spearman-Brown corrected reliability
        corrAB = corr(out.roiA, out.roiB)

        out.SB = (2*corrAB)./(1+corrAB);
        
        figure; scatter(out.roiA, out.roiB); lsline;
        xlabel('First Half'); ylabel('Second Half')
        title(sprintf('Internal Consistency for %s in the %s', task, roi));
end
    
    
end
