%% Construct dataset for IMT prediction paper

%  variables needed:
%       ID and IMT
%       images for 3 tasks: IAPS OF NF

% Produces output for Figure 1 (demographics)

% requries CANLab 'filenames.m' function to list files

%% cd to drive
cd /Volumes/ProjectDrive/AHAB_II/ML_projects/IMT_LassoPCR/datasets

%% load AHAB behavioral data
ahab = readtable('AHAB2_ID_IMT.csv');
%ahab_all = readtable('AHAB2_mega_490_11_27_2018_numeric_only.csv');
%ahab = ahab_all(:, {'LABID', 'AGE', 'SEX', 'mavg', 'BMI_V1', 'WAIST', 'HDL', 'TRIG', 'GLU', 'INSLN', 'SBPMS', 'DBPMS'})
%ahab.Properties.VariableNames = {'ID', 'age', 'sex', 'imt', 'BMI', 'waist_circumference', 'HDL', 'triglycerides', 'glucose', 'insulin', 'SBP', 'DBP'};
ahab.Properties.VariableNames = {'ID', 'IMT'};

%% load PIP behavioral data
pip = readtable('PIP_ID_IMT.csv');
%pip_all = readtable('PIP_n330_03_26_2019.csv', 'TreatAsEmpty', 'NA');
%pip_all.Properties.VariableNames{1} = 'id';
%pip = pip_all(:, {'id', 'age', 'gender', 'mavgimt', 'BMI', 'waist', 'hdl', 'trig', 'glucose', 'insulin', 'sbp_sit', 'dbp_sit'})
for i = 1:length(pip.id) % fix IDs (remove initials)
    pip.id{i} = pip.id{i}(1:4);
end
pip.id = str2double(pip.id);
%pip.Properties.VariableNames = {'ID', 'age', 'sex', 'imt', 'BMI', 'waist_circumference', 'HDL', 'triglycerides', 'glucose', 'insulin', 'SBP', 'DBP'};
pip.Properties.VariableNames = {'ID', 'IMT'};

%% join PIP & AHAB
dat = [ahab; pip];

%% make study labels
studies = {'PIP' 'AHAB'};
dat.study = dat.ID>10000;
dat.study = studies(dat.study+1)';

%% IAPS
filepath_IAPS = filenames('/Volumes/ProjectDrive/*/SPM12/First_Level/*/ER/con_0001.nii');
clear ID;
for i = 1:length(filepath_IAPS)
    tmp = strsplit(filepath_IAPS{i}, '/');
    ID(i) = str2double(tmp{7});
end
ID = ID';
IAPS = table(ID, filepath_IAPS);
dat = outerjoin(dat, IAPS);
dat.Properties.VariableNames{1} = 'ID';

%% Old Faces
filepath_OF = filenames('/Volumes/ProjectDrive/AHAB_II/SPM12/First_Level/*/OldFaces/con_0001.nii');
clear ID;
for i = 1:length(filepath_OF)
    tmp = strsplit(filepath_OF{i}, '/');
    ID(i) = str2double(tmp{7});
end
ID = ID';
OF = table(ID, filepath_OF);
dat = outerjoin(dat, OF);
dat.Properties.VariableNames{1} = 'ID';

%% New Faces
filepath_NF = filenames('/Volumes/ProjectDrive/AHAB_II/SPM12/First_Level/*/Faces/con_0001.nii');
clear ID;
for i = 1:length(filepath_NF)
    tmp = strsplit(filepath_NF{i}, '/');
    ID(i) = str2double(tmp{7});
end
ID = ID';
NF = table(ID, filepath_NF);
dat = outerjoin(dat, NF);
dat.Properties.VariableNames{1} = 'ID';

%% interim summary
disp('==== INTERIM SUMMARY (image counts) ====')
fprintf('\nSubjects across 2 studies: %d \n\n', size(dat, 1))
fprintf('AHAB subjects with OldFaces images: %d \n\n', length(filepath_OF))
fprintf('AHAB subjects with NewFaces images: %d \n\n', length(filepath_NF))
fprintf('Subjects with IAPS images: %d \n', length(filepath_IAPS))
fprintf('    AHAB: %d \n', sum(strcmp(dat.study, 'AHAB') & ~cellfun(@isempty, dat.filepath_IAPS)))
fprintf('    PIP: %d \n', sum(strcmp(dat.study, 'PIP') & ~cellfun(@isempty, dat.filepath_IAPS)))

%% remove bad subjects

% remove subject with brain abnormality
dat(isnan(dat.ID), :) % they have OF data only
dat(isnan(dat.ID), :) = []; %                                               -1 OF

% 48803 has poor coverage on IAPS and NF, but not OF  
dat{dat.ID == 48803,'filepath_IAPS'} = {''}; %                              -1 IAPS
dat{dat.ID == 48803,'filepath_NF'} = {''}; %                                -1 NF

% remove people who don't have IMT
dat(isnan(dat.IMT), :) % 2 subjects, have OF & NF but not IAPS
dat(isnan(dat.IMT), :) = []; %                                              -2 OF, -2 NF

%dat(dat.ID == 48803, {'BMI', 'waist_circumference', 'HDL', 'triglycerides', 'glucose', 'insulin', 'SBP', 'DBP'}) = []; % poor IAPS coverage (OF and NF is fine)

%% remove PIP subjects who also had AHAB data
dups = readtable('AHAB_PIP_EMA_BP_n44_04_25_2016.csv');
dups.Properties.VariableNames{1} = 'ID';
for i = 1:length(dups.ID)
    dups.ID{i} = dups.ID{i}(1:4);
end
dups.ID = str2double(dups.ID);
dat(ismember(dat.ID, dups.ID), :) % 15 PIP participants with IAPS
dat(ismember(dat.ID, dups.ID), :) = []; %                                   -15 IAPS

%% remove subjects who don't have any 1 of the 3 tasks
dat(cellfun(@isempty, dat.filepath_IAPS) & cellfun(@isempty, dat.filepath_NF) & cellfun(@isempty, dat.filepath_OF), :) = [];

%% final summary
disp('==== FINAL SUMMARY (after removing missing / problematic data) ====')
fprintf('\nSubjects across 2 studies: %d \n\n', size(dat, 1))
fprintf('AHAB subjects with OldFaces images: %d \n\n', sum(~cellfun(@isempty, dat.filepath_OF)))
fprintf('AHAB subjects with NewFaces images: %d \n\n', sum(~cellfun(@isempty, dat.filepath_NF)))
fprintf('Subjects with IAPS images: %d \n', sum(~cellfun(@isempty, dat.filepath_IAPS)))
fprintf('    AHAB: %d \n', sum(strcmp(dat.study, 'AHAB') & ~cellfun(@isempty, dat.filepath_IAPS)))
fprintf('    PIP: %d \n', sum(strcmp(dat.study, 'PIP') & ~cellfun(@isempty, dat.filepath_IAPS)))

%% write dataset to file
writetable(dat, 'ID_IMT_filepaths.csv');

