function [allTrials, stats] = readeprime_emoreg(fname)
% requires Greg Siegle's readeprime.m function
% also requires the output txt is named as below:
% fname = sprintf('ER pittsburgh task-%d-%d.txt', id, id);

%% get ID - asssumes similar file naming conventions
tmp = strsplit(fname, '-');
id = str2double(tmp{2});

%% load data
lookTrials = readeprime(fname, 'WatchProc', ...
    {'Order', 'PicValence', 'Procedure', 'Picture1', 'Rating.RESP', 'Rating.RT'}, ...
    0, -5, 30);
lookTrials = rmfield(lookTrials, 'fname');
lookTrials = struct2table(lookTrials);

regTrials = readeprime(fname, 'RegulateProc', ...
    {'Order', 'PicValence', 'Procedure', 'Picture1', 'RegRating.RESP', 'RegRating.RT'}, ...
    0, -5, 30);
regTrials = rmfield(regTrials, 'fname');
regTrials = struct2table(regTrials);
regTrials.Properties.VariableNames{5} = 'Rating_RESP';
regTrials.Properties.VariableNames{6} = 'Rating_RT';

%% make 'ID' vector
ID = repmat(id, 45, 1);


%% combine, add ID, and sort by order
allTrials = [regTrials; lookTrials];
allTrials = [table(ID) allTrials];
allTrials.Properties.VariableNames{2} = 'Trial';
allTrials = sortrows(allTrials, 'Trial');
allTrials.Trial = [1:45]';

%% make some edits
allTrials.Properties.VariableNames{5} = 'IAPS';
allTrials.Procedure = nominal(allTrials.Procedure);
allTrials.PicValence = nominal(allTrials.PicValence);

%% remove -999 values
allTrials.Rating_RESP(allTrials.Rating_RESP < 0) = NaN;
allTrials.Rating_RT(isnan(allTrials.Rating_RESP)) = NaN;

%% calculate means
stats.ID = id;
stats.LookNeg_rating = nanmean(table2array(allTrials(allTrials.PicValence=='negative' & allTrials.Procedure=='WatchProc', 'Rating_RESP')));
stats.LookNeut_rating = nanmean(table2array(allTrials(allTrials.PicValence=='neutral' & allTrials.Procedure=='WatchProc', 'Rating_RESP')));
stats.RegNeg_rating = nanmean(table2array(allTrials(allTrials.Procedure=='RegulateProc', 'Rating_RESP')));

