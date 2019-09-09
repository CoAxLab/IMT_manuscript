%% loop through a list of IDs to create tables of raw ratings and averaged ratings

% needs list of files

cd ~/Box/@Predict_IMT/ER_LNeg_LNeu_IMT/eprime

%% ahab

allT = table; allstats = table;
f = filenames('ahab/*txt')
for i = 1:length(f)
    disp(i)
    [all, stats] = readeprime_emoreg(f{i});
    allT = [allT; all];
    allstats = [allstats; struct2table(stats)];
end

%% add subject 40496 (manually extracted from txt file)
newsub.ID = 40496;
newsub.LookNeg_rating = 1.857;
newsub.LookNeut_rating = 1;
newsub.RegNeg_rating = 2.466;
allstats = [allstats; struct2table(newsub)];

%% write out to stats files

writetable(allT, 'ahab_emoreg_ratings_raw.csv')
writetable(allstats, 'ahab_emoreg_ratings_averages.csv')

%% pip
allT = table; allstats = table;
f = filenames('pip/*')
for i = 1:length(f)
    disp(i)
    [all, stats] = readeprime_emoreg(f{i});
    allT = [allT; all];
    allstats = [allstats; struct2table(stats)];
end

%% write out to stats files

writetable(allT, 'pip_emoreg_ratings_raw.csv')
writetable(allstats, 'pip_emoreg_ratings_averages.csv')

%% load and combine stats files

ahab = readtable('ahab_emoreg_ratings_averages.csv')
ahab_ids = dlmread('ER_correctLABIDS.txt') % some ahab IDs are wrong
for i = 1:length(ahab_ids)
    ahab.ID(ahab.ID == ahab_ids(i, 1)) = ahab_ids(i, 2);
end

pip = readtable('pip_emoreg_ratings_averages.csv')

all = [pip; ahab];
writetable(all, 'all_emoreg_ratings_averages.csv')

%% combine with IMT data
load ~/Box/@Predict_IMT/ER_LNeg_LNeu_IMT/data_merged_n338_ER_IMT_IL6_CRF.mat
all.Properties.VariableNames{1} = 'id'
dat = outerjoin(mergedn338(:, {'id', 'mavg'}), all)
dat.Properties.VariableNames{1} = 'id_mergedn338'
dat = sortrows(dat, 3, 'descend')
writetable(dat, 'all_mavg_ratings_crosscheck.csv')

%% correlate mavg with LookNeg-LookNeut rating
dat.LookNeg_LookNeut_rating = dat.LookNeg_rating - dat.LookNeut_rating

[r, pval] = corr(dat.LookNeg_LookNeut_rating, dat.mavg, 'rows', 'pairwise')
% r = -.028 p = .606

%%
