# generate Table 1 (demographics) given IDs used in the paper
# output csv file containing all demographic info

library(haven)
library(openxlsx)
library(psych)
library(sjPlot)

# load file listing participant IDs used in paper

dat <- read.csv("~/Box/@Predict_IMT/Datasets/ID_IMT_filepaths.csv")
dim(dat) # should be 617

# how many AHAB participants? 
sum(dat$study == 'AHAB') # 458

# how many PIP participants?
sum(dat$study == 'PIP')  # 159

# load AHAB and PIP master files

ahab <- read_sav("~/Box/@Predict_IMT/Datasets/AHAB2_mega_490_11_27_2018.sav")
ahab$ID = ahab$LABID

pip <- read_sav("~/Box/@Predict_IMT/Datasets/PIP_n330_03_26_2019.sav")
pip$ID = as.double(substr(pip$id, 1, 4))

# get MRI and URL dates, calculate days between
pip$days_between = as.numeric(pip$date.mri - pip$date.url)
ahabdays = read.xlsx('~/Box/@Predict_IMT/Datasets/490 dates 081919.xlsx')
ahabdays$URL = as.Date(ahabdays$URL, '%m/%d/%y')
ahabdays$MRI = as.Date(ahabdays$MRI, '%m/%d/%y')
ahabdays$days_between = as.numeric(ahabdays$MRI - ahabdays$URL)
ahab = merge(ahab, ahabdays[,c('LABID', 'days_between')])

# for each study:
#   get variables, 
#   merge with all_data, reduce,
#   double-check distributions, and 
#   generate descriptive statistics

ahab <- ahab[, c('ID', 'SEX', 'AGE', 'RACE', 'WAIST', 'BMI_V1', 'CHOL', 
                'TRIG', 'HDL', 'GLU', 'SBPMS', 'DBPMS','clinic_HR', 
                'mavg', 'YRS_SCH', 'OldFaces_Faces_ACC', 
                'Faces_Average_ACC', 'ER_LookNegative_rating', 
                'ER_LookNeutral_rating', 'SMK_STAT', 'days_between')]
ahab <- merge(dat[,c('ID', 'study')], ahab, by = 'ID')
dim(ahab) # should be 458
psych::describe(ahab)
mean(ahab$mavg)
sd(ahab$mavg)
# write table to file (ignore stats of categorical variables)
write.xlsx(data.frame(psych::describe(ahab)), file = '~/Box/@Predict_IMT/Tables/Table_1_Demographics_AHAB_TEK.xlsx', rowNames = T)


## need to circle back to categorical variables (sex, race, smoking)
sjt.frq(ahab$SEX)
sjt.frq(ahab$RACE)
sjt.frq(ahab$SMK_STAT)

pip <- pip[, c('ID', 'gender', 'age', 'race', 'waist', 'BMI', 'choles', 
                 'trig', 'hdl', 'glucose', 'sbp.sit', 'dbp.sit','hr.sit', 
                 'mavgimt', 'Yrs_School', 'LookNeg_Rating', 
                 'LookNeutral_Rating', 'smoke', 'days_between')]
pip <- merge(dat[,c('ID', 'study')], pip, by = 'ID')
dim(pip) # should be 159
psych::describe(pip)

## need to circle back to categorical variables (sex, race, smoking)
sjt.frq(pip$gender)
sjt.frq(pip$race)
sjt.frq(pip$smoke)

# write table to file (ignore stats of categorical variables)
write.xlsx(data.frame(psych::describe(pip)), file = '~/Box/@Predict_IMT/Tables/Table_1_Demographics_PIP_TEK.xlsx', rowNames = T)

## combine AHAB and PIP dataframes and save
ahab$Faces_Average_ACC <- NULL
ahab$OldFaces_Faces_ACC <- NULL
colnames(ahab) = c(c('ID', 'study', 'sex', 'age', 'race', 'waist', 'BMI', 'cholesterol', 
                     'triglycerides', 'HDL', 'glucose', 'SBP', 'DBP','HR', 
                     'IMT', 'years_school', 'LookNeg_Rating', 
                     'LookNeut_Rating', 'smoke', 'days_between'))
colnames(pip) = c(c('ID', 'study', 'sex', 'age', 'race', 'waist', 'BMI', 'cholesterol', 
                     'triglycerides', 'HDL', 'glucose', 'SBP', 'DBP','HR', 
                     'IMT', 'years_school', 'LookNeg_Rating', 
                     'LookNeut_Rating', 'smoke', 'days_between'))

alldat = rbind(pip, ahab)
write.csv(alldat, file = '~/Box/@Predict_IMT/Datasets/ID_IMT_demographics.csv')


