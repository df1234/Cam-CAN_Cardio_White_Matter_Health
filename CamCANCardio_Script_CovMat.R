################################################################################
############## MODELLING CARDIO-VASCULAR AND WHITE MATTER HEALTH ###############
############################# IN THE CAM-CAN COHORT ############################           
################################################################################

# This script produces key analyses for the paper

# Fuhrmann, D., Nesbitt, D., Shafto, M., Rowe, J., Price, D., Gadie, A., 
#    Cam-CAN & Kievit, R. A.(2018). Strong and specific associations between 
#    cardiovascular risk factors and brain white matter microstructure- and 
#    macrostructure in healthy aging. Neuobiology of Aging, 
#    doi: 10.1016/j.neurobiolaging.2018.10.005

# The script uses the covariance matrix. Estimates and test statistics will differ 
# somewhat from values reported in the paper because of missingness. To reproduce 
# results exactly, please request raw data from 
# https://camcan-archive.mrc-cbu.cam.ac.uk/dataaccess/

# Script author: Delia Fuhrmann, January 2017, using RStudio 1.0.143 and R 
# version 3.4.0 (You Stupid Darkness)

################################################################################

### Load lavaan
# install.packages("lavaan") # this only needs to be done once
library(lavaan)

### Load covariance matrix
cov.matrix = read.csv("CamCANCardio_CovarianceMatrix.csv", row.names = 1)

################################################################################

### Measuremment model of cardiovascular health
model_ThreeFactor<-
  '
# LVs
DPB =~ bp_dia1 + bp_dia2 + bp_dia3
SPB =~ bp_sys1 + bp_sys2 + bp_sys3
HR  =~ pulse1  + pulse2  + pulse3
'

fit_ThreeFactor <- cfa(model_ThreeFactor, sample.cov=as.matrix(cov.matrix), 
                       sample.nobs=667, std.lv=T)

summary(fit_ThreeFactor, standardized=T, rsquare=T, fit.measures=T)

################################################################################

### Model of of cardiovascular health and lesion burden
model_lesion <-
  '
# LVs
DBP =~ bp_dia1 + bp_dia2 + bp_dia3
SBP =~ bp_sys1 + bp_sys2 + bp_sys3
HR  =~ pulse1  + pulse2  + pulse3

# regressions
lst_07_tlv ~ DBP + SBP + HR + age
lst_07_n   ~ DBP + SBP + HR + age

# covariances
age ~~ DBP + SBP + HR
'

fit_lesion <- sem(model_lesion, sample.cov=as.matrix(cov.matrix), 
                  sample.nobs=272, std.lv=T)

summary(fit_lesion, standardized=T, rsquare=T, fit.measures=T)

################################################################################

### Model of of cardiovascular health and white matter microstructure (here MD)
model_MD <-
  '
# LVs
DBP =~ bp_dia1 + bp_dia2 + bp_dia3
SBP =~ bp_sys1 + bp_sys2 + bp_sys3
HR  =~ pulse1  + pulse2  + pulse3

# Regressions
md_UF   ~ DBP + SBP + HR + age
md_SLF  ~ DBP + SBP + HR + age
md_IFOF ~ DBP + SBP + HR + age
md_ATR  ~ DBP + SBP + HR + age
md_CST  ~ DBP + SBP + HR + age
md_FMaj ~ DBP + SBP + HR + age
md_FMin ~ DBP + SBP + HR + age
md_CG   ~ DBP + SBP + HR + age
md_CH   ~ DBP + SBP + HR + age
md_ILF  ~ DBP + SBP + HR + age

# Covariances
age ~~ DBP + SBP + HR
'

fit_MD <- sem(model_MD, sample.cov=as.matrix(cov.matrix), sample.nobs=667, std.lv=T)

summary(fit_MD, standardized=T, rsquare=T, fit.measures=T)

