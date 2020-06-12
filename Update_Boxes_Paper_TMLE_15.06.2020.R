###################################################################################################################################################
# Supplementary Material to Luque-Fernandez MA, Schomaker M, Rachet B, Schnitzer ME.                                                              #
# Targeted maximum likelihood estimation for a binary treatment: A tutorial.                                                                      #
# Statistics in Medicine. 2018;37(16):2530-46.                                                                                                    #
#                                                                                                                                                 #
#  BOXES FROM PAPER                                                                                                                               #
###################################################################################################################################################

### BOX 1 ###

# Function to generate data (data-generating process)
generateData<- function(n){
  w1 <- rbinom(n, size=1, prob=0.5)
  w2 <- rbinom(n, size=1, prob=0.65)
  w3 <- round(runif(n, min=0, max=4), digits=0)
  w4 <- round(runif(n, min=0, max=5), digits=0)
  A <- rbinom(n, size=1, prob= plogis(-5 +  0.05*w2 + 0.25*w3 + 0.6*w4 + 0.4*w2*w4))
  # counterfactual outcome
  Y.1 <- rbinom(n, size=1, prob= plogis(-1 + 1 -0.1*w1 + 0.35*w2 + 0.25*w3 + 0.20*w4 + 0.15*w2*w4))
  Y.0 <- rbinom(n, size=1, prob= plogis(-1 + 0 -0.1*w1 + 0.35*w2 + 0.25*w3 + 0.20*w4 + 0.15*w2*w4))
  # Observed outcome
  Y <- Y.1*A + Y.0*(1 - A)
  # return data.frame
  data.frame(w1, w2, w3, w4, A, Y, Y.1, Y.0)
}

# True ATE and MOR
set.seed(7777)
ObsDataTrueATE <- generateData(n = 5000000)
True_EY.1 <- mean(ObsDataTrueATE$Y.1)
True_EY.0 <- mean(ObsDataTrueATE$Y.0)
True_ATE  <- True_EY.1-True_EY.0 ;True_ATE
True_MOR  <- (True_EY.1*(1-True_EY.0))/((1-True_EY.1)*True_EY.0);True_MOR

cat("\n True_ATE:", abs(True_ATE))
cat("\n True_ATE:", True_MOR)

# Observed data for example
set.seed(7722)
ObsData <- generateData(n = 10000)

## TMLE implementated manually

### BOX 2 ###
# Step 1: estimation and prediction of the model for the outcome (G-computation)
gm  <- glm(Y ~ A + w1 + w2 + w3 + w4, family = binomial, data = ObsData)
# Prediction for A, A = 1 and, A = 0
QAW_0 <- predict(gm, type = "response")
Q1W_0 = predict(gm, newdata=data.frame(A = 1, ObsData [,c("w1","w2","w3","w4")]), type = "response")
Q0W_0 = predict(gm, newdata=data.frame(A = 0, ObsData [,c("w1","w2","w3","w4")]), type = "response")
# Estimated mortality risk difference
mean(Q1W_0 - Q0W_0)
# Estimated MOR
mean(Q1W_0)*(1-mean(Q0W_0))/((1-mean(Q1W_0))*mean(Q0W_0))

### BOX 3 ###
# Step 2: estimation and prediction of the propensity score (ps)
psm <- glm(A ~ w1 + w2 + w3 + w4, family = binomial, data = ObsData)
gW = predict(psm, type = "response")
summary(gW)

### BOX 4 ###
# Step 3: computation of H (clever covariate) and estimation of epsilon
H1W = ObsData$A / gW
H0W = (1-ObsData$A)  / (1 - gW)
epsilon <- coef(glm(ObsData$Y ~ -1 + H0W + H1W + offset(qlogis(QAW_0)), family = binomial))
epsilon


### BOX 5 ###
# Step 4: TARGETED estimate of the ATE and MOR
Q0W_1 <-  plogis(qlogis(Q0W_0) +  epsilon[1] / (1 - gW))
Q1W_1 <-  plogis(qlogis(Q1W_0) +  epsilon[2] / gW)

EY1tmle<-mean(Q1W_1)
EY0tmle<-mean(Q0W_1)

ATEtmle1 <-  mean(Q1W_1 - Q0W_1); ATEtmle1
tmle1.MOR <- mean(Q1W_1) * (1 - mean(Q0W_1)) / ((1 - mean(Q1W_1)) * mean(Q0W_1)); tmle1.MOR


# TABLE 1 from paper: visualize data (note, datatable requires package DT)
psi <- Q1W_1 - Q0W_1 
library(DT)
df <- round(cbind(Q1W_0, Q0W_0, gW, eps1=epsilon[1], eps2=epsilon[2], psi), digits = 4)
datatable(head(df, n = nrow(df)), options = list(pageLength = 5, digits = 3))

### BOX 6 ###
# Step 5: statistical inference (efficient influence curve)
# Efficient influence curve: ATE
d1 <- ((ObsData$A) * (ObsData$Y - Q1W_1)/gW) + Q1W_1 - EY1tmle
d0 <- ((1 - ObsData$A) * (ObsData$Y - Q0W_1))/(1 - gW)  + Q0W_1 - EY0tmle

IC <- d1 - d0
n <- nrow(ObsData)
varHat.IC <- var(IC) / n
ATEtmle1CI <- c(ATEtmle1 - 1.96 * sqrt(varHat.IC), ATEtmle1 + 1.96 * sqrt(varHat.IC)); ATEtmle1; ATEtmle1CI

# Efficient influence curve: MOR
ICmor_tmle <- (1 - EY0tmle) / EY0tmle / (1 - EY1tmle)^2  * d1 - EY1tmle / (1 - EY1tmle) / EY0tmle^2  * d0
varHat2.IC <- var(ICmor_tmle) / n
tmle1_ORCI <- tmle1.MOR + c(-1.96,1.96)*sqrt(varHat2.IC); tmle1.MOR; tmle1_ORCI


### Box 7 ###
# Naive approach: CONDITIONAL odds ratio
naive <-glm(data = ObsData, Y ~ A + w1 + w2 + w3 + w4, family = binomial)
summary(naive)
exp(naive$coef[2])
exp(confint(naive))

### Box 8 ###
# Augmented inverse probability treatment weighting (AIPTW) estimator
EY1aiptw <- mean((ObsData$A) * (ObsData$Y - Q1W_0) / gW + Q1W_0)
EY0aiptw <- mean((1 - ObsData$A) * (ObsData$Y - Q0W_0) / (1 - gW) + Q0W_0)

AIPTW_ATE <- EY1aiptw - EY0aiptw; AIPTW_ATE
D1 <- (ObsData$A) * (ObsData$Y - Q1W_0) / gW + Q1W_0 - EY1aiptw 
D0 <- (1 - ObsData$A) * (ObsData$Y - Q0W_0) / (1 - gW) + Q0W_0 - EY0aiptw
varHat_AIPTW <- var(D1 - D0) / n

# AIPTW ATE 95%CI
ATEaiptw_CI <- c(AIPTW_ATE - 1.96 * sqrt(varHat_AIPTW), AIPTW_ATE + 1.96 * sqrt(varHat_AIPTW)); AIPTW_ATE; ATEaiptw_CI

# AIPTW MOR 95%CI
AIPTW_MOR <- (EY1aiptw * (1 - EY0aiptw))/((1 - EY1aiptw) * EY0aiptw);AIPTW_MOR
ICmor_aiptw <- (1 - EY0aiptw) / EY0aiptw / (1 - EY1aiptw)^2 * D1 - EY1aiptw / (1 - EY1aiptw) / EY0aiptw^2 * D0
varHat_AIPTW2 <- var(ICmor_aiptw) / n
MORaiptw_CI <-c(AIPTW_MOR - 1.96*sqrt(varHat_AIPTW2), AIPTW_MOR + 1.96*sqrt(varHat_AIPTW2)); AIPTW_MOR; MORaiptw_CI


### Box 9 ###
# R-package tmle 
# Note: in 2018, the base super learner library was SL.step, SL.glm and SL.glm.interaction
# this changed and thus we explicitely specify the library here; still the results are not exactly identical to those reported in the paper
library(tmle)
library(SuperLearner)
TMLE2 <- tmle(Y = ObsData$Y, A = ObsData$A, W = ObsData[,c("w1", "w2", "w3", "w4")], family = "binomial", 
              Q.SL.library = c("SL.step", "SL.glm", "SL.glm.interaction"), g.SL.library = c("SL.step", "SL.glm", "SL.glm.interaction"))

#Note that the tmle function default bounds the probablities in the clever covariate denominators at 0.025.
#You can change this bound by specifying gbound=0.01 (or similar)

ATEtmle2 <- TMLE2$estimates$ATE$psi;ATEtmle2
TMLE2$estimates$ATE$CI
MORtmle2 <- TMLE2$estimates$OR$psi;MORtmle2
TMLE2$estimates$OR$CI

### Box 10 ###
# R-package tmle with user-selected Super learner libraries
# note: you need the following packagesto be installed (not necessarily loaded): randomForest, gam, rpart
# as in Box 9, there are minor changes in the results due to updates of packages
SL.library <- c("SL.glm","SL.step","SL.step.interaction", "SL.glm.interaction","SL.gam",
                "SL.randomForest", "SL.rpart") 

TMLE3 <- tmle(Y = ObsData$Y,A = ObsData$A,W = ObsData [,c("w1", "w2", "w3", "w4")], 
              family = "binomial", Q.SL.library = SL.library,g.SL.library = SL.library)#, gbound=0)

ATEtmle3 <- TMLE3$estimates$ATE$psi;ATEtmle3
TMLE3$estimates$ATE$CI
MORtmle3 <- TMLE3$estimates$OR$psi;MORtmle3
TMLE3$estimates$OR$CI

