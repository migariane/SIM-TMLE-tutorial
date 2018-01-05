######################################################################################
# Miguel Angel Luque Fernandez, Michael Schomaker, Bernard Rachet, Mireille Schnitzer
# Targeted Maximum Likelihood Estimation for a Binary Treatment: A tutorial
# R-syntax included in the boxes of the manuscript
######################################################################################

# Function to generate data (DGP)
generateData<- function(n){
  w1 <- rbinom(n, size=1, prob=0.5)
  w2 <- rbinom(n, size=1, prob=0.65)
  w3 <- round(runif(n, min=0, max=4), digits=0)
  w4 <- round(runif(n, min=0, max=5), digits=0)
  A <- rbinom(n, size=1, prob= plogis(-5 +  0.05*w2 + 0.25*w3 + 0.6*w4 + 0.4*w2*w4))
  # counterfactual
  Y.1 <- rbinom(n, size=1, prob= plogis(-1 + 1 -0.1*w1 + 0.35*w2 + 0.25*w3 + 0.20*w4 + 0.15*w2*w4))
  Y.0 <- rbinom(n, size=1, prob= plogis(-1 + 0 -0.1*w1 + 0.35*w2 + 0.25*w3 + 0.20*w4 + 0.15*w2*w4))
  # Observed outcome
  Y <- Y.1*A + Y.0*(1 - A)
  # return data.frame
  data.frame(w1, w2, w3, w4, A, Y, Y.1, Y.0)
}

# True ATE in the population
set.seed(7777)
ObsDataTrueATE <- generateData(n = 5000000)
True_EY.1 <- mean(ObsDataTrueATE$Y.1)
True_EY.0 <- mean(ObsDataTrueATE$Y.0)
True_ATE  <- True_EY.1-True_EY.0 ;True_ATE
True_MOR  <- (True_EY.1*(1-True_EY.0))/((1-True_EY.1)*True_EY.0);True_MOR

cat("\n True_ATE:", abs(True_ATE))

# Data for analysis
set.seed(7722)
ObsData <- generateData(n = 10000)
write.csv(ObsData, "ObsData.csv")

# Naive approach: conditional odds ratio
naive <-glm(data = ObsData, Y ~ A + w1 + w2 + w3 + w4, family = binomial)
summary(naive)
exp(naive$coef[2])
exp(confint(naive))

## TMLE implementation by hand
# Step 1 estimation and prediction of the model for the outcome (G-computation)
gm  <- glm(Y ~ A + w1 + w2 + w3 + w4, family = binomial, data = ObsData)
# Prediction for A, A = 1 and, A = 0
QAW_0 <- predict(gm, type = "response")
Q1W_0 = predict(gm, newdata=data.frame(A = 1, ObsData [,c("w1","w2","w3","w4")]), type = "response")
Q0W_0 = predict(gm, newdata=data.frame(A = 0, ObsData [,c("w1","w2","w3","w4")]), type = "response")
# Estimated mortality risk difference
mean(Q1W_0 - Q0W_0)
# Estimated MOR
mean(Q1W_0)*(1-mean(Q0W_0))/((1-mean(Q1W_0))*mean(Q0W_0))

# Step 2 estimation and prediction of the propensity score (ps)
psm <- glm(A ~ w1 + w2 + w3 + w4, family = binomial, data = ObsData)
gW = predict(psm, type = "response")
summary(gW)

# Step 3 computation of H (clever covariates) and estimation of epsilon
H1W = ObsData$A / gW
H0W = (1-ObsData$A)  / (1 - gW)
epsilon <- coef(glm(ObsData$Y ~ -1 + H0W + H1W + offset(qlogis(QAW_0)), family = binomial))

# Step 4 Targeted estimate of the ATE and Marginal Odds Ratio
Q0W_1 <-  plogis(qlogis(Q0W_0) +  epsilon[1] / (1 - gW))
Q1W_1 <-  plogis(qlogis(Q1W_0) +  epsilon[2] / gW)

# ATE
ATEtmle1 <-  mean(Q1W_1 - Q0W_1); ATEtmle1
cat("\n ATEtmle1_bias:", abs(True_ATE - ATEtmle1))
cat("\n ATEtmle1_rel_bias:",abs(True_ATE - ATEtmle1)/True_ATE,"%")

# Marginal OR
tmle1.MOR <- mean(Q1W_1) * (1 - mean(Q0W_1)) / ((1 - mean(Q1W_1)) * mean(Q0W_1)); tmle1.MOR

# Table to visualize the data
psi <- Q1W_1 - Q0W_1 
library(DT)
df <- round(cbind(Q1W_0, Q0W_0, gW, eps1=epsilon[1], eps2=epsilon[2], psi), digits = 4)
datatable(head(df, n = nrow(df)), options = list(pageLength = 5, digits = 3))

# Step 5 statistical inference (efficient influence curve)
# Efficient influence curve ATE
EY1tmle<-mean(Q1W_1)
EY0tmle<-mean(Q0W_1)

d1 <- ((ObsData$A) * (ObsData$Y - Q1W_1)/gW) + Q1W_1 - EY1tmle
d0 <- ((1 - ObsData$A) * (ObsData$Y - Q0W_1))/(1 - gW)  + Q0W_1 - EY0tmle

IC <- d1 - d0
n <- nrow(ObsData)
varHat.IC <- var(IC) / n
ATEtmle1CI <- c(ATEtmle1 - 1.96 * sqrt(varHat.IC), ATEtmle1 + 1.96 * sqrt(varHat.IC)); ATEtmle1; ATEtmle1CI

# Efficient influence curve MOR
ICmor_tmle <- (1 - EY0tmle) / EY0tmle / (1 - EY1tmle)^2  * d1 - EY1tmle / (1 - EY1tmle) / EY0tmle^2  * d0
varHat2.IC <- var(ICmor_tmle) / n
tmle1_ORCI <- tmle1.MOR + c(-1.96,1.96)*sqrt(varHat2.IC); tmle1.MOR; tmle1_ORCI

# Augmented inverse probability treatment weighting (AIPTW) estimator
EY1aiptw <- mean((ObsData$A) * (ObsData$Y - Q1W_0) / gW + Q1W_0)
EY0aiptw <- mean((1 - ObsData$A) * (ObsData$Y - Q0W_0) / (1 - gW) + Q0W_0)

AIPTW_ATE <- EY1aiptw - EY0aiptw; AIPTW_ATE
cat("\n AIPTW_bias:", abs(True_ATE - AIPTW_ATE))
cat("\n AIPTW_rel_bias:",abs(True_ATE - AIPTW_ATE) / True_ATE,"%")

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

# R-package tmle (base implementation includes SL.step, SL.glm and SL.glm.interaction)
library(tmle)
library(SuperLearner)
TMLE2 <- tmle(Y = ObsData$Y, A = ObsData$A, W = ObsData[,c("w1", "w2", "w3", "w4")], family = "binomial")

#Note that the tmle function default bounds the probablities in the clever covariate denominators at 0.025.
#You can remove this bound by specifying gbound=0

ATEtmle2 <- TMLE2$estimates$ATE$psi;ATEtmle2
TMLE2$estimates$ATE$CI
MORtmle2 <- TMLE2$estimates$OR$psi;MORtmle2
TMLE2$estimates$OR$CI

cat("\n ATEtmle2_bias:", abs(True_ATE - ATEtmle2))
cat("\n ATEtmle2_Rel_bias:",abs(True_ATE - ATEtmle2) / True_ATE,"%")

# R-package tmle with user-selected Super learner libraries
library(tmle)
library(SuperLearner)

SL.library <- c("SL.glm","SL.step","SL.step.interaction", "SL.glm.interaction","SL.gam",
                "SL.randomForest", "SL.rpart") 

TMLE3 <- tmle(Y = ObsData$Y,A = ObsData$A,W = ObsData [,c("w1", "w2", "w3", "w4")], 
              family = "binomial", Q.SL.library = SL.library,g.SL.library = SL.library)#, gbound=0)

ATEtmle3 <- TMLE3$estimates$ATE$psi;ATEtmle3
TMLE3$estimates$ATE$CI
MORtmle3 <- TMLE3$estimates$OR$psi;MORtmle3
TMLE3$estimates$OR$CI

cat("\n ATEtmle3_bias:", abs(True_ATE - ATEtmle3))
cat("\n ATEtmle3_rel_bias:", abs(True_ATE - ATEtmle3) / True_ATE,"%")
