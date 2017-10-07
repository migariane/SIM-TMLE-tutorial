######################################################################################
# Miguel Angel Luque Fernandez, Micheal Schomaker, Bernard Rachet, Mireille Schnitzer
# Targeted Maximum Likelihood Estimation for a Binary Outcome: A tutorial
# Table-2: R-syntax for simulations
######################################################################################

# Super Learner libraries
SL.library <- c("SL.glm","SL.step","SL.step.interaction","SL.glm.interaction",
                "SL.gam","SL.randomForest","SL.glmnet")

# Data generation A: dual misspecification for the model of the outcome and treatment
set.seed(7777)
generateData<- function(n){
    w1 <- rbinom(n, size=1, prob=0.5)
    w2 <- rbinom(n, size=1, prob=0.65)
    w3 <- round(runif(n, min=0, max=4), digits=0)
    w4 <- round(runif(n, min=0, max=5), digits=0)
    A <- rbinom(n, size=1, prob= plogis(-5 +  0.5*w2 + 0.25*w3 + 0.6*w4 + 0.4*w2*w4)) 
    
    # counterfactuals
    Y.1 <- rbinom(n, size=1, prob= plogis(-1 + 1 -0.1*w1 + 0.35*w2 + 0.25*w3 + 0.20*w4 + 0.15*w2*w4))
    Y.0 <- rbinom(n, size=1, prob= plogis(-1 + 0 -0.1*w1 + 0.35*w2 + 0.25*w3 + 0.20*w4 + 0.15*w2*w4))
    # Observed outcome
    Y <- Y.1*A + Y.0*(1 - A)
    # return data.frame
    data.frame(w1, w2, w3, w4, A, Y, Y.1, Y.0)
}

# Data generation B: misspecification for the model of the outcome
# set.seed(7777)
# generateData<- function(n){
# w1 <- rbinom(n, size=1, prob=0.5)
# w2 <- rbinom(n, size=1, prob=0.65)
# w3 <- round(runif(n, min=0, max=4), digits=0)
# w4 <- round(runif(n, min=0, max=5), digits=0)
#   A <- rbinom(n, size=1, prob= plogis(-5 + 0.7*w1 + 0.5*w2 + 0.25*w3 + 0.6*w4))
# counterfactuals
# Y.1 <- rbinom(n, size=1, prob= plogis(-1 + 1 -0.1*w1 + 0.35*w2 + 0.25*w3 + 0.15*w2*w4))
# Y.0 <- rbinom(n, size=1, prob= plogis(-1 + 0 -0.1*w1 + 0.35*w2 + 0.25*w3 + 0.15*w2*w4))
# Observed outcome
# Y <- Y.1*A + Y.0*(1 - A)
# return data.frame
# data.frame(w1, w2, w3, w4, A, Y, Y.1, Y.0)
#}

# True ATE
ObsDataTrueATE   <- generateData(n=5000000)
True_ATE  <- mean(ObsDataTrueATE$Y.1 - ObsDataTrueATE$Y.0);True_ATE
True_EY.1 <- mean(ObsDataTrueATE$Y.1)
True_EY.0 <- mean(ObsDataTrueATE$Y.0)
True_MOR  <- (True_EY.1*(1-True_EY.0)) / ((1-True_EY.1)*True_EY.0);True_MOR

#Simulations
library(tmle)
library(SuperLearner)
R <- 1000  
#Empty vectors 
naive_OR  <- rep(NA,R)
ATEtmle1  <- rep(NA,R)
MORtmle1  <- rep(NA,R)
AIPTW     <- rep(NA,R)
MOR_AIPTW <- rep(NA,R)
ATEtmle2  <- rep(NA,R)
MORtmle2  <- rep(NA,R)
ATEtmle3  <- rep(NA,R)
MORtmle3  <- rep(NA,R)

for(r in 1:R){  
    print(paste("This is simulation run number",r)) 
    CancerData <- generateData(n=1000) 
    
    # ATE naive approach
    naive_OR[r] <- exp(glm(data = CancerData, Y ~ A + w1 + w2 + w3 + w4, family = "binomial")$coef[2])
    
    # TMLE implementation by hand
    # Step 1
    gm <- glm(Y ~ A + w1 + w2 + w3 + w4, family="binomial", data=CancerData)
    
    # Prediction for A, A=1 and, A=0
    QAW <- predict(gm)
    Q1W = predict(gm, newdata=data.frame(A = 1, CancerData[,c("w1","w2","w3","w4")]))
    Q0W = predict(gm, newdata=data.frame(A = 0, CancerData[,c("w1","w2","w3","w4")]))
    
    # Step 2 estimation of the propensity score (ps)
    psm <- glm(A ~ w1 + w2 + w3 + w4, family = binomial, data=CancerData)
    gW = predict(psm, type = "response")
    g1W = (1/gW)
    g0W = (-1/(1-gW))
    
    # Step 3 computation of H(A,W) and estimation of epsilon
    HAW <- (CancerData$A / gW -(1-CancerData$A) / (1 - gW))
    H1W = (1/gW)
    H0W = (-1 / (1 - gW))
    epsilon <- coef(glm(CancerData$Y ~ -1 + HAW + offset(QAW), family = "binomial")) 
    
    # Step 4 Targeted estimate of the ATE 
    ATEtmle1[r] <- mean(plogis(Q1W + epsilon * H1W) - plogis(Q0W + epsilon * H0W))
    
    # Step 5 Targeted estimate of the MOR
    T1.EY1 <- mean(plogis(Q1W + epsilon * H1W))
    T1.EY0 <- mean(plogis(Q0W + epsilon * H0W))
    MORtmle1[r] <- (T1.EY1*(1-T1.EY0)) / ((1-T1.EY1)*T1.EY0)
    
    # Augmented inverse probability treatment weighting (AIPTW) estimator
    AIPTW[r]  <- mean((HAW*(CancerData$Y - plogis(QAW)) + (plogis(Q1W)-plogis(Q0W))))
    AIPTW1 <- mean(CancerData$A * (CancerData$Y - plogis(Q1W)) / gW + plogis(Q1W))
    AIPTW0 <- mean((1- CancerData$A) * (CancerData$Y - plogis(Q0W)) / (1 - gW) + plogis(Q0W))
    MOR_AIPTW[r] <- (AIPTW1 * (1- AIPTW0)) / ((1- AIPTW1) * AIPTW0)  
    
    # R-package tmle (base implementation includes SL.step, SL.glm and SL.glm.interaction)
    ATE2 <- tmle(Y=CancerData$Y, A=CancerData$A, W=CancerData[,c("w1","w2","w3","w4")], family="binomial")
    ATEtmle2[r] <- ATE2$estimates$ATE$psi
    MORtmle2[r] <- ATE2$estimates$OR$psi
    # User-selected Super learner libraries
    ATE3 <- tmle(Y = CancerData$Y, A=CancerData$A, W=CancerData[,c("w1","w2","w3","w4")], family="binomial", Q.SL.library=SL.library, g.SL.library=SL.library)
    ATEtmle3[r] <- ATE3$estimates$ATE$psi
    MORtmle3[r] <- ATE3$estimates$OR$psi
}
# True ATE
True_ATE
# True MOR
True_MOR
# Mean naive
mean(naive_OR)
# Mean AIPTW
mean(AIPTW)
mean(MOR_AIPTW)
# Estimate of TMLE
mean(ATEtmle1)
mean(MORtmle1)
# Estimate of TMLE + SL
mean(ATEtmle2)
mean(MORtmle2)
# Estimate of TMLE + SL2
mean(ATEtmle3)
mean(MORtmle3)

save.image("your path\results.RData")
