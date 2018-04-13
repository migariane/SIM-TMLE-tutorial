//#######################################################################################
// Miguel Angel Luque Fernandez, Micheal Schomaker, Bernard Rachet, Mireille Schnitzer
// Targeted Maximum Likelihood Estimation for a Binary Treatment: A tutorial
// Stata simple TMLE implementation 
//#######################################################################################

cd "your path to the data"
import delimited ObsData.csv, clear
set more off

* Step 1: prediction model for the outcome Q0 (g-computation)
glm y a w1 w2 w3 w4, fam(binomial)
predict double QAW_0, mu
gen aa=a
replace a = 0
predict double Q0W_0, mu
replace a = 1
predict double Q1W_0, mu
replace a = aa
drop aa

// Q to logit scale
gen logQAW = log(QAW / (1 - QAW))
gen logQ1W = log(Q1W / (1 - Q1W))
gen logQ0W = log(Q0W / (1 - Q0W))

* Step 2: prediction model for the treatment g0 (IPTW)
glm a w1 w2 w3 w4, fam(binomial)
predict gw, mu
gen double H1W = a / gw
gen double H0W = (1 - a) / (1 - gw)

* Step 3: Computing the clever covariate H(A,W) and estimating the parameter (epsilon) (MLE)
glm y H1W H0W, fam(binomial) offset(logQAW) noconstant
mat a = e(b)
gen eps1 = a[1,1]
gen eps2 = a[1,2]

* Step 4: update from Q0 to Q1
gen double Q1W_1 = exp(eps1 / gw + logQ1W) / (1 + exp(eps1 / gw + logQ1W))
gen double Q0W_1 = exp(eps2 / (1 - gw) + logQ0W) / (1 + exp(eps2 / (1 - gw) + logQ0W))

* Step 5: Targeted estimate of the ATE 
gen ATE = (Q1W_1 - Q0W_1)
summ ATE
global ATE = r(mean)

* Step 6: Statistical inference (efficient influence curve)
qui sum(Q1W_1)
gen EY1tmle = r(mean)
qui sum(Q0W_1)
gen EY0tmle = r(mean)

gen d1 = ((a * (y - Q1W_1)/gw)) + Q1W_1 - EY1tmle
gen d0 = ((1 - a) * (y - Q0W_1)/(1 - gw))  + Q0W_1 - EY0tmle

gen IC = d1 - d0
qui sum IC
gen varIC = r(Var) / r(N)

global LCI =  $ATE - 1.96*sqrt(varIC)
global UCI =  $ATE + 1.96*sqrt(varIC)
display "ATE:"  %05.4f  $ATE _col(15) "95%CI: " %05.4f  $LCI "," %05.4f  $UCI
