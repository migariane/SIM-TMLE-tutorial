## Targeted Maximum Likelihood Estimation for a Binary treatment: A tutorial. Statistics in Medicine, 2017
### Repository for the code presented in the manuscript (reproducible research in epidemiology and biostatisitcs)

### Miguel Angel Luque Fernandez, Michael Schomaker, Bernard Rachet, Mireille Schnitzer
### ABSTRACT
When estimating the average treatment effect for a binary outcome, methods that incorporate propensity scores, the g-formula, or targeted maximum likelihood estimation (TMLE) are preferred over naïve regression approaches which are typically biased. TMLE is a semiparametric and double-robust substitution estimator, suitable for the estimation of causal effects under standard causal assumptions. However, TMLE requires weaker assumptions than its competitors. We provide a step-by-step guided implementation of TMLE and illustrate it in a realistic scenario based on cancer epidemiology where assumptions about correct model specification and positivity (i.e., when a study participant had zero probability of receiving the treatment) are nearly violated. This article provides a concise and reproducible educational introduction to TMLE for a binary outcome and exposure. The reader should gain sufficient understanding of TMLE from this introductory tutorial to be able to apply the method in practice. Extensive R-code is provided in easy-to-read boxes throughout the article for replicability. Stata users will find a testing implementation of TMLE and additional material in the appendix and at the following GitHub repository: https://github.com/migariane/SIM-TMLE-tutorial


KEYWORDS: causal inference, machine learning, observational studies, targeted maximum likelihood estimation, Super Learner
