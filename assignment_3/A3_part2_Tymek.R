# File Description -------------------------------------------------------------
#
#   Tymoteusz Barcinski - s221937
#   
#   Advanced Dataanalysis and Statistical Modelling
#   Assignment 3 - Mixed effects and hierarchical models
#   Part 2 - Hierarchical Models: Random variance
#
#_______________________________________________________________________________
rm(list=ls()) 
print(utils::getSrcDirectory(function(){}))
print(utils::getSrcFilename(function(){}, full.names = TRUE))
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()
options(scipen=0)
Sys.setenv(LANG = "en")

library(ggplot2)
library(MASS)
library(dplyr)
# library(tsibble)
library(forecast)
# library(matlib)

# miexed effects stuff
library(nlme)
library(lme4)

library(lmtest)
library(influence.ME)
library(lmerTest)

library(patchwork)
library(stringr)
library(data.table)


# library(glmmTMB)
library(mvtnorm)
library(numDeriv)
library(LaplacesDemon)
library(TMB)

# Reading data -----------------------------------------------------------------
dat <- read.table("clothingFullAss03.csv", sep=',', head=TRUE)
df <- data.frame(dat) %>% select(-X) 
cols <- c("sex", "subjId", "day", "subDay")
df[cols] <- lapply(df[cols], as.factor)
summary(df)
sum(is.na(df)) # no missing values
df["sex_optimization"] = as.numeric(df$sex == "male")
attach(df)

# Model_1 - likelihood implementation -----------------------------------------
## regular likelihood - multivariate normal ------------------------------------
fit_lm = lm(y ~ X)
fit_lm_coefficeitns = as.numeric(fit_lm$coefficients[c(1, 3)])


n = dim(df)[1]
X = cbind(rep(1, n), df$sex_optimization)
p = dim(X)[2]
Z = cbind(rep(1, n))
u_number = dim(Z)[2]
y = log(clo)

obj_1 = function(beta){
  result = 0
  for(subject_i in unique(df$subjId)){
    X_i = X[df$subjId == subject_i, , drop = F]
    n_i = dim(X_i)[1]
    y_i = y[df$subjId == subject_i]
    y_hat_i = as.numeric(X_i %*% beta[1:p])
    Psi_i = exp(beta[p+1])*diag(u_number)
    Z_i = Z[df$subjId == subject_i, , drop = F]
    Sigma_full_i = exp(beta[p+2])*diag(n_i) + Z_i %*% Psi_i %*% t(Z_i)
    result = result +
      sum(dmvnorm(y_i, mean = y_hat_i, sigma = Sigma_full_i, log = TRUE))
  }
  return(-result)
}

beta_initial1 = c(rep(1, p), "sigma^2_log" = log(1), "sigma_u^2_log" = log(1))
# sanity check
obj_1(beta_initial1)

opt_1 <- nlminb(beta_initial1, obj_1)
opt_1$par
exp(opt_1$par[3:4])

hessian_1 = hessian(obj_1, opt_1$par)
sqrt(diag(solve(hessian_1)))

# reference
fit1 = lmer(log(clo) ~ as.factor(sex_optimization) +
              (1|subjId), REML = F)
summary(fit1)

## TMB - Laplace approximation -------------------------------------------------
### Hierarchical appraoch ------------------------------------------------------
n = dim(df)[1]
X = cbind(rep(1, n), df$sex_optimization)
p = dim(X)[2]
Z = cbind(rep(1, n))
u_number = dim(Z)[2]
y = log(clo)
n_subjects = length(unique(df$subjId))

df_subject = data.frame(table(df$subjId))
subjectId_frequencies = df_subject$Freq

data <- list(y = y, sex = X[, 2],
             subjectId_factor = as.factor((df$subjId)),
             subjectId_frequencies = subjectId_frequencies)

# data <- list(y = y, mu = X[, 1], sex = X[, 2],
#              subjectId_factor = as.factor(as.numeric(df$subjId)))

compile("clothing_1.cpp")
dyn.load(dynlib("clothing_1"))
is.loaded("clothing_1")

# parameters <- list(u = rep(1, n_subjects),
#                    beta = opt_1$par)

parameters <- list(u = rep(1, n_subjects),
                   beta = c(fit_lm_coefficeitns, -1, -1))

obj_test <- MakeADFun(data = data,
                 parameters = parameters,
                 random = "u",
                 DLL = "clothing_1")


## Test eval function and gradient
obj$fn()
obj$gr()

## Fit model
opt <- nlminb(obj_test$par, obj_test$fn, obj_test$gr)
opt$par
opt$objective

c(opt$par, opt_1$par)
# WORKS PERFECTLY !!! 
opt$par
# opt_1$par

## report on result
sdreport(obj_test)

### Multivariate distribution approach -----------------------------------------
# range = seq(100)
# 
# n = dim(df)[1]
# X = cbind(rep(1, n), df$sex_optimization)[range, ]
# p = dim(X)[2]
# Z = cbind(rep(1, n))[range, ]
# u_number = dim(Z)[2]
# y = log(clo)
# n_subjects = length(unique(df$subjId))
# 
# df_subject = data.frame(table(df$subjId))
# subjectId_frequencies = df_subject$Freq
# 
# compile("clothing_1_multivariate.cpp")
# dyn.load(dynlib("clothing_1_multivariate"))
# 
# data <- list(y = y, sex = X[, 2],
#              subjectId_factor = as.factor((df$subjId))[range],
#              subjectId_frequencies = subjectId_frequencies[range],
#              nsubjects_size = length(subjectId_frequencies))
# 
# # parameters <- list(beta = fit_lm_coefficeitns,
# #                    sigma2_u_log = -1,
# #                    sigma2_log = -1)
# 
# parameters <- list(beta = c(0, 0),
#                    sigma2_u_log = -1,
#                    sigma2_log = -1)
# 
# obj <- MakeADFun(data = data,
#                  parameters = parameters,8
#                  DLL = "clothing_1_multivariate")
# 
# 
# ## Test eval function and gradient
# obj$fn()
# obj$gr()
# 
# ## Fit model
# opt <- nlminb(obj$par, obj$fn, obj$gr)
# opt$par
# opt$objective
# 
# c(opt$par, opt_1$par)
# # WORKS PERFECTLY !!!
# opt$par
# opt_1$par
# 
# ## report on result
# sdreport(obj)

# Model_2 - likelihood implementation -----------------------------------------
## regular likelihood - multivariate normal ------------------------------------
n = dim(df)[1]
X = cbind(rep(1, n), df$sex_optimization)
p = dim(X)[2]
Z1 = cbind(rep(1, n))
u_number1 = dim(Z1)[2]
Z2 = cbind(rep(1, n))
u_number2 = dim(Z2)[2]
y = log(clo)

obj_2 = function(beta){
  result = 0
  for(subject_i in unique(df$subjId)){
    X_i = X[df$subjId == subject_i, , drop = F]
    n_i = dim(X_i)[1]
    y_i = y[df$subjId == subject_i]
    y_hat_i = as.numeric(X_i %*% beta[1:p])
    
    Psi1_i = exp(beta[p+1])*diag(u_number)
    Z1_i = Z1[df$subjId == subject_i, , drop = F]
    Z2_i = Z2[df$subjId == subject_i, , drop = F]
    
    subject_i_days_length = length(unique(df$day[df$subjId == subject_i]))
    day_matrix_list = vector("list", subject_i_days_length)
    for(index_day in seq(subject_i_days_length)){
      subject_i_days = df$day[df$subjId == subject_i]
      day_j = unique(subject_i_days)[index_day]
      Psi2_ij_tmp = exp(beta[p+3])*diag(u_number2)
      Z2_ij = Z2_i[subject_i_days == day_j, , drop = F]
      matrix_tmp = Z2_ij %*% Psi2_ij_tmp %*% t(Z2_ij)
      day_matrix_list[[index_day]] = matrix_tmp
    }
    Psi2_ij = as.matrix(bdiag(day_matrix_list))
    Sigma_full_i = exp(beta[p+2])*diag(n_i) +
      Z1_i %*% Psi1_i %*% t(Z1_i) + Psi2_ij
    
    result = result +
      sum(dmvnorm(y_i, mean = y_hat_i, sigma = Sigma_full_i, log = TRUE))
  }
  return(-result)
}

beta_initial2 = c(opt_1$par, "sigma_v^2_log" = log(1))
# sanity check
obj_2(beta_initial2)

opt_2 <- nlminb(beta_initial2, obj_2)
opt_2$par
opt_2$objective

hessian_2 = hessian(obj_2, opt_2$par)
sqrt(diag(solve(hessian_2)))
exp(opt_2$par[3:5])

# reference
fit1 = lmer(log(clo) ~ as.factor(sex_optimization) +
              (1|subjId) + (1|subjId:day), REML = F)
plot(fit1)
# summary(fit1)
summary(fit1)

resid_male = residuals(fit1)[df$sex == "male"]
resid_female = residuals(fit1)[df$sex == "female"]
c(var(resid_male), var(resid_female), var(resid_female)/var(resid_male))

## TMB - Laplace approximation -------------------------------------------------
range = seq(n)
# the R session is aborted for more observations 
# the optimization is over 46 + 136 random parameters, quite a lot

n = dim(df)[1]
X = cbind(rep(1, n), df$sex_optimization)[range, ]
p = dim(X)[2]
Z = cbind(rep(1, n))[range, ]
u_number = dim(Z)[2]
y = log(clo)[range]
n_subjects = length(unique(df$subjId))
n_days = length(unique(df$subDay))

data <- list(y = y, sex = X[, 2],
             subjectId_factor = as.factor(df$subjId)[range],
             subjectId_day_factor = as.factor(df$subDay)[range]
             )

compile("clothing_2.cpp")
dyn.load(dynlib("clothing_2"))

# parameters <- list(u = rep(1, n_subjects),
#                    beta = opt_1$par)

parameters <- list(u = rep(0, n_subjects),
                   v = rep(0, n_days),
                   beta = fit_lm_coefficeitns,
                   sigma2_u_log = log(1),
                   sigma2_v_log = log(1),
                   sigma2_log = log(1))

obj_2_Laplace <- MakeADFun(data = data,
                 parameters = parameters,
                 random = c("u","v"),
                 DLL = "clothing_2")


## Test eval function and gradient
obj_2_Laplace$fn()
obj_2_Laplace$gr()

## Fit model
opt <- nlminb(obj_2_Laplace$par, obj_2_Laplace$fn,
                        obj_2_Laplace$gr)
opt$par
opt$objective

# c(opt$par, opt_1$par)
# # WORKS PERFECTLY !!! 

## report on result
rap <- sdreport(obj_2_Laplace,getJointPrecision = TRUE)
print(rap)

# Model_3 - likelihood implementation -----------------------------------------
# one alpha as a weight for each variance parameter 
## regular likelihood - multivariate normal ------------------------------------
n = dim(df)[1]
X = cbind(rep(1, n), df$sex_optimization)
p = dim(X)[2]
Z1 = cbind(rep(1, n))
u_number1 = dim(Z1)[2]
Z2 = cbind(rep(1, n))
u_number2 = dim(Z2)[2]
X_alphas = X # note that it's the case in our model
y = log(clo)

obj_3 = function(beta){
  
  result = 0
  for(subject_i in unique(df$subjId)){
    X_i = X[df$subjId == subject_i, , drop = F]
    n_i = dim(X_i)[1]
    y_i = y[df$subjId == subject_i]
    y_hat_i = as.numeric(X_i %*% beta[1:p])
    
    X_alphas_i = unique(X_alphas[df$subjId == subject_i, , drop = F])
    Psi1_i = as.numeric(
      exp(X_alphas_i %*% beta[c("sigma_u^2_log", "alpha")])
    )*diag(u_number1)
    
    Z1_i = Z1[df$subjId == subject_i, , drop = F]
    Z2_i = Z2[df$subjId == subject_i, , drop = F]

    subject_i_days_length = length(unique(df$day[df$subjId == subject_i]))
    day_matrix_list = vector("list", subject_i_days_length)
    for(index_day in seq(subject_i_days_length)){
      subject_i_days = df$day[df$subjId == subject_i]
      day_j = unique(subject_i_days)[index_day]
      X_alphas_ij = X_alphas_i
      Psi2_ij_tmp = as.numeric(
        exp(X_alphas_ij %*% beta[c("sigma_v^2_log", "alpha")])
      )*diag(u_number2)
      Z2_ij = Z2_i[subject_i_days == day_j, , drop = F]
      
      matrix_tmp = Z2_ij %*% Psi2_ij_tmp %*% t(Z2_ij)
      day_matrix_list[[index_day]] = matrix_tmp
    }
    Psi2_ij = as.matrix(bdiag(day_matrix_list))
    
    Sigma_i = as.numeric(
      exp(X_alphas_i %*% beta[c("sigma^2_log", "alpha")])
    )*diag(n_i)
    
    Sigma_full_i = Sigma_i + Z1_i %*% Psi1_i %*% t(Z1_i) + Psi2_ij
    
    result = result +
      sum(dmvnorm(y_i, mean = y_hat_i, sigma = Sigma_full_i, log = TRUE))
  }
  return(-result)
}

beta_initial3 = c(opt_2$par, "alpha" = log(1))

# sanity check
obj_3(beta_initial3)

opt_3 <- nlminb(beta_initial3, obj_3)
opt_3$par
exp(opt_3$par["alpha"])
opt_3$objective

hessian_3 = hessian(obj_3, opt_3$par)
sqrt(diag(solve(hessian_3)))
exp(opt_3$par[3:6])

### Functions
# spline()
# splinefun()
# uniroot()

## TMB - Laplace approximation -------------------------------------------------
range = seq(n)
n = dim(df)[1]
X = cbind(rep(1, n), df$sex_optimization)[range, ]
p = dim(X)[2]
Z = cbind(rep(1, n))[range, ]
u_number = dim(Z)[2]
y = log(clo)[range]
n_subjects = length(unique(df$subjId))
n_days = length(unique(df$subDay))

data <- list(y = y, sex = X[, 2],
             subjectId_factor = as.factor(df$subjId)[range],
             subjectId_day_factor = as.factor(df$subDay)[range]
             )

compile("clothing_3.cpp")
dyn.load(dynlib("clothing_3"))

parameters <- list(u = rep(0, n_subjects),
                   v = rep(0, n_days),
                   beta = fit_lm_coefficeitns,
                   alpha = log(1),
                   sigma2_u_log = log(1),
                   sigma2_v_log = log(1),
                   sigma2_log = log(1))

obj_3_Laplace <- MakeADFun(data = data,
                           parameters = parameters,
                           random = c("u","v"),
                           DLL = "clothing_3")


## Test eval function and gradient
obj_3_Laplace$fn()
obj_3_Laplace$gr()

## Fit model
opt_3_Laplace <- nlminb(obj_3_Laplace$par, obj_3_Laplace$fn,
              obj_3_Laplace$gr)
opt_3_Laplace$par
opt_3_Laplace$objective

# c(opt$par, opt_1$par)
# # WORKS PERFECTLY !!! 

## report on result
rap_3 <- sdreport(obj_3_Laplace, getJointPrecision = TRUE)
print(rap_3)  


# Model_5 - likelihood implementation -----------------------------------------
# one alpha as a weight for each variance parameter 
## regular likelihood - multivariate normal ------------------------------------
n = dim(df)[1]
X = cbind(rep(1, n), df$sex_optimization)
p = dim(X)[2]
Z1 = cbind(rep(1, n))
u_number1 = dim(Z1)[2]
Z2 = cbind(rep(1, n))
u_number2 = dim(Z2)[2]
X_alphas = X # note that it's the case in our model
y = log(clo)

obj_5 = function(beta){
  result = 0
  for(subject_i in unique(df$subjId)){
    X_i = X[df$subjId == subject_i, , drop = F]
    n_i = dim(X_i)[1]
    y_i = y[df$subjId == subject_i]
    y_hat_i = as.numeric(X_i %*% beta[1:p])
    
    X_alphas_i = unique(X_alphas[df$subjId == subject_i, , drop = F])
    Psi1_i = as.numeric(
      exp(X_alphas_i %*% beta[c("sigma_u^2_log", "alpha")])
    )*diag(u_number1)
    
    Z1_i = Z1[df$subjId == subject_i, , drop = F]
    Z2_i = Z2[df$subjId == subject_i, , drop = F]
    
    subject_i_days_length = length(unique(df$day[df$subjId == subject_i]))
    day_matrix_list = vector("list", subject_i_days_length)
    for(index_day in seq(subject_i_days_length)){
      subject_i_days = df$day[df$subjId == subject_i]
      day_j = unique(subject_i_days)[index_day]
      X_alphas_ij = X_alphas_i
      Psi2_ij_tmp = as.numeric(
        exp(X_alphas_ij %*% beta[c("sigma_v^2_log", "alpha")])
      )*diag(u_number2)
      Z2_ij = Z2_i[subject_i_days == day_j, , drop = F]
      
      matrix_tmp = Z2_ij %*% Psi2_ij_tmp %*% t(Z2_ij)
      day_matrix_list[[index_day]] = matrix_tmp
    }
    Psi2_ij = as.matrix(bdiag(day_matrix_list))
    
    Sigma_i = as.numeric(
      exp(X_alphas_i %*% beta[c("sigma^2_log", "alpha")])
    )*diag(n_i)
    Sigma_full_i = Sigma_i + Z1_i %*% Psi1_i %*% t(Z1_i) + Psi2_ij
    
    result = result +
      sum(dmvt(y_i, y_hat_i, sigma = as.matrix(Sigma_full_i),
               df = 2*exp(beta["lambda_log"]),
               log = TRUE))
  }
  return(-result)
}

beta_initial5= c(opt_3$par, "lambda_log" = log(1))

# sanity check
obj_5(beta_initial5)

 

opt_5 <- nlminb(beta_initial5, obj_5)

opt_5$par
opt_3$par
exp(opt_5$par["lambda_log"])

opt_5$objective

hessian_5 = hessian(obj_5, opt_5$par)
sqrt(diag(solve(hessian_5)))
exp(opt_5$par)

## Wishart distribution - gamma -------------------------------------------------

## TMB - Laplace approximation -------------------------------------------------

range = seq(n)
n = dim(df)[1]
X = cbind(rep(1, n), df$sex_optimization)[range, ]
p = dim(X)[2]
Z = cbind(rep(1, n))[range, ]
u_number = dim(Z)[2]
y = log(clo)[range]
n_subjects = length(unique(df$subjId))
n_days = length(unique(df$subDay))

subjectId_factor = as.factor(df$subjId)[range]
subjectId_day_factor = as.factor(df$subDay)[range]

# testing stuff
X_test = df[, c("subjId", "subDay")]
subjectId_day_factor_gamma = as.factor(unique(X_test)[, 1])

data <- list(y = y, sex = X[, 2],
             subjectId_factor = subjectId_factor,
             subjectId_day_factor = subjectId_day_factor,
             subjectId_day_factor_gamma = subjectId_day_factor_gamma)


compile("clothing_5_hierarchical.cpp")
dyn.load(dynlib("clothing_5_hierarchical"))

parameters <- list(u = rep(0, n_subjects),
                   v = rep(0, n_days),
                   gamma = rep(1, n_subjects),
                   beta = fit_lm_coefficeitns,
                   alpha = log(1),
                   sigma2_u_log = log(1),
                   sigma2_v_log = log(1),
                   lambda = log(1),
                   sigma2_log = log(1))

obj_5_Laplace <- MakeADFun(data = data,
                           parameters = parameters,
                           random = c("u", "v", "gamma"),
                           DLL = "clothing_5_hierarchical")

compile("clothing_5_hierarchical_regular.cpp")
dyn.load(dynlib("clothing_5_hierarchical_regular"))

parameters <- list(u = rep(0, n_subjects),
                   v = rep(0, n_days),
                   gamma = rep(1, n_subjects),
                   beta = fit_lm_coefficeitns,
                   alpha = log(1),
                   sigma2_u_log = log(1),
                   sigma2_v_log = log(1),
                   lambda = log(1),
                   sigma2_log = log(1))

obj_5_Laplace <- MakeADFun(data = data,
                           parameters = parameters,
                           random = c("u", "v", "gamma"),
                           DLL = "clothing_5_hierarchical_regular")

## Test eval function and gradient
obj_5_Laplace$fn()
obj_5_Laplace$gr()

## Fit model
obj <- nlminb(obj_5_Laplace$par, obj_5_Laplace$fn,
              obj_5_Laplace$gr)

obj$par
obj$objective

opt_5$par

## report on result
rap_5 <- sdreport(obj_5_Laplace)
print(rap_5)

rap$diag.cov.random
names(rap)
rap$par.fixed

# Model_6 - Laplace approximation -----------------------------------------
# one alpha as a weight for each variance parameter 

range = seq(n)

n = dim(df)[1]
X = cbind(rep(1, n), df$sex_optimization)[range, ]
p = dim(X)[2]
Z = cbind(rep(1, n))[range, ]
u_number = dim(Z)[2]
y = log(clo)[range]
n_subjects = length(unique(df$subjId))
n_days = length(unique(df$subDay))

# subjectId_factor = as.factor(as.numeric(df$subjId))
# subjectId_day_factor = as.factor(as.numeric(df$subDay))
# 
# # testing stuff
# X_test = df[, c("subjId", "subDay")]
# unique(X_test) # R is just to good :p
# subjectId_day_factor_gamma = as.factor(as.numeric(unique(X_test)[, 1]))

subjectId_factor = as.factor(df$subjId)[range]
subjectId_day_factor = as.factor(df$subDay)[range]

# testing stuff
X_test = df[, c("subjId", "subDay")]
# unique(X_test) # R is just to good :p
subjectId_day_factor_gamma = as.factor(unique(X_test)[, 1])[range]

data <- list(y = y, sex = X[, 2],
             subjectId_factor = subjectId_factor,
             subjectId_day_factor = subjectId_day_factor,
             subjectId_day_factor_gamma = subjectId_day_factor_gamma)


# compile("clothing_6.cpp")
dyn.load(dynlib("clothing_6"))


# dyn.load("clothing_6.dll")
# is.loaded(dynlib("clothing_6"))
# 
# is.loaded("clothing_6")
# getLoadedDLLs()
# 
# dyn.unload(dynlib("clothing6")) # First attempt
# gc() # Garbage Collection
# dyn.unload(dynlib("clothing6")) # Second attempt
# getLoadedDLLs() #Verify that second attempt works
# 
# -John

# compile("clothing6testing.cpp")
# dyn.load(dynlib("clothing6testing"))

parameters <- list(u = rep(0, n_subjects),
                   v = rep(0, n_days),
                   gamma = rep(0, n_subjects),
                   beta = fit_lm_coefficeitns,
                   alpha = log(1),
                   sigma2_u_log = log(1),
                   sigma2_v_log = log(1),
                   sigma2_G_log = log(1),
                   sigma2_log = log(1)
                   )

obj_6_Laplace <- MakeADFun(data = data,
                           parameters = parameters,
                           random = c("u", "v", "gamma"),
                           DLL = "clothing_6"
                           )

## Test eval function and gradient
obj_6_Laplace$fn()
obj_6_Laplace$gr()

## Fit model
opt_6 <- nlminb(obj_6_Laplace$par, obj_6_Laplace$fn,
              obj_6_Laplace$gr)
opt_6$par
opt_6$objective

## report on result
rap_6 <- sdreport(obj_6_Laplace)
print(rap_6)


rap$par.random
rap$diag.cov.random
names(rap)
rap$par.fixed

# 7 - Comparison of the posterior distributions --------------------------------
## estimation of random effects in t-distribution hierarchical model -----------

random_normal_gamma = rap_6$par.random[names(rap_6$par.random) == "gamma"]
random_hierarchical_gamma = rap_5$par.random[names(rap_5$par.random) == "gamma"]

par(mfrow = c(1, 2))
plot(exp(random_normal), col = "red", type = "p")
lines(random_hierarchical, col = "blue", type = "p")
plot(exp(random_normal), random_hierarchical)
abline(0, 1)

n_days = length(unique(df$subDay))

# 136 fitted values
mean_fun_normal <- function(){
  u = rap_6$par.random[names(rap_6$par.random) == "u"]
  v = rap_6$par.random[names(rap_6$par.random) == "v"]
  beta = rap_6$par.fixed
  result = numeric(n_days)
  index = 1
  for(subject_i in unique(df$subjId)){
    subject_i_days_length = length(unique(df$day[df$subjId == subject_i]))
    for(index_day in seq(subject_i_days_length)){
      result[index] = beta[1] + df$sex_optimization[index]*beta[2] +
        u[as.numeric(subject_i) + 1] + v[as.numeric(index_day) + 1]
      index = index + 1
    }
  }
  return(result)
}

mean_fun_normal_3 <- function(){
  u = rap_3$par.random[names(rap_3$par.random) == "u"]
  v = rap_3$par.random[names(rap_3$par.random) == "v"]
  beta = rap_3$par.fixed
  result = numeric(n_days)
  index = 1
  for(subject_i in unique(df$subjId)){
    subject_i_days_length = length(unique(df$day[df$subjId == subject_i]))
    print(subject_i)
    for(index_day in seq(subject_i_days_length)){
      result[index] = beta[1] + df$sex_optimization[index]*beta[2] +
        u[as.numeric(subject_i) + 1] + v[as.numeric(index_day)]
      index = index + 1
    }
    result[index] = beta[1] + df$sex_optimization[index]*beta[2] +
      u[as.numeric(subject_i) + 1] + v[as.numeric(index_day)]
    index = index + 1
    
  }
  return(result)
}

X = cbind(rep(1, n), df$sex_optimization)
mean_fun_hierarchical <- function(){
  u = rap_5$par.random[names(rap_5$par.random) == "u"]
  v = rap_5$par.random[names(rap_5$par.random) == "v"]
  u_all = rep(u, indexing_helper_subjects)
  v_all = rep(u, indexing_helper_days)
  result = numeric(n_days)
  index = 1
  result = X %*% beta[1:2] + u_all + v_all
  return(result)
  # for(subject_i in unique(df$subjId)){
  #   subject_i_days_length = length(unique(df$day[df$subjId == subject_i]))
  #   for(index_day in seq(subject_i_days_length)){
  #     result[index] = beta[1] + df$sex_optimization[index]*beta[2] +
  #       u[as.numeric(subject_i) + 1] + v[as.numeric(index_day) + 1]
  #     index = index + 1
  #   }
  # }
  return(result)
}

mean_fun <- function(rap_tmp){
  u = rap_tmp$par.random[names(rap_tmp$par.random) == "u"]
  v = rap_tmp$par.random[names(rap_tmp$par.random) == "v"]
  u_all = rep(u, indexing_helper_subjects)
  v_all = rep(v, indexing_helper_days)
  
  beta = rap_tmp$par.fixed
  result = numeric(n_days)
  index = 1
  result = X %*% beta[1:2] + u_all + v_all
  return(result)
}

indexing_helper_subjects = data.frame(table(df$subjId))$Freq
indexing_helper_days = data.frame(table(df$subDay))$Freq

# fitted_hierarchical = mean_fun_hierarchical()
# fitted_normal = mean_fun_normal()
# fitted_normal_3 = mean_fun_normal_3()

fitted_hierarchical = mean_fun(rap_5)
fitted_normal = mean_fun(rap_6)
fitted_normal_3 = mean_fun(rap_3)

# fitted_hierarchical_all = rep(fitted_hierarchical, indexing_helper_days)
# fitted_normal_all = rep(fitted_normal, indexing_helper_days)
# fitted_normal_3_all = rep(fitted_normal_3, indexing_helper_days)

resid_hierarchical = y - fitted_hierarchical
resid_normal = y - fitted_normal
resid_normal_3 = y - fitted_normal_3

gamma_5 = rap_5$par.random[names(rap_5$par.random) == "gamma"]
gamma_6 = rap_6$par.random[names(rap_6$par.random) == "gamma"]

var_hierarchical = exp(rap_5$par.fixed["sigma2_log"]) *
  exp(df$sex_optimization * rap_5$par.fixed["alpha"]) /
  rep(gamma_5, indexing_helper_subjects)
var_normal = exp(rap_6$par.fixed["sigma2_log"]) *
  exp(df$sex_optimization * rap_6$par.fixed["alpha"]) *
  exp(-rep(gamma_6, indexing_helper_subjects))
var_normal_3 = exp(rap_3$par.fixed["sigma2_log"]) *
  exp(df$sex_optimization * rap_3$par.fixed["alpha"])

resid_hierarchical_standarized = resid_hierarchical/sqrt(var_hierarchical)
resid_normal_standarized = resid_normal/sqrt(var_normal)
resid_normal_3_standarized = resid_normal_3/sqrt(var_normal_3)


#### Comparing the variances
vertical_lines = c(240)

par(mfrow = c(3, 1))
plot(var_hierarchical, col = "red", type = "p")
lines(var_normal, col = "blue", type = "p")
lines(var_normal_3, col = "green", type = "p")
legend("topright",
       legend=c("Hierarchical", "Normal", "Normal_3"),
       col=c("red", "blue", "green"), pch = c(1,1))
abline(v = vertical_lines)

plot(resid_hierarchical, col = "red", type = "p")
lines(resid_normal, col = "blue", type = "p")
lines(resid_normal_3, col = "green", type = "p")
legend("topright",
       legend=c("Hierarchical", "Normal", "Normal_3"),
       col=c("red", "blue", "green"), pch = c(1,1))
abline(v = vertical_lines)

plot(resid_hierarchical_standarized, col = "red", type = "p")
lines(resid_normal_standarized, col = "blue", type = "p")
lines(resid_normal_3_standarized, col = "green", type = "p")
# lines(resid(fit1), col = "orange", type = "p")
legend("topright",
       legend=c("Hierarchical", "Normal", "Normal_3"),
       col=c("red", "blue", "green"), pch = c(1,1))
abline(v = vertical_lines)

par(mfrow = c(1, 3))
plot(fitted_hierarchical, col = "red", type = "p")
lines(fitted_normal, col = "blue", type = "p")
lines(fitted_normal_3, col = "green", type = "p")
legend("topright",
       legend=c("Hierarchical", "Normal", "Normal_3"),
       col=c("red", "blue", "green"), pch = c(1,1))

plot(resid_hierarchical, col = "red", type = "p")
lines(resid_normal, col = "blue", type = "p")
lines(resid_normal_3, col = "green", type = "p")
legend("topright",
       legend=c("Hierarchical", "Normal", "Normal_3"),
       col=c("red", "blue", "green"), pch = c(1,1))

# par(mfrow = c(1,1))
plot(resid_hierarchical_standarized, col = "red", type = "p")
lines(resid_normal_standarized, col = "blue", type = "p")
lines(resid_normal_3_standarized, col = "green", type = "p")
# lines(resid(fit1), col = "orange", type = "p")
legend("topright",
       legend=c("Hierarchical", "Normal", "Normal_3"),
       col=c("red", "blue", "green"), pch = c(1,1))

## -----------------------------------------------------------------------------
par(mfrow = c(1, 1))
plot(resid_hierarchical, col = "red", type = "p")
lines(resid_normal, col = "blue", type = "p")
lines(resid_normal_3, col = "green", type = "p")
legend("topright",
       legend=c("Hierarchical", "Normal", "Normal_3"),
       col=c("red", "blue", "green"), pch = c(1,1))


par(mfrow = c(1, 1))
plot(resid_normal_3, col = "green", type = "p")
lines(resid_hierarchical, col = "red", type = "p")
legend("topright",
       legend=c("fit1", "Hierarchical"),
       col=c("green", "red"), pch = c(1,1))

names(rap)
rap$par.fixed
mean.fun(rap$par.fixed[1:4],rap$par.random[1],
         data$t[1:7],data$season[1:7])

################################################################################
model = lmer(log(clo) ~ sex * poly(tInOp, 2) * poly(tOut, 2) +
                (1 + poly(tInOp, 2) + poly(tOut, 2) | subjId/subDay))


