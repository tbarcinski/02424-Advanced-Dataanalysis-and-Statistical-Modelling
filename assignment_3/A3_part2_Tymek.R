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

library(TMB)

# Reading data -----------------------------------------------------------------
dat <- read.table("clothingFullAss03.csv", sep=',', head=TRUE)
df <- data.frame(dat) %>% select(-X) 
cols <- c("sex", "subjId", "day", "subDay")
df[cols] <- lapply(df[cols], as.factor)
summary(df)
sum(is.na(df)) # no missing values
df["sex_optimization"] = as.numeric(df$sex == "female")
attach(df)

# Model_1 - likelihood implementation -----------------------------------------
## regular likelihood - multivariate normal ------------------------------------
n = dim(df)[1]
X = cbind(rep(1, n), df$sex_optimization)
p = dim(X)[2]
Z = cbind(rep(1, n))
u_number = dim(Z)[2]
y = log(clo)

fit_lm = lm(y ~ X)
fit_lm_coefficeitns = as.numeric(fit_lm$coefficients[c(1, 3)])


obj_1 = function(beta){
  result = 0
  for(subject_i in unique(df$subjId)){
    X_i = X[df$subjId == subject_i, , drop = F]
    n_i = dim(X_i)[1]
    y_i = y[df$subjId == subject_i]
    y_hat_i = X_i %*% beta[1:p]
    Psi_i = exp(beta[p+1])*diag(u_number)
    Z_i = Z[df$subjId == subject_i, , drop = F]
    Sigma_full_i = exp(beta[p+2])*diag(n_i) + Z_i %*% Psi_i %*% t(Z_i)
    result = result +
      sum(dmvnorm(y_i, mean = y_hat_i, sigma = Sigma_full_i, log = TRUE))
  }
  return(-result)
}

beta_initial1 = c(rep(1, p), "sigma^2_log" = log(1), "sigma_subject^2_log" = log(1))
# sanity check
obj_1(beta_initial1)

opt_1 <- nlminb(beta_initial1, obj_1)
opt_1$par
exp(opt_1$par)

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

obj <- MakeADFun(data = data,
                 parameters = parameters,
                 random = "u",
                 DLL = "clothing_1")


## Test eval function and gradient
obj$fn()
obj$gr()

## Fit model
opt <- nlminb(obj$par, obj$fn, obj$gr)
opt$par
opt$objective

c(opt$par, opt_1$par)
# WORKS PERFECTLY !!! 
opt$par
# opt_1$par

## report on result
sdreport(obj)

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
    y_hat_i = X_i %*% beta[1:p]
    
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

beta_initial2 = c(opt_1$par, "sigma_day^2_log" = log(1))
# sanity check
obj_2(beta_initial2)

opt_2 <- nlminb(beta_initial2, obj_2)
opt_2$par
opt_2$objective

# reference
fit1 = lmer(log(clo) ~ as.factor(sex_optimization) +
              (1|subjId) + (1|subjId:day), REML = F)
# summary(fit1)
fit1

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
obj_2_Laplace <- nlminb(obj_2_Laplace$par, obj_2_Laplace$fn,
                        obj_2_Laplace$gr)
obj_2_Laplace$par
obj_2_Laplace$objective

# c(opt$par, opt_1$par)
# # WORKS PERFECTLY !!! 

## report on result
rap <- sdreport(obj_2_Laplace,getJointPrecision = TRUE)
sdreport(obj_2_Laplace)


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

X_alphas = X # note that it's the case in our model formulation

y = log(clo)

obj_3 = function(beta){
  
  result = 0
  for(subject_i in unique(df$subjId)){
    X_i = X[df$subjId == subject_i, , drop = F]
    n_i = dim(X_i)[1]
    
    y_i = y[df$subjId == subject_i]
    y_hat_i = X_i %*% beta[1:p]
    
    X_alphas_i = unique(X_alphas[df$subjId == subject_i, , drop = F])
    Psi1_i = as.numeric(
      exp(X_alphas_i %*% beta[c("sigma_subject^2_log", "alpha")])
    )*diag(u_number)
    
    Z1_i = Z1[df$subjId == subject_i, , drop = F]
    Z2_i = Z2[df$subjId == subject_i, , drop = F]
    
    subject_i_days_length = length(unique(df$day[df$subjId == subject_i]))
    day_matrix_list = vector("list", subject_i_days_length)
    for(index_day in seq(subject_i_days_length)){
      subject_i_days = df$day[df$subjId == subject_i]
      day_j = unique(subject_i_days)[index_day]
      
      X_alphas_ij = X_alphas_i
      Psi2_ij_tmp = as.numeric(
        exp(X_alphas_ij %*% beta[c("sigma_day^2_log", "alpha")])
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

beta_initial3 = c(opt_2$par, "alpha" = 1)

# sanity check
obj_3(beta_initial3)

opt_3 <- nlminb(beta_initial3, obj_3)
opt_3$par
exp(opt_3$par["alpha"])

### Functions
# spline()
# splinefun()
# uniroot()

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

X_alphas = X # note that it's the case in our model formulation

y = log(clo)

obj_5 = function(beta){
  
  result = 0
  for(subject_i in unique(df$subjId)){
    X_i = X[df$subjId == subject_i, , drop = F]
    n_i = dim(X_i)[1]
    
    y_i = y[df$subjId == subject_i]
    y_hat_i = X_i %*% beta[1:p]
    
    X_alphas_i = unique(X_alphas[df$subjId == subject_i, , drop = F])
    Psi1_i = as.numeric(
      exp(X_alphas_i %*% beta[c("sigma_subject^2_log", "alpha")])
    )*diag(u_number)
    
    Z1_i = Z1[df$subjId == subject_i, , drop = F]
    Z2_i = Z2[df$subjId == subject_i, , drop = F]
    
    subject_i_days_length = length(unique(df$day[df$subjId == subject_i]))
    day_matrix_list = vector("list", subject_i_days_length)
    for(index_day in seq(subject_i_days_length)){
      subject_i_days = df$day[df$subjId == subject_i]
      day_j = unique(subject_i_days)[index_day]
      
      X_alphas_ij = X_alphas_i
      Psi2_ij_tmp = as.numeric(
        exp(X_alphas_ij %*% beta[c("sigma_day^2_log", "alpha")])
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
      sum(dmvt(y_i, delta = y_hat_i, type = "shifted",
               sigma = Sigma_full_i,
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


compile("clothing_6.cpp")
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
obj_6_Laplace <- nlminb(obj_6_Laplace$par, obj_6_Laplace$fn,
                        obj_6_Laplace$gr)
obj_6_Laplace$par
obj_6_Laplace$objective

## report on result
rap <- sdreport(obj_6_Laplace)
rap$par.random
rap$diag.cov.random
names(rap)
rap$par.fixed

dlls <- getLoadedDLLs()
isTMBdll <- function(dll)!is(try(getNativeSymbolInfo("MakeADFunObject",dll),TRUE),"try-error")
TMBdll <- sapply(dlls, isTMBdll)
if(sum(TMBdll) == 0) stop("There are no TMB models loaded (use 'dyn.load').")
if(sum(TMBdll) >1 ) stop("Multiple TMB models loaded. Failed to guess DLL name.")
names(dlls[TMBdll])

gc()



