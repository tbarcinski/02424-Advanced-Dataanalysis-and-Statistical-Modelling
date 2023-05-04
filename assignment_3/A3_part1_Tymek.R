# File Description -------------------------------------------------------------
#
#   Tymoteusz Barcinski - s221937
#   
#   Advanced Dataanalysis and Statistical Modelling
#   Assignment 3 - Mixed effects and hierarchical models
#   Part1 - Mixed effects models
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

# Reading data -----------------------------------------------------------------
dat <- read.table("clothingFullAss03.csv", sep=',', head=TRUE)
df <- data.frame(dat) %>% select(-X) 
cols <- c("sex", "subjId", "day", "subDay")
df[cols] <- lapply(df[cols], as.factor)
summary(df)
sum(is.na(df)) # no missing values
attach(df)

# Linear model from assignment 2 ----------------------------------------------
linear_model_formula = log(clo) ~ poly(tInOp, 2) + poly(tOut, 2) +
  sex + poly(tInOp, 1):poly(tOut, 1)
mlog = lm(linear_model_formula)
summary(mlog)
drop1(mlog, test = "F")

par(mfrow = c(2, 2))
plot(mlog)
shapiro.test(mlog$residuals)

# problems with residuals - do not follow a normal distribution
# skewed distribution

# Simple mixed effect model ----------------------------------------------------
plot_residuals_nlme = function(model){
  
  par(mfrow = c(2, 2))
  residuals_tmp = residuals(model) / sd(residuals(model))
  
  # plot of the residuals themselves
  plot(residuals_tmp)
  
  # QQ plot for the residuals
  qqnorm(residuals_tmp, main = "Q-Q Plot for the residuals")
  qqline(residuals_tmp, col = "steelblue", lwd = 2)
  
  plot(fitted.values(model)[sex == "female"],
       residuals_tmp[sex == "female"], col = "red")
  lines(fitted.values(model)[df$sex == "male"],
        residuals_tmp[sex == "male"], col = "blue", type = "p")
  legend("topright", 95,
         legend=c("female", "male"),
         col=c("red", "blue"), pch=c(1,1))
  
  # QQ plot for the random effects
  qqnorm(random.effects(model)[, 1],
         main = "Q-Q Plot for the random intercept")
  qqline(random.effects(model)[, 1], col = "steelblue", lwd = 2)
}
plot_residuals_lme4 = function(model){
  
  par(mfrow = c(2, 2))
  residuals_tmp = residuals(model)
  
  # plot of the residuals themselves
  plot(residuals_tmp)
  
  # QQ plot for the residuals
  qqnorm(residuals_tmp, main = "Q-Q Plot for the residuals")
  qqline(residuals_tmp, col = "steelblue", lwd = 2)
  
  plot(fitted.values(model)[sex == "female"],
       residuals_tmp[sex == "female"], col = "red")
  lines(fitted.values(model)[df$sex == "male"],
        residuals_tmp[sex == "male"], col = "blue", type = "p")
  legend("topright", 95,
         legend=c("female", "male"),
         col=c("red", "blue"), pch=c(1,1))
  
  # QQ plot for the random effects
  qqnorm(random.effects(model)$subjId[, 1],
         main = "Q-Q Plot for the random intercept")
  qqline(random.effects(model)$subjId[, 1], col = "steelblue", lwd = 2)
}
plot_histograms_gender = function(model){
  par(mfrow = c(1, 1))
  hist(model$residuals[, 2], breaks = 30,
       col = rgb(1, 0, 0, alpha = 0.5))
  hist(model$residuals, breaks = 30, add=T,
       col = rgb(0, 0, 1, alpha = 0.5))
  legend("topright", 95,
         legend=c("Mixed effects", "Linear model"),
         col=c("red", "blue"), pch=c(1,1))
}

## MEF based on LM -------------------------------------------------------------

# Initial full model
# "Singularity in backsolve at level 0, block 1": when sex:tOut interaction
# was present

# rescailing
df[c("tOut_rescaled", "tInOp_rescaled")] = cbind(
  (df$tOut - mean(df$tOut)) / sd(df$tOut),
  (df$tInOp - mean(df$tInOp)) / sd(df$tInOp)
)
attach(df)

### Log transformation ---------------------------------------------------------
mlog_1 = lme(log(clo) ~ poly(tInOp, 2)*poly(tOut, 2)*sex,
            random = ~1|subjId, method = "ML")
summary(mlog_1)
plot_residuals_nlme(mlog_1)

mlog_2 = lmer(log(clo) ~ poly(tInOp, 2)*poly(tOut, 2)*sex + 
            (sex + tOut + tInOp|subjId), REML = F)
summary(mlog_2)

random.effects(mlog_2)
anova(mlog_2)

### Model reduction -------------------------------------------------------------
drop1(mef_1, test = "Chisq")
mef_1 = update(mef_1, . ~ . -  poly(tInOp, 1):poly(tOut, 1))

drop1(mef_1, test = "Chisq")
mef_1 = update(mef_1, . ~ . -  sex)

drop1(mef_1, test = "Chisq")

summary(mef_1)

### sqrt transformation --------------------------------------------------------

## Creation of the Hat matrix --------------------------------------------------

## Derivatives for the transformations -----------------------------------------
sqrt_derivative <- function(y_input){
  return(1/(2*sqrt(y_input)))
}
to_subtract <- sum(log(abs(sqrt_derivative(clo))))
AIC_sqrt_original <- AIC(lm2) - 2*to_subtract
BIC_sqrt_original <- BIC(lm2) - 2*to_subtract

log_derivative <- function(y_input){
  return(1/y_input)
}
to_subtract_log <- sum(log(abs(log_derivative(clo))))
AIC_log_original <- AIC(lm2_log) - 2*to_subtract_log
BIC_log_original <- BIC(lm2_log) - 2*to_subtract_log


# reestimation of the model
mef_1 = lme(formula(mef_1), random = ~1|subjId, method = "REML")
summary(mef_1)

# even though the p-value is large the likelihood might be not symetric and
# the LRT should be prefered, considering the high number of observations
plot_residuals_nlme(mef_1)

lm_test = lm(formula(mef_1))
anova(mef_1, lm_test)
# random effects should be kept, however it's not a fair comparison
# when no random effects are presetnt the variabtion is explained by other 
# variable as it's done with mlog model

## MEF full -------------------------------------------------------------
mef_2 = lme(log(clo) ~ poly(tInOp, 2)*poly(tOut, 2)*sex,
            random = ~1|subjId, method = "ML")
summary(mef_2)
plot_residuals_nlme(mef_2)
# again, heavy tails :(

### Model reduction ------------------------------------------------------------

drop1(mef_2, test = "Chisq")
# I can't drop anything :(

# Single term deletions using Satterthwaite's method:
# accounting for effective degrees of freedom
mef_2_lme4 = lmer(log(clo) ~ poly(tInOp, 2)*poly(tOut, 2)*sex
                 + (1|subjId))
drop1(mef_2_lme4, test = "Chisq")

mef_2_lme4 = update(mef_2_lme4, . ~ . -  poly(tInOp, 2):poly(tOut, 2):sex)
drop1(mef_2_lme4, test = "Chisq")

mef_2_lme4 = update(mef_2_lme4, . ~ . -  poly(tInOp, 2):sex)
drop1(mef_2_lme4, test = "Chisq")

mef_2_lme4 = update(mef_2_lme4, . ~ . -  poly(tInOp, 2):poly(tOut, 2))
drop1(mef_2_lme4, test = "Chisq")
summary(mef_2_lme4)
AIC(mef_2_lme4)

plot_residuals_lme4(mef_2_lme4)

# Nested random effects --------------------------------------------------------
mef_nested_1 = lme(log(clo)~ poly(tInOp, 2)*poly(tOut, 2)*sex,
                   random=~1|subjId/day, method = "REML")
summary(mef_nested_1)
plot_residuals_nlme(mef_nested_1)

sd(residuals(mef_nested_1))
mean(residuals(mef_nested_1))

random.effects(mef_nested_1)

AIC(mef_nested_1)
# -1120.284

random = random.effects(mef_nested_1)
rownames(random$day["(Intercept)"])
mylist = str_split(rownames(random$day["(Intercept)"]), "/")
# rbindlist(mylist, fill=TRUE)

array_tmp = array(NA, dim = c(length(mylist), 2))
for (list_elem in seq(length(mylist))){
  array_tmp[list_elem, ] = as.numeric(unlist(mylist[list_elem]))
}

subjId_length_vector = rep(NA, dim(random$subjId)[1])
for(list_elem in seq(dim(random$subjId)[1])){
  subjId_length = sum(random_nested_df$subjId == (list_elem - 1))
  subjId_length_vector[list_elem] = subjId_length
}

list_tmp = c()
values_random = as.numeric(unlist(random$subjId))
for(list_elem in seq(length(values_random))){
  list_tmp = c(list_tmp, rep(values_random[list_elem],
                             each = subjId_length_vector[list_elem]))
}

random_nested_df = data.frame(array_tmp, random$day["(Intercept)"], list_tmp)
colnames(random_nested_df) = c("subjId", "day", "nested_effect", "random_intercept")
random_nested_df[c("subjId", "day")] <- lapply(random_nested_df[c("subjId", "day")], as.factor)

par(mfrow = c(1, 1))
points(random_nested_df$subjId, random_nested_df$nested_effect, type = "p")

## Approximation with SubDay ---------------------------------------------------
mef_nested_approximation = lme(log(clo)~ poly(tInOp, 2)*poly(tOut, 2)*sex,
                   random=~1|subjId/subDay, method = "REML")
summary(mef_nested_approximation)
plot_residuals_nlme(mef_nested_approximation)
AIC(mef_nested_approximation)
# -1092.927

# Repeated measurments set up --------------------------------------------------
mef_gau = lme(log(clo)~ poly(tInOp, 2)*poly(tOut, 2)*sex,
             random=~1|subDay, method = "ML",
             correlation = corGaus(form = ~time2|subDay, nugget = F))
summary(mef_gau)

mef_ar1 = lme(log(clo)~ poly(tInOp, 2)*poly(tOut, 2)*sex,
              random=~1|subDay, method = "ML",
              correlation = corAR1(form = ~time2|subDay),
              data = df)
summary(mef_ar1)

# note time vs time2
mef_exp = lme(log(clo)~ poly(tInOp, 2)*poly(tOut, 2)*sex,
              random=~1|subDay, method = "ML",
              correlation = corExp(form = ~time2|subDay, nugget = F))
summary(mef_exp)
plot_residuals_nlme(mef_exp)

anova(mef_gau, mef_ar1, mef_exp)

par(mfrow = c(1, 2))
plot(Variogram(mef_gau), main = "Gaussian")
plot(Variogram(mef_ar1), main = "AR1")
plot(Variogram(mef_exp), main = "Exp")

anova()

# conversion to canocial parameters - LECTURE 9
nu.sq<-0.1404056^2
sigma.sq<-0.2171559^2*0.2186743
tau.sq<-0.2171559^2-sigma.sq
rho.sq<-2.3863954
c(nu.sq=nu.sq, sigma.sq=sigma.sq, tau.sq=tau.sq, rho.sq=rho.sq)

# TEST FOR THE NUGGET EFFECT AS IN LECTURE 9 R CODE


# Graphical presentation -------------------------------------------------------
install.packages("data.table")

fm1 <- lme(distance ~ age, Orthodont, random = ~ age | Subject)
## normal plot of standardized residuals by gender
qqnorm(fm1, ~ resid(fm1, type = "p")|Sex, abline = c(0, 1))
## normal plots of random effects
qqnorm(fm1, ~ranef(.))

summary(fm1)

ranef(fm1)

# FROM LECTURE 8 - NICE PLOTTING
# Group over subjects
plot(mef_2_lme4,form=resid(.,type='p')~fitted(.)|subjId,abline=0)
# Goup over Type of chair
plot(fit1,form=resid(.,type='p')~fitted(.)|Type,abline=0)


