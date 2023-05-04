# File Description -------------------------------------------------------------
#
#   Tymoteusz Barcinski - s221937
#   
#   Advanced Dataanalysis and Statistical Modelling
#   Assignment 2 - Part A
#
#_______________________________________________________________________________
rm(list=ls()) 
print(utils::getSrcDirectory(function(){}))
print(utils::getSrcFilename(function(){}, full.names = TRUE))
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()
Sys.setenv(LANG = "en")

library(MASS)
library(dplyr)
library(ggplot2)
library(car)
library(patchwork)
library(mvtnorm)
library(nlme)
library(numDeriv)
library(lmtest)
library(GGally) 
library(VGAM)

# Reading data
dat <- read.table("data/clothing.csv", sep=',', head=TRUE)
df <- data.frame(dat) %>% select(-X) 
cols <- c("sex", "subjId", "day")
df[cols] <- lapply(df[cols], as.factor)
summary(df)
sum(is.na(df)) # no missing values
attach(df)

# EDA ###################
dim(df)
table(sex)
table(day)
length(unique(subjId))

#investigating time
plot(time, type = "p")
df %>% group_by(subjId, day) %>% summarize(n = n()) %>% 
  print(n = 136)
df %>% group_by(subjId) %>% summarize(n = n()) %>%
  group_by(n) %>% summarize(sum = n())

p_continous <- ggpairs(df, aes(color = sex),
                 columns = c("tOut", "tInOp", "clo"))
p_continous

table(day, sex)
length(unique(clo))
attach(df)

# Modelling #
library("car")
qqPlot(residuals.glm(m1, type = "deviance"))

plot_residuals <- function(glm_object){
  
  deviance = residuals.glm(glm_object, type = "deviance")
  pearson = residuals.glm(glm_object, type = "pearson")
  fitted_values = 
  
  par(mfrow = c(3, 2))
  # plot(cumsum(df$time), deviance, type = "p", col=df$sex)
  # plot(cumsum(df$time), pearson, type = "p", col=df$sex)
  
  # plot(deviance, type = "p", col=df$sex)
  # plot(pearson, type = "p", col=df$sex)
  
  # hist(deviance, col=df$sex)
  # hist(pearson, col=df$sex)
  
  qqnorm(deviance)
  qqline(deviance, col = "steelblue", lwd = 2)
  qqnorm(pearson)
  qqline(pearson, col = "steelblue", lwd = 2)
  
  # acf_deviance = acf(deviance, plot=F)
  # acf_pearson = acf(pearson, plot=F)
  # plot(acf_deviance)
  # plot(acf_pearson)
  
  plot(fitted.values(glm_object)[df$sex == "female"],
       deviance[df$sex == "female"], col = "red")
  lines(fitted.values(glm_object)[df$sex == "male"],
        deviance[df$sex == "male"], col = "blue", type = "p")
  legend("topright", 95,
         legend=c("female", "male"),
         col=c("red", "blue"), pch=c(1,1))
  
  plot(fitted.values(glm_object)[df$sex == "female"],
       pearson[df$sex == "female"], col = "red")
  lines(fitted.values(glm_object)[df$sex == "male"],
        pearson[df$sex == "male"], col = "blue", type = "p")
  legend("topright", 95,
         legend=c("female", "male"),
         col=c("red", "blue"), pch=c(1,1))
  
  # breaks_number = 50
  upper_limit = 150
  
  hist(deviance[df$sex == "male"], col = rgb(0, 0, 1, alpha = 0.5),
       ylim = c(0, upper_limit), xlim = range(deviance))
  hist(deviance[df$sex == "female"], add=T, col = rgb(1, 0, 0, alpha = 0.5),
       ylim = c(0, upper_limit), xlim = range(deviance))
  legend("topright", 95,
         legend=c("female", "male"),
         col=c("red", "blue"), pch=c(1,1))
  
  hist(pearson[df$sex == "male"], col = rgb(0, 0, 1, alpha = 0.5),
       ylim = c(0, upper_limit), xlim = range(pearson))
  hist(pearson[df$sex == "female"], add=T, col = rgb(1, 0, 0, alpha = 0.5),
       ylim = c(0, upper_limit), xlim = range(pearson))
  legend("topright", 95,
         legend=c("female", "male"),
         col=c("red", "blue"), pch=c(1,1))
  
  par(mfrow = c(1, 1))
}

# Gamma models ###################################
m1 = glm(clo ~ tInOp + tOut + sex, family = Gamma(link = "inverse"))
summary(m1)
plot_residuals(m1)
par(mfrow = c(2, 2))
plot(m1)


m3 = glm(clo ~ sex*poly(tInOp, 3) + sex*poly(tOut, 3) + sex +
           poly(tInOp, 2)*poly(tOut, 2),
         family = Gamma(link = "log"))
summary(m3)
# plot_residuals(m3)

## log link #########
### 1
drop1(m3, test = "Chisq")
m3_reduction = glm(clo ~ sex*poly(tInOp, 2) + sex*poly(tOut, 3) + sex +
           poly(tInOp, 2)*poly(tOut, 2),
         family = Gamma(link = "log"))
anova(m3, m3_reduction, test = "LRT")
m3 <- m3_reduction

### 2
drop1(m3, test = "Chisq")
m3_reduction = glm(clo ~ sex*poly(tInOp, 1) + sex*poly(tOut, 3) + sex +
                     poly(tInOp, 2)*poly(tOut, 2),
                   family = Gamma(link = "log"))
anova(m3, m3_reduction, test = "LRT")
m3 <- m3_reduction

### 3
drop1(m3, test = "Chisq")
m3_reduction = glm(clo ~ tInOp + sex*poly(tOut, 3) + sex +
                     poly(tInOp, 2)*poly(tOut, 2),
                   family = Gamma(link = "log"))
anova(m3, m3_reduction, test = "LRT")
m3 <- m3_reduction

### 4
drop1(m3, test = "Chisq")
m3_reduction = glm(clo ~ tInOp + sex*poly(tOut, 3) + sex +
                     poly(tInOp, 2)*poly(tOut, 1) + poly(tInOp, 1):poly(tOut, 2),
                   family = Gamma(link = "log"))
anova(m3, m3_reduction, test = "LRT")
m3 <- m3_reduction

### 5
drop1(m3, test = "Chisq")
m3_reduction = glm(clo ~ tInOp + sex*poly(tOut, 3) + sex +
                     poly(tInOp, 2)*poly(tOut, 1),
                   family = Gamma(link = "log"))
anova(m3, m3_reduction, test = "LRT")
m3 <- m3_reduction

### 6
drop1(m3, test = "Chisq")
m3_reduction = glm(clo ~ tInOp + sex*poly(tOut, 2) + sex +
                     poly(tInOp, 2)*poly(tOut, 1),
                   family = Gamma(link = "log"))
anova(m3, m3_reduction, test = "LRT")
m3 <- m3_reduction

### 7
drop1(m3, test = "Chisq")
m3_reduction = glm(clo ~ tInOp + sex:poly(tOut, 1) + I(tOut^2) + sex +
                     poly(tInOp, 2)*poly(tOut, 1),
                   family = Gamma(link = "log"))
anova(m3, m3_reduction, test = "LRT")
m3 <- m3_reduction
### no further reductions possible
summary(m3)

### 8
drop1(m3, test = "Chisq")
m3_reduction = glm(clo ~ tInOp + tOut + sex:poly(tOut, 1) + I(tOut^2) + sex +
                     poly(tInOp, 1):poly(tOut, 1),
                   family = Gamma(link = "log"))
anova(m3, m3_reduction, test = "LRT")
m3 <- m3_reduction
summary(m3)

### 9
drop1(m3, test = "Chisq")
m3_reduction = glm(clo ~ sex*poly(tInOp, 2) + sex*poly(tOut, 3) + sex +
                     poly(tInOp, 2)*poly(tOut, 1),
                   family = Gamma(link = "log"))
anova(m3, m3_reduction, test = "LRT")
m3 <- m3_reduction

m4 <- glm(clo ~ poly(tInOp, 2) + poly(tOut, 2) + sex + sex:tOut+
            tInOp:tOut, family = Gamma(link = "log"))
summary(m4)
par(mfrow = c(2, 2))
plot(m4)
## heavy tails unfortunately :(
anova(m4, m1, test = "Chisq")

############## inverse link ####################
m3_inverse = glm(clo ~ sex*poly(tInOp, 3) + sex*poly(tOut, 3) + sex +
           poly(tInOp, 2)*poly(tOut, 2),
         family = Gamma(link = "inverse"))

### 1
drop1(m3_inverse, test = "Chisq")
m3_inverse_reduction = glm(clo ~ sex*poly(tInOp, 2) + sex*poly(tOut, 3) +
                             sex +
                     poly(tInOp, 2)*poly(tOut, 2),
                   family = Gamma(link = "inverse"))
anova(m3_inverse, m3_inverse_reduction, test = "LRT")
m3_inverse <- m3_inverse_reduction

### 2
drop1(m3_inverse, test = "Chisq")
m3_inverse_reduction = glm(clo ~ sex:poly(tInOp, 1) + poly(tInOp, 2) +
                             sex*poly(tOut, 3) + sex +
                             poly(tInOp, 2)*poly(tOut, 2),
                           family = Gamma(link = "inverse"))
anova(m3_inverse, m3_inverse_reduction, test = "LRT")
m3_inverse <- m3_inverse_reduction

### 3
drop1(m3_inverse, test = "Chisq")
m3_inverse_reduction = glm(clo ~ poly(tInOp, 2) + sex*poly(tOut, 3) + sex +
                             poly(tInOp, 2)*poly(tOut, 2),
                           family = Gamma(link = "inverse"))
anova(m3_inverse, m3_inverse_reduction, test = "LRT")
m3_inverse <- m3_inverse_reduction

### 4
drop1(m3_inverse, test = "Chisq")
m3_inverse_reduction = glm(clo ~ poly(tInOp, 2) + sex*poly(tOut, 3) + sex +
                             poly(tInOp, 1):poly(tOut, 1) + 
                             poly(tInOp, 2):poly(tOut, 1) +
                             poly(tInOp, 1):poly(tOut, 2),
                           family = Gamma(link = "inverse"))
anova(m3_inverse, m3_inverse_reduction, test = "LRT")
m3_inverse <- m3_inverse_reduction

### 5
drop1(m3_inverse, test = "Chisq")
m3_inverse_reduction = glm(clo ~ poly(tInOp, 2) + sex*poly(tOut, 3) + sex +
                             poly(tInOp, 1):poly(tOut, 1),
                           family = Gamma(link = "inverse"))
anova(m3_inverse, m3_inverse_reduction, test = "LRT")
m3_inverse <- m3_inverse_reduction

### 6
drop1(m3_inverse, test = "Chisq")
m3_inverse_reduction = glm(clo ~ poly(tInOp, 2) +
                             sex:poly(tOut, 1) + sex:poly(tOut, 2) + 
                             poly(tOut,3) + sex +
                             poly(tInOp, 1):poly(tOut, 1),
                           family = Gamma(link = "inverse"))
anova(m3_inverse, m3_inverse_reduction, test = "LRT")
m3_inverse <- m3_inverse_reduction

### 7
drop1(m3_inverse, test = "Chisq")
m3_inverse_reduction = glm(clo ~ poly(tInOp, 2) +
                             sex:poly(tOut, 1) + sex:poly(tOut, 2) + 
                             poly(tOut,2) + sex +
                             poly(tInOp, 1):poly(tOut, 1),
                           family = Gamma(link = "inverse"))
anova(m3_inverse, m3_inverse_reduction, test = "LRT")
m3_inverse <- m3_inverse_reduction

### 8
drop1(m3_inverse, test = "Chisq")
m3_inverse_reduction = glm(clo ~ poly(tInOp, 2) +
                             sex:poly(tOut, 1) + 
                             poly(tOut,2) + sex +
                             poly(tInOp, 1):poly(tOut, 1),
                           family = Gamma(link = "inverse"))
anova(m3_inverse, m3_inverse_reduction, test = "LRT")
m3_inverse <- m3_inverse_reduction

### 9
drop1(m3_inverse, test = "Chisq")

m4_inverse <- glm(clo ~ poly(tInOp, 2) + poly(tOut, 2) + sex + sex:tOut+
                    tInOp:tOut, family = Gamma(link = "inverse"))
summary(m4_inverse)
par(mfrow = c(2, 2))
plot(m4)

c(AIC(m4), AIC(m4_inverse))
c(BIC(m4), BIC(m4_inverse))

# model_lm_final <- model_transformed_months_only_sin
# to_subtract <- sum(log(abs(transformation_1_derivative(optimal_lambda, df$pow_obs))))
# (AIC_transformed_domain <- AIC(model_lm_final) - 2*to_subtract)




###### 4. Including subject ID
weights_vector[df$sex == "female"] <- 0.4
m8 = glm(clo ~ tInOp*tOut*subjId*sex, family = Gamma(link = "inverse"),
         weights = weights_vector)
AIC(m8)
summary(m8)
par(mfrow = c(2, 2))
plot(m8)
plot_residuals(m8)
shapiro.test(m8$residuals)

# weights_vector[df$sex == "female"] <- 0.7
m9 = glm(clo ~ poly(tInOp, 2) + poly(tOut, 2) +
           tInOp:tOut + subjId + subjId:tInOp + subjId:tOut,
         family = Gamma(link = "inverse"))
summary(m9)
par(mfrow = c(2, 2))
plot(m9)
plot_residuals(m9)
shapiro.test(m9$residuals)

drop1(m9, test = "Chisq")
m9 <- update


# Linear model #########################
lm1 <- lm(sqrt(clo) ~ poly(tInOp, 1) + poly(tOut, 1) + sex)
summary(lm1)
par(mfrow = c(2, 2))
plot(lm1)
shapiro.test(lm1$residuals)

## sqrt transformation ##########
lm2 <- lm(sqrt(clo) ~ sex*poly(tInOp, 2) + sex*poly(tOut, 2) + sex +
            poly(tInOp, 2)*poly(tOut, 2))
summary(lm2)
# par(mfrow = c(2, 2))
# plot(lm2)

### 1
drop1(lm2, test = "F")
lm2_reduction <- lm(sqrt(clo) ~ sex*poly(tInOp, 2) + sex*poly(tOut, 2) + sex +
            poly(tInOp, 1):poly(tOut, 1) +
            poly(tInOp, 2):poly(tOut, 1) +
            poly(tInOp, 1):poly(tOut, 2))
anova(lm2, lm2_reduction, test = "LRT")
lm2 <- lm2_reduction

### 2
drop1(lm2, test = "F")
lm2_reduction <- lm(sqrt(clo) ~ sex*poly(tInOp, 2) + sex*poly(tOut, 2) + sex +
                      poly(tInOp, 1):poly(tOut, 1))
anova(lm2, lm2_reduction, test = "LRT")
lm2 <- lm2_reduction

### 3
drop1(lm2, test = "F")
lm2_reduction <- lm(sqrt(clo) ~ sex*poly(tInOp, 2) + sex:poly(tOut, 1) +
                      sex + poly(tOut, 2) +
                      poly(tInOp, 1):poly(tOut, 1))
anova(lm2, lm2_reduction, test = "LRT")
lm2 <- lm2_reduction

### 4
drop1(lm2, test = "F")
lm2_reduction <- lm(sqrt(clo) ~ sex:poly(tInOp, 1) + sex:poly(tOut, 1) +
                      sex + poly(tOut, 2) + poly(tInOp, 2) +
                      poly(tInOp, 1):poly(tOut, 1))
anova(lm2, lm2_reduction, test = "LRT")
lm2 <- lm2_reduction

### 5
drop1(lm2, test = "F")
lm2_reduction <- lm(sqrt(clo) ~  + sex:poly(tOut, 1) +
                      sex + poly(tOut, 2) + poly(tInOp, 2) +
                      poly(tInOp, 1):poly(tOut, 1))
anova(lm2, lm2_reduction, test = "LRT")
lm2 <- lm2_reduction

### 6
drop1(lm2, test = "F")
lm2_reduction <- lm(sqrt(clo) ~  + sex:poly(tOut, 1) +
                      sex + poly(tOut, 2) + poly(tInOp, 1) +
                      poly(tInOp, 1):poly(tOut, 1))
anova(lm2, lm2_reduction, test = "LRT")
# no further reduction is posiible
summary(lm2)
AIC(lm2) # transformed domain

sqrt_derivative <- function(y_input){
  return(1/(2*sqrt(y_input)))
}
to_subtract <- sum(log(abs(sqrt_derivative(clo))))
AIC_sqrt_original <- AIC(lm2) - 2*to_subtract
BIC_sqrt_original <- BIC(lm2) - 2*to_subtract

lm_full <- lm(sqrt(clo) ~ sex*poly(tInOp, 2)*poly(tOut, 2))
summary(lm_full)

### 4
drop1(lm_full, test = "F")
lm_full_reduction <- lm(sqrt(clo) ~ sex:poly(tInOp, 1) +  sex:poly(tInOp, 2) +
                          sex:poly(tOut, 1) + sex:poly(tOut, 2) +
                      sex + poly(tOut, 2) + poly(tInOp, 2) +
                      poly(tInOp, 1):poly(tOut, 1) +
                      sex:poly(tInOp, 1):poly(tOut, 1))
anova(lm_full, lm_full_reduction, test = "LRT")
lm_full <- lm_full_reduction

### 4
drop1(lm_full, test = "F")
lm_full_reduction <- lm(sqrt(clo) ~ sex:poly(tInOp, 1) +  sex:poly(tInOp, 2) +
                          sex:poly(tOut, 1) +
                          sex + poly(tOut, 2) + poly(tInOp, 2) +
                          poly(tInOp, 1):poly(tOut, 1) +
                          sex:poly(tInOp, 1):poly(tOut, 1))
anova(lm_full, lm_full_reduction, test = "LRT")
lm_full <- lm_full_reduction

### 4
drop1(lm_full, test = "F")
# no further reductions
summary(lm_full)

## log transformation ##########
lm2_log <- lm(log(clo) ~ sex*poly(tInOp, 2) + sex*poly(tOut, 2) + sex +
                poly(tInOp, 2)*poly(tOut, 2))
summary(lm2_log)
# par(mfrow = c(2, 2))
# plot(lm2)

### 1
drop1(lm2_log, test = "F")
lm2_log_reduction <- lm(log(clo) ~ sex*poly(tInOp, 2) + sex*poly(tOut, 2) +
                          sex +
                          poly(tInOp, 1):poly(tOut, 1) +
                          poly(tInOp, 2):poly(tOut, 1) +
                          poly(tInOp, 1):poly(tOut, 2))
anova(lm2_log, lm2_log_reduction, test = "LRT")
lm2_log <- lm2_log_reduction

### 2
drop1(lm2_log, test = "F")
lm2_log_reduction <- lm(log(clo) ~ sex*poly(tInOp, 2) + sex*poly(tOut, 2) + sex +
                          poly(tInOp, 1):poly(tOut, 1) +
                          poly(tInOp, 2):poly(tOut, 1))
anova(lm2_log, lm2_log_reduction, test = "LRT")
lm2_log <- lm2_log_reduction

### 3
drop1(lm2_log, test = "F")
lm2_log_reduction <- lm(log(clo) ~ poly(tInOp, 2) + sex:poly(tInOp, 1) +
                          sex*poly(tOut, 2) + sex +
                          poly(tInOp, 1):poly(tOut, 1) +
                          poly(tInOp, 2):poly(tOut, 1))
anova(lm2_log, lm2_log_reduction, test = "LRT")
lm2_log <- lm2_log_reduction

### 4
drop1(lm2_log, test = "F")
lm2_log_reduction <- lm(log(clo) ~ poly(tInOp, 2) + sex:poly(tInOp, 1) +
                          sex:poly(tOut, 1) + poly(tOut, 2) +sex +
                          poly(tInOp, 1):poly(tOut, 1) +
                          poly(tInOp, 2):poly(tOut, 1))
anova(lm2_log, lm2_log_reduction, test = "LRT")
lm2_log <- lm2_log_reduction

### 5
drop1(lm2_log, test = "F")
lm2_log_reduction <- lm(log(clo) ~ poly(tInOp, 2) +
                          sex:poly(tOut, 1) + poly(tOut, 2) +sex +
                          poly(tInOp, 1):poly(tOut, 1) +
                          poly(tInOp, 2):poly(tOut, 1))
anova(lm2_log, lm2_log_reduction, test = "LRT")
lm2_log <- lm2_log_reduction

### 6
drop1(lm2_log, test = "F")
lm2_log_reduction <- lm(log(clo) ~ poly(tInOp, 2) +
                          sex:poly(tOut, 1) + poly(tOut, 2) +sex +
                          poly(tInOp, 1):poly(tOut, 1))
anova(lm2_log, lm2_log_reduction, test = "LRT")
lm2_log <- lm2_log_reduction

### 6
drop1(lm2_log, test = "F")
lm2_log_reduction <- lm(log(clo) ~ poly(tInOp, 1) +
                          sex:poly(tOut, 1) + poly(tOut, 2) +sex +
                          poly(tInOp, 1):poly(tOut, 1))
anova(lm2_log, lm2_log_reduction, test = "LRT")
# no further reductions are possible

summary(lm2_log)
AIC(lm2_log) # transformed domain

log_derivative <- function(y_input){
  return(1/y_input)
}
to_subtract_log <- sum(log(abs(log_derivative(clo))))
AIC_log_original <- AIC(lm2_log) - 2*to_subtract_log
BIC_log_original <- BIC(lm2_log) - 2*to_subtract_log

### Results
# results_models <- cbind(
#   c(AIC(m4), AIC(m4_inverse), AIC_sqrt_original, AIC_log_original),
#   c(BIC(m4), BIC(m4_inverse), BIC_sqrt_original, BIC_log_original)
# )

results_models <- cbind(
  c(AIC_sqrt_original, AIC_log_original),
  c(BIC_sqrt_original, BIC_log_original)
)
rownames(results_models) <- c("lm_sqrt", "lm_log")
colnames(results_models) <- c("AIC", "BIC")

# Gamma full model - 3 way interactions - estimation ###############################
gamma_full <- glm(clo ~ sex*poly(tInOp, 2)*poly(tOut, 2),
                  family = Gamma(link = "inverse"))
summary(gamma_full)
par(mfrow = c(2, 2))
plot(gamma_full)
plot_residuals(gamma_full)

### 4
drop1(gamma_full, test = "F")
gamma_full_reduction <- glm(clo ~ sex:poly(tInOp, 1) +  sex:poly(tInOp, 2) +
                              sex:poly(tOut, 1) + sex:poly(tOut, 2) +
                              sex + poly(tOut, 2) + poly(tInOp, 2) +
                              poly(tInOp, 1):poly(tOut, 1) +
                              sex:poly(tInOp, 1):poly(tOut, 1),
                            family = Gamma(link = "inverse"))
anova(gamma_full, gamma_full_reduction, test = "LRT")
gamma_full <- gamma_full_reduction

### 4
drop1(gamma_full, test = "F")
gamma_full_reduction <- glm(clo ~ sex:poly(tInOp, 1) +  sex:poly(tInOp, 2) +
                              sex:poly(tOut, 1) + 
                              sex + poly(tOut, 2) + poly(tInOp, 2) +
                              poly(tInOp, 1):poly(tOut, 1) +
                              sex:poly(tInOp, 1):poly(tOut, 1),
                            family = Gamma(link = "inverse"))
anova(gamma_full, gamma_full_reduction, test = "LRT")
# no further reduction
summary(gamma_full)

gamma_without_na <- glm(clo ~ tInOp + tOut + sex + sex:tOut,
                        family = Gamma(link = "inverse"))

gamma_testing <- glm(clo ~ tInOp + tOut + sex + sex:tOut,
                  family = Gamma(link = "inverse"))
summary(gamma_testing)

c(AIC(gamma_full), AIC_sqrt_full, BIC(gamma_full), BIC_sqrt_full)

AIC_sqrt_full <- AIC(lm_full) - 2*to_subtract
BIC_sqrt_full <- BIC(lm_full) - 2*to_subtract

par(mfrow = c(2, 2))
plot(gamma_full)


# breaks_number = 40
# par(mfrow = c(1, 1))
# hist(lm2$residuals[df$sex == "female"], col = "red", breaks = breaks_number,
#      ylim = c(0, upper_limit))
# hist(lm2$residuals[df$sex != "female"], add=T, breaks = breaks_number,
#      ylim = c(0, upper_limit), col = rgb(0, 0, 0.5, alpha = 0.1))


# we should pick Gamma_inverse model both by AIC, BIC

summary(m4_inverse)
par(mfrow = c(2, 2))
plot(m4_inverse)

m4_subject <- glm(clo ~ poly(tInOp, 2) + poly(tOut, 2) + subjId + subjId:tOut +
                    tInOp:tOut, family = Gamma(link = "inverse"))
summary(m4_subject)
drop1(m4_subject, test = "LRT")
AIC(m4_subject) # much better AIC


# Gamma dispersion parameter #############

m1 = glm(clo ~ tInOp + tOut, family = Gamma(link = "inverse"),
         data = df)
logLik(m1)

### Weighted version
weights_vector <- rep(1, dim(df)[1])
# weights_vector[df$sex == "female"] <- exp(-0.485523736)
weights_vector[df$sex == "female"] <- exp(0.099)
m1 = glm(clo ~ tInOp + tOut, family = Gamma(link = "inverse"),
         data = df, weights = weights_vector)
summary(m1)
# plot_residuals(m1)
par(mfrow = c(2, 2))
plot(m1)

logLik(m1)
# 'log Lik.' 448.0734 (df=5)

X = model.matrix(m1)
X_dispersion = cbind(rep(1, dim(X)[1]), as.integer(df$sex == "female"))
# predictions_scratch = X %*% as.numeric(m1$coefficients)
# norm(1/m1$fitted.values - predictions_scratch)

weights_vector <- rep(1, dim(df)[1])
theta_initial = c(rep(1, length(m1$coefficients)), "dispersion" = -1, "female_weight" = -0.48)

# theta_initial = c(rep(1, length(m1$coefficients)), 0.05)
### own optimization
gamma_own <- function(theta){
  print(as.numeric(theta[4:5]))
  
  theta_linear_space = X %*% theta[1:dim(X)[2]]
  # inverse transformation
  # mean = 1/theta_linear_space
  # log transformation
  
  # regressions
  mean = exp(theta_linear_space)
  dispersion = exp(X_dispersion %*% theta[4:5])
  
  # weights_vector[df$sex == "female"] <- exp(theta[length(theta)])
  # modelling the dispersion parameter rather than precision
  # scale_parameter = rep(mean, length(shape_parameter))*shape_parameter
  
  weights_vector[df$sex == "female"] <- exp(theta[length(theta)])
  shape_parameter = 1/dispersion
  scale_parameter = mean*dispersion

  # shape_parameter = 1/theta[length(theta) - 1]
  # scale_parameter = mean*theta[length(theta) - 1]

  return(-sum(dgamma(df$clo, shape = shape_parameter,
                     scale = scale_parameter, log = TRUE)))
}
gamma_own(theta_initial)
opt <- nlminb(theta_initial, gamma_own)
# -3.36093943    1.06208529
# 0.03470264    2.89239620

opt$par
m1$coefficients
gamma_own(as.numeric(opt$par))

### with the package

library(Gammareg)
X2 <- df$tInOp
X3 <- df$tOut
Z1 <- as.integer(df$sex == "female")
Z2 <- as.integer(df$sex == "male")
formula.mean = df$clo ~ X2 + X3
formula.shape = ~ Z1
a=Gammareg(formula.mean, formula.shape, meanlink="log")
summary(a)
opt$par
exp(opt$par[4:5])d:

### Optimal stuff: 1/16.369001352  -0.485523736 

library(Gammareg)

X1 <- rep(1, 500)
X2 <- runif(500, 0, 30)
X3 <- runif(500, 0, 15)
X4 <- runif(500, 10, 20)
mui <- 15 + 2*X2 + 3*X3
alphai <- exp(0.2 + 0.1*X2 + 0.3*X4)
Y <- rgamma(500, shape=alphai, scale=mui/alphai)
X <- cbind(X1, X2, X3)
Z <- cbind(X1, X2, X4)
formula.mean = Y ~ X2 + X3
formula.shape = ~ X2 + X4
a=Gammareg(formula.mean, formula.shape, meanlink="ide")
summary(a)


## Linear model working !!!! ######################
objective_weights_lm = function(weight_female_function){
  weights_vector[df$sex == "female"] <- exp(weight_female_function)
  m1_tmp = lm(1/(clo) ~ tInOp + tOut, data = df, weights = weights_vector)
  return(-as.numeric(logLik(m1_tmp)))
}
objective_weights_lm(1)
opt <- nlminb(-0.5, objective_weights_lm)

# NOT WORKING
# model_optimization <- m1
weights_vector <- rep(1, dim(df)[1])
objective_weights_not_working = function(theta){
  print(theta)
  weights_vector[df$sex == "female"] <-exp(theta)
  # weights_vector <- theta[2]*weights_vector
  m1_tmp = glm(clo ~ tInOp + tOut, family = Gamma(link = "inverse"),
               data = df, weights = weights_vector)
  return(-as.numeric(logLik(m1_tmp)))
}

# theta_initial = c(1, "dispersion" = 0.05)
theta_initial = -0.5
objective_weights_not_working(theta_initial)

opt <- nlminb(theta_initial, objective_weights)
exp(opt$par)
opt

# Profile Likelihood ####################
model_optimization <- gamma_full
X = model.matrix(model_optimization)
n = dim(X)[1]
p = dim(X)[2]
y = df$clo
X_dispersion = cbind(rep(1, n), as.integer(df$sex == "female"))

# weights_vector <- rep(1, dim(df)[1])
objective <- function(theta){
  theta_linear_space = X %*% theta[1:p]
  # regressions: log link for mu and log link for dispersion
  mean = exp(theta_linear_space)
  weights_vector[df$sex == "female"] <- theta[p+2]
  # dispersion = exp(as.matrix(X_dispersion[, 1]) %*% theta[p+1])*weights_vector
  dispersion = exp(X_dispersion %*% theta[(p+1):(p+2)])
  shape_parameter = 1/dispersion
  scale_parameter = mean*dispersion
  return(-sum(dgamma(y, shape = shape_parameter,
                     scale = scale_parameter, log = TRUE)))
}
theta_initial = c(rep(1, length(model_optimization$coefficients)), "dispersion_male" = -3,
                  "female_multiplicative" = 1)
objective(theta_initial)
opt <- nlminb(theta_initial, objective)
opt
exp(opt$par[4:5])

residuals <- y - X %*% opt$par[1:p]

breaks_number = 70
upper_limit = 20
par(mfrow = c(1, 1))
hist(residuals[df$sex == "female"], col = "red", breaks = breaks_number,
     ylim = c(0, upper_limit))
hist(residuals[df$sex != "female"], add=T, breaks = breaks_number,
     ylim = c(0, upper_limit), col = rgb(0, 0, 0.5, alpha = 0.1))
shapiro.test(residuals)

## wald uncertentiy #####
hessian_matrix = hessian(objective, opt$par)
standard_error = sqrt(diag(solve(hessian_matrix)))

Wald_lower_CI = opt$par["female_multiplicative"] -
  qnorm(0.975)*standard_error[length(standard_error)]
Wald_upper_CI = opt$par["female_multiplicative"] +
  qnorm(0.975)*standard_error[length(standard_error)]
c(Wald_lower_CI, Wald_upper_CI)

profile_objective <- function(weight){
  fun.tmp <- function(theta_inner, weight_input){
    objective(c(theta_inner, weight_input))
  }
  theat_initial_inner_opt = c(model_optimization$coefficients,
                              "dispersion_male" = -3)
  nlminb(theat_initial_inner_opt, fun.tmp, weight_input = weight)$objective
}
(profile_objective1)

p1 <- seq(Wald_lower_CI - 0.5, Wald_upper_CI + 0.5, by = 0.1)
logLp1 <- sapply(p1, profile_objective) ## note sapply!
logLp1 <- logLp1 - min(logLp1) ## normalization
L_CI_lower = min(p1[exp(-logLp1) > exp(-qchisq(0.95,df=1)/2)])
L_CI_upper =  max(p1[exp(-logLp1) > exp(-qchisq(0.95,df=1)/2)])

obs_hessian <- hessian(profile_objective, opt$par["female_multiplicative"])[1, 1]
quadratic_approximation <-
  exp( -0.5* obs_hessian * (p1 - opt$par["female_multiplicative"])^2)
quadratic_approximation <- quadratic_approximation / max(quadratic_approximation)

par(mfrow = c(1,1))
plot(p1, exp(-logLp1), type = "l",
     xlab="Ratio", ylab="Profile Likelihood",
     main="The comparison between Wald's and Likelihood based Confidence Intervals")
# axis(side=1, at=seq(0, 10, by=1))
lines(p1, rep(exp(-qchisq(0.95,df=1)/2), length(p1)), col = 2)
rug(L_CI_lower, ticksize = 0.1, lwd = 2, col = "red")
rug(L_CI_upper, ticksize = 0.1, lwd = 2, col = "red")
lines(p1, quadratic_approximation, col = "blue")
abline(v = opt$par["female_multiplicative"], lty = 2)
rug(Wald_lower_CI, ticksize = 0.1, lwd = 2, col = "blue")
rug(Wald_upper_CI, ticksize = 0.1, lwd = 2, col = "blue")
legend("topright", 95,
       legend=c("Profile likelihood", "Quadratic approximation",
                "95% confidence interval"),
       col=c("black", "blue", "red"), lty = 1:1, cex=0.8,
       inset = 0.02)

c(L_CI_lower, L_CI_upper)

wald_statistic = (opt$par["weight"] - 1)/standard_error[length(standard_error)]

# LRT
ll_full <- -profile_objective(as.numeric(opt$par["female_multiplicative"]))
ll_test <- -profile_objective(1)
LRT <- -2*(ll_test - (ll_full))
p_value <- 1 - pchisq(LRT, df = 1)
p_value


