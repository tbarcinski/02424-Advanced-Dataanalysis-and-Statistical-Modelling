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
dat <- read.table("clothing.csv", sep=',', head=TRUE)
df <- data.frame(dat) %>% select(-X) 
cols <- c("sex", "subjId", "day")
df[cols] <- lapply(df[cols], as.factor)
summary(df)
sum(is.na(df)) # no missing values
attach(df)

########### EDA ##################
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

########## Modelling #################################
boxcox_object = boxcox(df$clo ~ 1)
lambda = boxcox_object$x[which.max(boxcox_object$y)]
df["clo_trans"] = (df$clo^lambda - 1)/lambda
hist(df$clo_trans, breaks = 50)
attach(df)

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

#### Simplest model 
m1 = glm(clo ~ tInOp + tOut + sex, family = Gamma(link = "inverse"))
summary(m1)
plot_residuals(m1)
par(mfrow = c(2, 2))
plot(m1)

# m2 = glm(clo ~ tInOp*tOut*sex, family = Gamma(link = "log"))
# summary(m2)
# par(mfrow = c(2, 2))
# plot(m2)
# 
# drop1(m2, test = "Chisq")
# m2 <- update(m2, . ~ . - tInOp:tOut:sex)
# drop1(m2, test = "Chisq")
# m2 <- update(m2, . ~ . - tInOp:tOut)
# drop1(m2, test = "Chisq")
# m2 <- update(m2, . ~ . - tInOp:sex)
# drop1(m2, test = "Chisq")
# m2 <- update(m2, . ~ . - tInOp)
# drop1(m2, test = "Chisq")
# summary(m2)
# plot_residuals(m2)


#### Full model
m3 = glm(clo ~ sex*poly(tInOp, 3) + sex*poly(tOut, 3) + sex +
           poly(tInOp, 2)*poly(tOut, 2),
         family = Gamma(link = "log"))
summary(m3)
# plot_residuals(m3)

########### log link ###############################
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



################## Linear model #########################33
lm1 <- lm(sqrt(clo) ~ poly(tInOp, 1) + poly(tOut, 1) + sex)
summary(lm1)
par(mfrow = c(2, 2))
plot(lm1)
shapiro.test(lm1$residuals)

##### sqrt transformation ##########
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

#######################################################
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

#################### Gamma full ###############################
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

##### sqrt transformation ##########
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
results_models <- cbind(
  c(AIC(m4), AIC(m4_inverse), AIC_sqrt_original, AIC_log_original),
  c(BIC(m4), BIC(m4_inverse), BIC_sqrt_original, BIC_log_original)
)
rownames(results_models) <- c("Gamma_log", "Gamma_inverse", "lm_sqrt", "lm_log")
colnames(results_models) <- c("AIC", "BIC")
# we should pick Gamma_inverse model both by AIC, BIC

summary(m4_inverse)
par(mfrow = c(2, 2))
plot(m4_inverse)

m4_subject <- glm(clo ~ poly(tInOp, 2) + poly(tOut, 2) + subjId + subjId:tOut +
                    tInOp:tOut, family = Gamma(link = "inverse"))
summary(m4_subject)
drop1(m4_subject, test = "LRT")
AIC(m4_subject) # much better AIC

# ################### Profile Likelihood ####################
# model_optimization <- lm2
# 
# design_matrix = model.matrix(model_optimization)
# n = dim(design_matrix)[1]
# p = dim(design_matrix)[2]
# y = sqrt(clo)
# 
# objective = function(theta){
#   y_hat = design_matrix %*% theta[1:p]
#   Sigma_weighted = diag(n)
#   diag(Sigma_weighted)[design_matrix[, "sexmale"] == 0] = 1/theta[p+2]
#   # print(diag(Sigma_weighted))
#   # diag(Sigma_weighted)[design_matrix[, "LABUSA"] == 1] = (1 - theta[p+2])
#   Sigma_weighted = theta[p+1]*Sigma_weighted
#   result = sum(dmvnorm(y, mean = y_hat, sigma = Sigma_weighted, log = TRUE))
#   return(-result)
# }
# theat_initial = c(model_optimization$coefficients, "sigma_squared" = 1,
#                   "weight" = 1)
# opt_profile <- nlminb(theat_initial, objective)
# opt_profile$par
# # 
# predictions <- y - design_matrix %*% opt_profile$par[1:p]
# 
# par(mfrow = c(1, 1))
# hist(predictions[df$sex == "female"], col = "red", breaks = breaks_number,
#      ylim = c(0, upper_limit))
# hist(predictions[df$sex != "female"], add=T, breaks = breaks_number,
#      ylim = c(0, upper_limit), col = rgb(0, 0, 0.5, alpha = 0.1))
# 
# 
# ###############################################
# # reestimation for wegiths
# weights_vector = rep(1, dim(df)[1])
# weights_vector[df$sex == "female"] = as.numeric(opt_profile$par["weight"])
# 
# lm2_weights <- lm(sqrt(clo) ~ poly(tInOp, 2) + poly(tOut, 2) + sex +
#             tInOp:tOut, weights = weights_vector)
# summary(lm2_weights)
# 
# par(mfrow = c(2, 2))
# plot(lm2_weights)
# par(mfrow = c(1, 1))
# hist(lm2_weights$residuals[df$sex == "female"], col = "red", breaks = breaks_number,
#      ylim = c(0, upper_limit))
# hist(lm2_weights$residuals[df$sex != "female"], add=T, breaks = breaks_number,
#      ylim = c(0, upper_limit), col = rgb(0, 0, 0.5, alpha = 0.1))
# 
# 
# ########################33
# hessian_matrix = hessian(objective,opt$par)
# standard_error = sqrt(diag(solve(hessian_matrix)))
# 
# c(opt$par["weight"] - qnorm(0.975)*standard_error[length(standard_error)],
#   opt$par["weight"] + qnorm(0.975)*standard_error[length(standard_error)])
# 
# ### PROFILE LIKELIHOOD CI
# profile_objective <- function(weight){
#   fun.tmp <- function(theta_inner, weight_input){
#     objective(c(theta_inner, weight_input))
#   }
#   theat_initial_inner_opt = c(rep(1, length(model_optimization$coefficients)),
#                               "sigma_squared" = 1)
#   nlminb(theat_initial_inner_opt, fun.tmp, weight_input = weight)$objective
# }
# profile_objective(10)
# 
# p1 <- seq(-0.5, opt$par["weight"]*2.7, by = 0.01)
# logLp1 <- sapply(p1, profile_objective) ## note sapply!
# logLp1 <- logLp1 - min(logLp1) ## normalization
# L_CI_lower = min(p1[exp(-logLp1) > exp(-qchisq(0.95,df=1)/2)])
# L_CI_upper =  max(p1[exp(-logLp1) > exp(-qchisq(0.95,df=1)/2)])
# 
# obs_hessian <- hessian(profile_objective, opt$par["weight"])[1, 1]
# quadratic_approximation <-
#   exp( -0.5* obs_hessian * (p1 - opt$par["weight"])^2)
# quadratic_approximation <- quadratic_approximation / max(quadratic_approximation)
# 
# par(mfrow = c(1,1))
# plot(p1, exp(-logLp1), type = "l",
#      xlab="Weight", ylab="Profile Likelihood",
#      main="The comparison between Wald's and Likelihood based Confidence Intervals")
# axis(side=1, at=seq(0, 10, by=1))
# lines(p1, rep(exp(-qchisq(0.95,df=1)/2), length(p1)), col = 2)
# rug(L_CI_lower, ticksize = 0.1, lwd = 2, col = "red")
# rug(L_CI_upper, ticksize = 0.1, lwd = 2, col = "red")
# c(L_CI_lower, L_CI_upper)
# lines(p1, quadratic_approximation, col = "blue")
# abline(v = opt$par["weight"], lty = 2)
# rug(opt$par["weight"] - qnorm(0.975)*standard_error[length(standard_error)],
#     ticksize = 0.1, lwd = 2, col = "blue")
# rug(opt$par["weight"] + qnorm(0.975)*standard_error[length(standard_error)],
#     ticksize = 0.1, lwd = 2, col = "blue")
# legend("topright", 95,
#        legend=c("Profile likelihood", "Quadratic approximation",
#                 "95% confidence interval"),
#        col=c("black", "blue", "red"), lty = 1:1, cex=0.8,
#        inset = 0.02)
# 
# c(opt$par["weight"] - qnorm(0.975)*standard_error[length(standard_error)],
#   opt$par["weight"] + qnorm(0.975)*standard_error[length(standard_error)])
# c(L_CI_lower, L_CI_upper)
# 
# wald_statistic = (opt$par["weight"] - 1)/standard_error[length(standard_error)]
# 
# # LRT
# ll_full <- -profile_objective(as.numeric(opt$par["weight"]))
# ll_test <- -profile_objective(1)
# LRT <- -2*(ll_test - (ll_full))
# p <- 1 - pchisq(LRT, df = 1)
# p
# 
# 
# 
# lm2_testing <- lm(sqrt(clo) ~ tInOp + I(tInOp^2) + sex:tInOp + sex:I(tInOp^2) + 
#                     tOut + I(tOut^2) + sex:tOut + sex:I(tOut^2) + 
#                     sex + tInOp:tOut)
# summary(lm2_testing)
# par(mfrow = c(2, 2))
# plot(lm2_testing)
# 
# lm3 <- lm(sqrt(clo) ~ poly(tInOp, 2) + poly(tOut, 2) + sex +
#             tInOp:tOut)
# summary(lm2)
# par(mfrow = c(2, 2))
# plot(lm2)
# 
# 
# 
# 
# lm2 <- lm(sqrt(clo) ~ sex*poly(tInOp, 2) + poly(tOut, 3) + sex + subjId +
#             subjId*poly(tInOp, 2) + subjId*poly(tOut, 2) + tInOp:tOut)
# summary(lm2)
# par(mfrow = c(2, 2))
# plot(lm2)
# shapiro.test(lm2$residuals)
# par(mfrow = c(1, 1))
# plot(lm2$residuals)
# 
# #### model reduction, comparison between models
# clo_female <- as.vector(df %>% filter(sex == "female") %>% select(clo))$clo
# clo_male <- as.vector(df %>% filter(sex == "male") %>% select(clo))$clo
# par(mfrow = c(1, 2))
# boxcox(clo_female ~ 1)
# boxcox(clo_male ~ 1)

########################################################################3
drop1(m3, test = "Chisq")
anova(m3, m3_reduction, test = "LRT")

m3 <- update(m3, . ~ . - sex:poly(tInOp, 3))
drop1(m3, test = "Chisq")
m3 <- update(m3, . ~ . - poly(tInOp, 3))
drop1(m3, test = "Chisq")
summary(m3)
plot_residuals(m3)


#### Including weights
weight_female = 100
weights_vector <- rep(1, dim(df)[1])
weights_vector[df$sex == "female"] <- weight_female

############################################
m1 = glm(clo ~ tInOp + tOut + sex, family = Gamma(link = "inverse"),
         data = df)
summary(m1)
plot_residuals(m1)
par(mfrow = c(2, 2))
plot(m1)

X = model.matrix(m1)
theta_initial = c(0.5, 0.5, 0.1, 0.5, 0.05)

### own optimization
gamma_own <- function(theta){
  mean = X %*% theta[1:dim(X)[2]]
  shape_parameter = theta[length(theta)]
  scale_parameter = mean/shape_parameter
  return(-sum(dgamma(1/df$clo, shape = shape_parameter,
                     scale = scale_parameter, log = TRUE)))
}
gamma_own(theta_initial)
opt <- nlminb(theta_initial, gamma_own)
opt


# model_optimization <- m1
weights_vector <- rep(1, dim(df)[1])
objective_weights = function(weight_female_function){
  # weights_vector <- rep(exp(weight_female_function), dim(df)[1])
  # weights_vector <- rep(1/weight_female_function, dim(df)[1])
  weights_vector[df$sex == "female"] <- exp(weight_female_function)
  # weights_vector <- weights_vector 
  # weights_vector = exp(weight_female_function)
  m1_tmp = glm(clo ~ tInOp + tOut + sex, family = Gamma(link = "inverse"),
               data = df, weights = weights_vector)
  return(-as.numeric(logLik(m1_tmp)) / sum(weights_vector))
}
objective_weights(1)
opt <- nlminb(4, objective_weights)
exp(opt$par)

weights_vector <- rep(1, dim(df)[1])
candidates <- seq(-2, 5, length.out = 100)
AIC_list <- rep(NA, length(candidates))
for (i in 1:length(candidates)){
  weights_vector[df$sex == "female"] <- exp(candidates[i])
  # weights_vector <- rep(exp(candidates[i]), dim(df)[1])
  # m1_checking <- update(m1,
  #                       weights = weights_vector)
  m1_tmp = glm(clo ~ tInOp + tOut + sex, family = Gamma(link = "inverse"),
               data = df, weights = weights_vector)
  AIC_list[i] = AIC(m1_tmp) / sum(weights_vector)
}


# weights_vector[df$sex == "female"] <- exp(2)
weights_vector <- rep(5, dim(df)[1])
m1_tmp = glm(clo ~ tInOp + tOut + sex, family = Gamma(link = "inverse"),
             data = df, weights = weights_vector)
summary(m1_tmp)

glmGamma <- glm(clo ~ tInOp, family = Gamma(link = "inverse"), data = df)
summary(glmGamma)

library(MASS)
myshape <- gamma.shape(glmGamma)
gamma.dispersion(glmGamma)
gampred <- predict(glmGamma , type = "response", se = TRUE, 
                   dispersion = 1/myshape$alpha) 
summary(glmGamma, dispersion = 1/myshape$alpha)

###############
fit2 <- vglm(clo ~ tInOp + tOut + sex, gamma2, data = df,
             trace = TRUE, crit = "coef")
coef(fit2, matrix = TRUE)

##### data generation stuff
binary <- sample(c(rep(1, 50), rep(0, 50)))
x_small <- 
  x_large <- 
  
  
  
  weights_vector[df$sex == "female"] <- exp(-0.5)
m1_checking <- update(m1,
                      weights = weights_vector)
summary(m1)
summary(m1_checking)
c(AIC(m1), AIC(m1_checking))
par(mfrow = c(2, 2))
plot(m1)
plot_residuals(m1)

  
  
  
  
  
  
  
  






