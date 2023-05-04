# Preamble #####################################################################


### File Description ###########################################################
##
##   Søren Skjernaa - s223316
##   31/03-2023
##   
##   Advanced Dataanalysis and Statistical Modelling
##   Assignment 2
##   
##   Note: The data file for the analysis is assumed to be found at the relative
##         paths "data/clothing.csv" and "data/earinfect.txt".
##
##
###############################################################################¤


### Clean up ###################################################################
rm(list=ls()) 
print(utils::getSrcDirectory(function(){}))
print(utils::getSrcFilename(function(){}, full.names = TRUE))
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()
Sys.setenv(LANG = "en")


### Libraries ##################################################################

# Plotting
library(ggplot2)     # Nice plots
library(patchwork)   # Layout of plots
library(gridExtra)

# Handling data frames
library(tidyr)       # Reshape data frames
library(dplyr)       # Rename column in data frames

# Statistics
library(nlme)
library(betareg)

# Post-hoc analysis
library(emmeans)
library(multcomp)

# Model diagnostics
library(boot)
library(MASS)
library(MESS)
library(nortest)
library(influence.ME)

# Other
library(xtable)      # Latex tables


# _ ############################################################################
# Part B: Ear infections #######################################################


### Load the data and get initial overview #####################################

# Load data
df_data <- read.csv("data/earinfect.txt", sep = " ", header = TRUE) 
str(df_data)
df_data <- mutate_at(df_data, c("swimmer", "location", "age", "sex"), 
                     as.factor)
summary(df_data)
attach(df_data)


# Table data to get picture if balanced
table(swimmer, location)
table(swimmer, age)
table(swimmer, sex)
table(location, age) 
table(location, sex)
table(age, sex)


### Plot of the data ###########################################################

# We plot the rates of infections for each group, to have a comparable Y-value
df_data$rate <- df_data$infections / df_data$persons 
comb <- t(combn(colnames(df_data)[1:4], 2))

plot_list <- list()
for (i in 1:nrow(comb)){
  plot_list[[i]] <- ggplot(df_data, 
                           aes_string(x=comb[i,1], y="rate", col=comb[i,2])) +
    geom_boxplot() +
    labs(x = comb[i,1],
         y = "Rate",
         col = comb[i,2]) +
    theme(plot.title = element_text(hjust=0.5, size=14, face="bold"),
          axis.title = element_text(size=12, face="bold"),
          legend.title = element_text(size=12, face="bold",),
          axis.text = element_text(size=12, face="bold"))
  
  
}
grid.arrange(grobs = plot_list, nrow = 3, ncol = 2)

df_data = subset(df_data, select = -c(rate))
rm(list=setdiff(ls(), "df_data"))


### Fit Poisson model ##########################################################

# Second order interaction model
m1 <- glm(infections ~ (swimmer + location + age + sex)^2, 
          family = poisson(link = "log"), offset = log(persons),
          data = df_data)
summary(m1)
par(mfrow = c(2, 2))
plot(m1)

# m1_quassi <- glm(infections ~ (swimmer + location + age + sex)^2, 
#           family = quasipoisson, offset = log(persons),
#           data = df_data)
# summary(m1_quassi)
# par(mfrow = c(2, 2))
# plot(m1_quassi)

# Full model
mf <- glm(infections ~ swimmer*age*sex*location, 
          family = poisson(link = "log"), offset = log(persons),
          data = df_data)
summary(mf)
par(mfrow = c(2, 2))
plot(mf)

# Initial Goodnes of fit test
dev_mf <- summary(mf)$deviance
df_mf <- summary(mf)$df.resid
p_mf <- 1 - pchisq(dev_mf, df_mf)
p_mf


### Model diagnostics ##########################################################

# Initial Goodnes of fit test
dev <- summary(m1)$deviance
df <- summary(m1)$df.resid
p <- 1 - pchisq(dev, df)
p # p = 0.98 so the data do not reject the initial model.

pred <- predict(model_final, type="response", 
                inteval="confidence",se.fit=TRUE)

# pred_2 <- predict(m1, type="link", inteval="cofidence",se.fit=TRUE)
# pred_2

index <- seq(1, length(pred$fit), 1)
par(mfrow = c(1, 1))
plot(pred$fit, ylim = c(0, 25), lwd = 2,
     xlab="Index", ylab="Number of ear infections",
     main="Fitted values vs observations with the measure of uncertainty")
grid()
segments(x0=index,x1=index,
         y0=pred$fit - qnorm(0.975)*pred$se.fit,
         y1=pred$fit + qnorm(0.975)*pred$se.fit, col="red", lty = 2, lwd = 2)
lines(df_data$infections, type = "p", col = "blue", lwd = 2)
legend("topright", 95,
       legend=c("Fitted values", "Observations",
                "95% confidence interval"),
       col=c("black", "blue", "red"), lty = c(1, 1, 2), cex=0.8,
       inset = 0.02)

## Model diagnostic plots ##
par(mfrow=c(2,2))
plot(m1, which = 1:4)
par(mfrow=c(1,1))




### Model reduction ############################################################

# Reduction of second order interaction model
drop1(m1, test = "Chisq")

m1 <- update(m1, . ~ . - swimmer:sex)
drop1(m1, test = "Chisq")

m1 <- update(m1, . ~ . - swimmer:age)
drop1(m1, test = "Chisq")

m1 <- update(m1, . ~ . - location:age)
drop1(m1, test = "Chisq")

m1 <- update(m1, . ~ . - age:sex)
drop1(m1, test = "Chisq")

m1 <- update(m1, . ~ . - age)
drop1(m1, test = "Chisq")

m1 <- update(m1, . ~ . - swimmer:location)
drop1(m1, test = "Chisq")

m1 <- update(m1, . ~ . - swimmer:location)
drop1(m1, test = "Chisq")

m1 <- update(m1, . ~ . - swimmer)
drop1(m1, test = "Chisq")

m1 <- update(m1, . ~ . - location:sex)
drop1(m1, test = "Chisq")

m1 <- update(m1, . ~ . - sex)

### Model reduction FULL ############################################################

# Reduction of second order interaction model
drop1(mf, test = "Chisq")
mf <- update(mf, . ~ . - swimmer:age:sex:location)
drop1(mf, test = "Chisq")
mf <- update(mf, . ~ . - swimmer:age:sex)
drop1(mf, test = "Chisq")
mf <- update(mf, . ~ . - swimmer:age:location)
drop1(mf, test = "Chisq")
mf <- update(mf, . ~ . - swimmer:sex:location)
drop1(mf, test = "Chisq")
mf <- update(mf, . ~ . - swimmer:sex)
drop1(mf, test = "Chisq")
mf <- update(mf, . ~ . - swimmer:age)
drop1(mf, test = "Chisq")
mf <- update(mf, . ~ . - age:sex:location)
drop1(mf, test = "Chisq")
mf <- update(mf, . ~ . - age:location)
drop1(mf, test = "Chisq")
mf <- update(mf, . ~ . - age:sex)
drop1(mf, test = "Chisq")
mf <- update(mf, . ~ . - age)
drop1(mf, test = "Chisq")
mf <- update(mf, . ~ . - swimmer:location)
drop1(mf, test = "Chisq")
mf <- update(mf, . ~ . - swimmer)
drop1(mf, test = "Chisq")
mf <- update(mf, . ~ . - sex:location)
drop1(mf, test = "Chisq")
mf <- update(mf, . ~ . - sex)
round(drop1(mf, test = "Chisq"), 3)

summary(mf)
mf <- update(mf, . ~ . - 1)
round(summary(mf), 3)

model_final <- mf
summary(model_final)

dev <- summary(model_final)$deviance
df <- summary(model_final)$df.resid
p <- 1 - pchisq(dev, df)
p # p = 0.98 so the data do not reject the initial model.

# Plot of deviance residuals against the different factors
stddevresid <- rstandard(model_final, type = "deviance")
par(mfrow=c(2,2))
plot(swimmer, stddevresid, xlab="", ylab = "Deviance residuals")
abline(h=0, col="red")
plot(location, stddevresid, xlab="", ylab = "Deviance residuals")
abline(h=0, col="red")
plot(age, stddevresid, xlab="", ylab = "Deviance residuals")
abline(h=0, col="red")
plot(sex, stddevresid, xlab="", ylab = "Deviance residuals")
abline(h=0, col="red")
par(mfrow=c(1,1))


round(confint(mf), 3)
round(confint.default(mf), 3)

par(mfrow = c(2, 2))
plot(mf)

# # Control of result
# test <- glm(infections ~ swimmer + age + sex + location +
#                         swimmer:age + swimmer:sex + swimmer:location + 
#                         age:sex + age:location, 
#                       family = poisson(link = "log"), offset = log(persons),
#                       data = df_data)
# tdev <- summary(test)$deviance
# tdf <- summary(test)$df.resid
# d <- tdev - dev
# 1 - pchisq(d, 1)


### Post-Hoc analysis ##########################################################

summary(m1)
m1 <- update(m1, . ~ . - 1)
summary(m1)


# Confidence intervals
confint(m1)
