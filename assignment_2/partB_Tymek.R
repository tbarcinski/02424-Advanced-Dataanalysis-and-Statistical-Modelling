# File Description -------------------------------------------------------------
#
#   Tymoteusz Barcinski - s221937
#   
#   Advanced Dataanalysis and Statistical Modelling
#   Assignment 2 - Part B
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

# Reading data
dat <- read.table("earinfect.txt", sep=' ', head=TRUE)
cols <- c("swimmer", "location", "age", "sex")
dat[cols] <- lapply(dat[cols], as.factor)
df <- data.frame(dat)
sum(is.na(df)) # no missing values
attach(df)

#### EDA -----------------------------------------------------------------------
table(swimmer)
table(location)
table(age)
table(sex)

table(swimmer, location)
table(swimmer, age)
table(swimmer, sex)
# balanced data in terms of study

#### Modelling -----------------------------------------------------------------
m1 = glm(infections ~ swimmer + location + age + sex, offset = log(persons),
          family = poisson(link = "log"))
summary(m1, cor = TRUE)

par(mfrow = c(2, 2))
plot(m1, which = 1:4)
(p_val <- 1 - pchisq(m1$deviance, df = m1$df.residual))
# 0.97 and hence the model does give a good fit to the data
# model reduction 
drop1(m1, test = "Chisq")
m1 <- update(m1, . ~ . - sex)
drop1(m1, test = "Chisq")
m1 <- update(m1, . ~ . - age)
drop1(m1, test = "Chisq")
m1 <- update(m1, . ~ . - swimmer)
drop1(m1, test = "Chisq")

# model with only location
summary(m1)
par(mfrow = c(2, 2))
plot(m1)

pearson <- residuals(m1, type = "pearson")
deviance <- residuals(m1, type = "deviance")
par(mfrow = c(1, 2))
plot(pearson)
plot(deviance)
persons - deviance # almost the same

pred <- predict(m1, type="response", 
                inteval="confidence",se.fit=TRUE)

pred_2 <- predict(m1, type="link", inteval="cofidence",se.fit=TRUE)
pred_2

index <- seq(1, length(pred$fit), 1)
par(mfrow = c(1, 1))
plot(pred$fit, ylim = c(0, 25))
segments(x0=index,x1=index,
         y0=pred$fit - qnorm(0.975)*pred$se.fit,
         y1=pred$fit + qnorm(0.975)*pred$se.fit, col="red", lty = 2)
lines(df$infections, type = "p", col = "green")


######## trying to include some interactions ... ##################
m2 = glm(infections ~ swimmer + location + age + sex +
           sex:location, offset = log(persons),
         family = poisson(link = "log"))
summary(m2)
(p_val <- 1 - pchisq(m2$deviance, df = m2$df.residual))

drop1(m2, test = "Chisq")
m2 <- update(m2, . ~ . - age)
drop1(m2, test = "Chisq")
m2 <- update(m2, . ~ . - location:sex)
drop1(m2, test = "Chisq")
m2 <- update(m2, . ~ . - sex)
drop1(m2, test = "Chisq")
m2 <- update(m2, . ~ . - swimmer)
drop1(m2, test = "Chisq")
summary(m2)

m3 <- glm(infections ~ swimmer + location + age + sex +
            swimmer:location + sex:location + age:location, offset = log(persons),
          family = poisson(link = "log"))
summary(m3)
(p_val <- 1 - pchisq(m3$deviance, df = m3$df.residual))
par(mfrow = c(2, 2))
plot(m3)

drop1(m3, test = "Chisq")
m3 <- update(m3, . ~ . - location:age)
drop1(m3, test = "Chisq")
m3 <- update(m3, . ~ . - age)
drop1(m3, test = "Chisq")
m3 <- update(m3, . ~ . - swimmer:location)
drop1(m3, test = "Chisq")
m3 <- update(m3, . ~ . - swimmer)
drop1(m3, test = "Chisq")
m3 <- update(m3, . ~ . - location:sex)
drop1(m3, test = "Chisq")
m3 <- update(m3, . ~ . - sex)

summary(m3)


#### RANDOM STUFF ###########
# m1_quasi = glm(infections ~ swimmer + location + age + sex, offset = log(persons),
#          family = quasipoisson(link = "log"))
# summary(m1_quasi)
# par(mfrow = c(2, 2))
# plot(m1_quasi, which = 1:4)
# ### polynomials ???
# m1_sqrt = glm(infections ~ swimmer + location + age + sex, offset = log(persons),
#          family = poisson(link = "sqrt"))
# summary(m1_sqrt)
# 
# par(mfrow = c(2, 2))
# plot(m1_sqrt, which = 1:4)
# (p_val <- 1 - pchisq(m1_sqrt$deviance, m1_sqrt$df.residual))

scope <- ~. + swimmer*location*age*sex
add1(m1, scope, test = "Chisq")


m2 = glm(infections ~ swimmer*location*age + sex, offset = log(persons),
         family = poisson(link = "log"))
summary(m2)
(p_val <- 1 - pchisq(m2$deviance, m2$df.residual))

drop1(m2, test = "Chisq")
m2 <- update(m2, . ~ . - sex)
drop1(m2, test = "Chisq")
m2 <- update(m2, . ~ . - swimmer:location:age)
drop1(m2, test = "Chisq")
m2 <- update(m2, . ~ . - location:age)
drop1(m2, test = "Chisq")
m2 <- update(m2, . ~ . - swimmer:age)
drop1(m2, test = "Chisq")
m2 <- update(m2, . ~ . - age)
drop1(m2, test = "Chisq")
m2 <- update(m2, . ~ . - swimmer:location)
drop1(m2, test = "Chisq")
m2 <- update(m2, . ~ . - swimmer)
drop1(m2, test = "Chisq")

summary(m2)
(p_val <- 1 - pchisq(m2$deviance, m2$df.residual))

########### trying more #################
m4 = glm(infections ~ sex + swimmer*location*age, offset = log(persons),
         family = poisson(link = "log"))
summary(m4)
plot(m4)
 (p_val <- 1 - pchisq(m4$deviance, df = m4$df.residual))



drop1(m4, test = "Chisq")
m4 <- update(m4, . ~ . - sex)
drop1(m4, test = "Chisq")
m4 <- update(m4, . ~ . - swimmer:location:age)
drop1(m4, test = "Chisq")
m4 <- update(m4, . ~ . - location:age)
drop1(m4, test = "Chisq")
m4 <- update(m4, . ~ . - swimmer:age)
drop1(m4, test = "Chisq")
m4 <- update(m4, . ~ . - age)
drop1(m4, test = "Chisq")
m4 <- update(m4, . ~ . - swimmer:location)
drop1(m4, test = "Chisq")
m4 <- update(m4, . ~ . - swimmer)
summary(m4)

(p_val <- 1 - pchisq(m4$deviance, df = m4$df.residual))


m4 = glm(infections ~ swimmer + location*age*sex, offset = log(persons),
         family = poisson(link = "log"))
summary(m4)
drop1(m4, test = "Chisq")
m4 <- update(m4, . ~ . - location:sex:age)
drop1(m4, test = "Chisq")
m4 <- update(m4, . ~ . - age:sex)
drop1(m4, test = "Chisq")
m4 <- update(m4, . ~ . - location:age)
drop1(m4, test = "Chisq")
m4 <- update(m4, . ~ . - age)
drop1(m4, test = "Chisq")
m4 <- update(m4, . ~ . - swimmer)
drop1(m4, test = "Chisq")
m4 <- update(m4, . ~ . - location:sex)
drop1(m4, test = "Chisq")
m4 <- update(m4, . ~ . - sex)
summary(m4)

m4 = glm(infections ~ age + location*swimmer*sex, offset = log(persons),
         family = poisson(link = "log"))
summary(m4)
drop1(m4, test = "Chisq")
m4 <- update(m4, . ~ . - age)
drop1(m4, test = "Chisq")
m4 <- update(m4, . ~ . - location:swimmer:sex)
drop1(m4, test = "Chisq")
m4 <- update(m4, . ~ . - swimmer:sex)
drop1(m4, test = "Chisq")
m4 <- update(m4, . ~ . - location:swimmer)
drop1(m4, test = "Chisq")
m4 <- update(m4, . ~ . - swimmer)
drop1(m4, test = "Chisq")
m4 <- update(m4, . ~ . - sex:location)
drop1(m4, test = "Chisq")
m4 <- update(m4, . ~ . - sex)
summary(m4)

m4 = glm(infections ~ location + age*swimmer*sex, offset = log(persons),
         family = poisson(link = "log"))
summary(m4)
drop1(m4, test = "Chisq")
m4 <- update(m4, . ~ . - age:swimmer:sex)
drop1(m4, test = "Chisq")
m4 <- update(m4, . ~ . - swimmer:sex)
drop1(m4, test = "Chisq")
m4 <- update(m4, . ~ . - sex:age)
drop1(m4, test = "Chisq")
m4 <- update(m4, . ~ . - sex)
drop1(m4, test = "Chisq")
m4 <- update(m4, . ~ . - age:swimmer)
drop1(m4, test = "Chisq")
m4 <- update(m4, . ~ . - age)
drop1(m4, test = "Chisq")
m4 <- update(m4, . ~ . - swimmer)
summary(m4)

###############################################
m5 = glm(infections ~ location, offset = log(persons),
         family = poisson(link = "log"))
m6 = glm(infections ~ location + swimmer, offset = log(persons),
         family = poisson(link = "log"))
m7 = glm(infections ~ location + swimmer + location:swimmer, offset = log(persons),
         family = poisson(link = "log"))
anova(m5, m6, m7, test = "Chisq")

###########################################



