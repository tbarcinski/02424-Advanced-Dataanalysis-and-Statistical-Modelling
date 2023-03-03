rm(list=ls()) 
print(utils::getSrcDirectory(function(){}))
print(utils::getSrcFilename(function(){}, full.names = TRUE))
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

library(MASS)
library(dplyr)
library(nortest)
library(GGally) 
library(ggplot2)
library(corrplot)
library(tree)
library(car)
library(patchwork)
library(mvtnorm)
library(nlme)
?qqPlot

dat <- read.table("dioxin.csv", sep=',', head=TRUE)
df <- data.frame(dat)
summary(df)
sum(is.na(df))
# do the imputation of the missing values
# examine and comment the full model vs the model with excluded NA variables

cols <- c("PLANT", "TIME", "LAB", "OXYGEN", "LOAD", "PRSEK")
df[cols] <- lapply(df[cols], as.factor)
str(df)
attach(df)

## checking for normality
boxcox <- boxcox(DIOX ~ 1)
qqnorm(log(df$DIOX))
qqline(log(df$DIOX), col = "red")
shapiro.test(log(df$DIOX))
lillie.test(log(df$DIOX))
# wallyplot R residuals

### EDA
# ggpairs(df) 
df_part <- df[,c("DIOX", "O2COR", "NEFFEKT", "QRAT")]
df_part["diox_trans"] <- log(df_part$DIOX)
(p <- pairs(df_part, panel=panel.smooth))
corrplot(cor(df_part))
cor(df_part)

# cor_matrix <- cor(df)
ggplot(df, aes(x = O2COR, y = DIOX)) +
  geom_point(aes(col = PLANT)) +
  facet_wrap(~LAB)

ggplot(df, aes(x = O2COR, y = log(DIOX))) +
  geom_point(aes(col = PLANT)) +
  facet_wrap(~LAB)

ggplot(df, aes(x = PLANT, y = DIOX)) + geom_boxplot() +
  facet_wrap(~LAB)

##### analysis of stuff ###########
df %>% group_by(TIME, PLANT) %>% summarise(n = n())
df %>% group_by(LAB, PLANT) %>% summarise(n = n())
# there are certian mesaurments that were sent just to one lab
df %>% filter(PLANT == "KARA")
df %>% filter(PLANT == "RENO_S") %>% 
  select(DIOX, LAB, OXYGEN, LOAD, PRSEK, O2COR, NEFFEKT, QRAT,
         TIME) %>% 
  filter(LAB == "KK")

#### Analysis of time ######
active_variables <- c("O2COR", "NEFFEKT", "QRAT")
p1 <- df %>% filter(PLANT == "RENO_N") %>% 
  ggplot(., aes(x = O2COR, y = log(DIOX), col = TIME)) +
  geom_point() + facet_wrap(~LAB)
p2 <- df %>% filter(PLANT == "RENO_N") %>% 
  ggplot(., aes(x = NEFFEKT, y = log(DIOX), col = TIME)) +
  geom_point() + facet_wrap(~LAB)
p3 <- df %>% filter(PLANT == "RENO_N") %>% 
  ggplot(., aes(x = QRAT, y = log(DIOX), col = TIME)) +
  geom_point() + facet_wrap(~LAB)
p1/p2/p3


model = tree(diox_trans~., data = df_part[, -1])
plot(model)
text(model)

# co, so2, prsek, we should keep them
# pm <- ggpairs(tips, mapping = aes(color = sex),
#               columns = c("total_bill", "time", "tip"))
## First simple model
lm <- lm(diox_trans ~. - DIOX, data = df_part)
summary(lm)

par(mfrow = c(2, 2))
plot(lm, which = 1:4)
par(mfrow = c(1, 1))
# 16 outlier???
# quadratic stuff in the residuals

# To do
# 1) interactions plots between continous and factor variables
# 2) interactions plots
# 3) box plots of the block variables for the diox response
# 4) Do the plots with factors as in the lecture
# MSW plants
# interactions

######## 1) SIMPLE ADDITIVE MODEL ############
model_1 <- lm(log(DIOX) ~ OXYGEN + LOAD + PRSEK + PLANT + TIME + LAB, data = df)
summary(model_1)
par(mfrow = c(2, 2))
plot(model_1)

acf(residuals(model_1))
rstudent(model_1)[which.max(abs(rstudent(model_1)))]
# 13
drop1(model_1, test = "F")

model_2 <- update(model_1,~. -PRSEK, data = df)
anova(model_2, model_1, test="Chisq")

## wykresy częściowych rezydów modyfikowanych względem zmiennych


### CHECKING OF ASSUMPTIONS

## outliers
par(mfrow = c(1, 1))
plot(rstudent(model_1), ylim = c(-3, 3))
abline(h = 0)
abline(h = -2)
abline(h = 2)
shapiro.test(residuals(model_1)) # well behaved residuals

######## 2) SIMPLE ADDITIVE MODEL WITH CONTINOUS VARIABLES #####
model_c1 <- lm(log(DIOX) ~ O2COR + NEFFEKT + QRAT + PLANT + TIME + LAB)
summary(model_c1)
par(mfrow = c(2, 2))
plot(model_c1)

drop1(model_c1, test = "F")

model_c2 <- update(model_c1, ~. - QRAT)
summary(model_c2)
anova(model_c1, model_c2, test = "Chisq")
df[13:16, ]
### reduction

rstudent(model_c2)[which.max(abs(rstudent(model_c2)))]
par(mfrow = c(1, 1))
plot(rstudent(model_c2), ylim = c(-3, 3))
abline(h = 0)
abline(h = -2)
abline(h = 2)
shapiro.test(residuals(model_c2)) 

#### Predictions
# can we just exp the residuals

# interpretation
exp(model_c2$coefficients$O2COR)
exp(0.18)
exp(2.9)

#### INTERACTIONS
model_c2_full <- lm(log(DIOX) ~ (O2COR + NEFFEKT + QRAT + PLANT + TIME + LAB)^2)
summary(model_c2_full, cor = T)

## with plant
## with labolatories


#### QUADRATIC TERMS
model_quadratic_1 <- lm(log(DIOX) ~ polym(O2COR, NEFFEKT, QRAT, 
                                          degree =2))
summary(model_quadratic_1)
par(mfrow = c(2, 2))
plot(model_1)


## An odd looking data-point? 13
sort(DIOX)
df[12:14, ]
n<-length(DIOX)

pairs(df[, c("DIOX", "O2COR", "NEFFEKT", "QRAT")],
      panel=panel.smooth,
      pch=c(rep(1,12),19,rep(1,n-13)),
      col=c(rep(1,12),2,rep(1,n-13)),
      cex=c(rep(1,12),2,rep(1,n-13)))

pairs(df_part[, c("diox_trans", "O2COR", "NEFFEKT", "QRAT")],
      panel=panel.smooth,
      pch=c(rep(1,12),19,rep(1,n-13)),
      col=c(rep(1,12),2,rep(1,n-13)),
      cex=c(rep(1,12),2,rep(1,n-13)))
####### 3) PREDICTION, DEPENDENCE ON OPERATING CONDITONS, PLANTS #####

######## LIKELIHOOD ESTIMATION OF WEIGHTS #################
### testing stuff
par(mfrow = c(1, 1))
model_test = lm(Volume ~ Girth + Height, data = trees)
n = dim(trees)[1]
p = dim(trees)[2]
y = trees$Volume

design_matrix = model.matrix(model_test)

objective <- function(theta){
  y_hat = design_matrix %*% theta[1:p]
  sum = 0
  for(i in 1:n){
    sum = sum + dnorm(y[i], mean = y_hat[i], sd = sqrt(theta[p+1]),log=TRUE)
  }
  return(-sum)
}
theta = c("theta_0" = 0, "theta_1" = 0, "theta_3" = 0, "sigma_squared" = 1)

opt <- nlminb(theta,objective)
opt$par

objective_2 <- function(theta){
  y_hat = design_matrix %*% theta[1:p]
  result = sum(dmvnorm(y, mean = y_hat, sigma = theta[p+1]*diag(n), log = TRUE))
  return(-result)
}
theta = c("theta_0" = 0, "theta_1" = 0, "theta_3" = 0, "sigma_squared" = 1)

opt_2 <- nlminb(theta, objective_2)
opt_2$par

#### doing it on our data



