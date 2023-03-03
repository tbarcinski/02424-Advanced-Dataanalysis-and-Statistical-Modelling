# File Description -------------------------------------------------------------
#
#   Søren Skjernaa - s223316
#   10/03-2023
#   
#   Advanced Dataanalysis and Statistical Modelling
#   Assignment 1
#   
#   Note: The data file for the analysis is assumed to be found at the relative
#         path "data/assignment_1_dioxin.csv".
#
#_______________________________________________________________________________

# Library imports --------------------------------------------------------------

# Plotting
library(ggplot2)
library(gridExtra)
# library(diagram)

# Estimation
# library(lmerTest)

# Post-hoc analysis
library(emmeans)
library(multcomp)
# library(xtable)

# Model diagnostics
library(MASS)
library(MESS)
library(nortest)
library(influence.ME)


# Functions --------------------------------------------------------------------


# . ----------------------------------------------------------------------------
# Initial data exploration -----------------------------------------------------

## Data preprocessing ----------------------------------------------------------
DATA <- read.table("data/assignment_1_dioxin.csv", sep=",", header=TRUE)

summary(DATA)
DATA$OBSERV <- factor(DATA$OBSERV)
DATA$TIME <- factor(DATA$TIME)
DATA$PLANT <- factor(DATA$PLANT)
DATA$LAB <- factor(DATA$LAB)
summary(DATA)

attach(DATA)


### Tabulating the data --------------------------------------------------------

# We have N = 57 observations

# Planned active variables
table(OXYGEN)
table(LOAD)
table(PRSEK, useNA ="ifany")

table(OXYGEN, LOAD)
table(OXYGEN, PRSEK)
table(LOAD, PRSEK)

# Block effects
table(PLANT)
table(TIME)
table(LAB)

table(PLANT, TIME)   # Only RENO_N have repeated measurements.
table(PLANT, LAB)
table(LAB, TIME)

# Measured active variables
table(O2COR, OXYGEN)
table(NEFFEKT, LOAD)
table(QRAT, PRSEK)

### Summary statistics for the outcome -----------------------------------------
mean(DIOX)
sd(DIOX)

tapply(DIOX, OXYGEN, mean)
tapply(DIOX, OXYGEN, sd)

tapply(DIOX, LOAD, mean)
tapply(DIOX, LOAD, sd)

tapply(DIOX, PRSEK, mean)
tapply(DIOX, PRSEK, sd)

tapply(DIOX, PLANT, mean)
tapply(DIOX, PLANT, sd)

tapply(DIOX, PRSEK, mean)
tapply(DIOX, PRSEK, sd)



### Other interisting summary statistics ---------------------------------------

tapply(QRAT, PLANT, mean)
tapply(QRAT, PLANT, sd)

tapply(NEFFEKT, PLANT, mean)
tapply(NEFFEKT, PLANT, sd)

### Plotting the data ----------------------------------------------------------

# Scatter plots of all variables
plot(DATA)


# Boxplot of the measured active variables vs the planned active variables
g1 <- ggplot(DATA, aes(x=OXYGEN, y=O2COR, colour=OXYGEN)) +
  geom_boxplot() +
  theme(axis.text=element_text(size=14,face="bold"),
        axis.title=element_text(size=16,face="bold"),
        legend.position="none") + 
  xlab("OXYGEN") + ylab("O2COR")

g2 <- ggplot(DATA, aes(x=LOAD, y=NEFFEKT, colour=LOAD)) +
  geom_boxplot() +
  theme(axis.text=element_text(size=14,face="bold"),
        axis.title=element_text(size=16,face="bold"),
        legend.position="none") + 
  xlab("LOAD") + ylab("NEFFEKT")

g3 <- ggplot(DATA, aes(x=PRSEK, y=QRAT, colour=PRSEK)) + geom_boxplot() +
  theme(axis.text=element_text(size=14,face="bold"),
        axis.title=element_text(size=16,face="bold"),
        legend.position="none") + 
  xlab("PRSEK") + ylab("QRAT")

grid.arrange(g1, g2, g3, nrow=1, ncol=3)


# Boxplots
g1 <- ggplot(DATA, aes(x=OXYGEN, y=DIOX, colour=OXYGEN)) +
      geom_boxplot() +
      theme(axis.text=element_text(size=14,face="bold"),
            axis.title=element_text(size=16,face="bold"),
            legend.position="none") + 
      xlab("OXYGEN") + ylab("DIOXIN")

g2 <- ggplot(DATA, aes(x=LOAD, y=DIOX, colour=LOAD)) +
      geom_boxplot() +
      theme(axis.text=element_text(size=14,face="bold"),
            axis.title=element_text(size=16,face="bold"),
            legend.position="none") + 
      xlab("LOAD") + ylab("DIOXIN")

g3 <- ggplot(DATA, aes(x=PRSEK, y=DIOX, colour=PRSEK)) + geom_boxplot() +
      theme(axis.text=element_text(size=14,face="bold"),
            axis.title=element_text(size=16,face="bold"),
            legend.position="none") + 
      xlab("PRSEK") + ylab("DIOXIN")


g4 <- ggplot(DATA, aes(x=PLANT, y=DIOX, colour=TIME)) +
      geom_boxplot() +
      theme(axis.text=element_text(size=14,face="bold"),
            axis.title=element_text(size=16,face="bold"),
            legend.position="right") + 
      xlab("PLANT") + ylab("DIOXIN")

g5 <- ggplot(DATA, aes(x=LAB, y=DIOX, colour=LAB)) +
      geom_boxplot() +
      theme(axis.text=element_text(size=14,face="bold"),
            axis.title=element_text(size=16,face="bold"),
            legend.position="none") + 
      xlab("LAB") + ylab("DIOXIN")

grid.arrange(g1, g2, g3, g4, g5, nrow=2, ncol=3)


# Interaction plot for the planned active variables
par(mfrow=c(2,2))
interaction.plot(OXYGEN, LOAD, DIOX,
                 legend=TRUE, bty="n", col=1:8, xtick = TRUE)
interaction.plot(OXYGEN, PRSEK, DIOX,
                 legend=TRUE, bty="n", col=1:8, xtick = TRUE)
interaction.plot(LOAD, PRSEK, DIOX,
                 bty="n", col=1:3, xtick = TRUE)
par(mfrow=c(1,1))

# Scatter plots of the measured active variables (investigate different slope for PLANT)
g1 <- ggplot(DATA, aes(x=O2COR, y=log(DIOX), colour=PLANT)) +
      geom_point() + 
      theme(axis.text=element_text(size=14,face="bold"),
            axis.title=element_text(size=16,face="bold"),
            legend.position="right") + 
      xlab("02COR") + ylab("DIOXIN")

g2 <- ggplot(DATA, aes(x=NEFFEKT, y=log(DIOX), colour=PLANT)) +
  geom_point() + 
  theme(axis.text=element_text(size=14,face="bold"),
        axis.title=element_text(size=16,face="bold"),
        legend.position="right") + 
  xlab("NEFFEKT") + ylab("DIOXIN")

g3 <- ggplot(DATA, aes(x=QRAT, y=log(DIOX), colour=PLANT)) +
  geom_point() + 
  theme(axis.text=element_text(size=14,face="bold"),
        axis.title=element_text(size=16,face="bold"),
        legend.position="right") + 
  xlab("QRAT") + ylab("DIOXIN")

grid.arrange(g1, g2, g3, nrow=2, ncol=2)


# Scatter plots of the measured active variables (investigate different slope for LAB)
g1 <- ggplot(DATA, aes(x=O2COR, y=log(DIOX), colour=LAB)) +
  geom_point() + 
  theme(axis.text=element_text(size=14,face="bold"),
        axis.title=element_text(size=16,face="bold"),
        legend.position="right") + 
  xlab("02COR") + ylab("DIOXIN")

g2 <- ggplot(DATA, aes(x=NEFFEKT, y=log(DIOX), colour=LAB)) +
  geom_point() + 
  theme(axis.text=element_text(size=14,face="bold"),
        axis.title=element_text(size=16,face="bold"),
        legend.position="right") + 
  xlab("NEFFEKT") + ylab("DIOXIN")

g3 <- ggplot(DATA, aes(x=QRAT, y=log(DIOX), colour=LAB)) +
  geom_point() + 
  theme(axis.text=element_text(size=14,face="bold"),
        axis.title=element_text(size=16,face="bold"),
        legend.position="right") + 
  xlab("QRAT") + ylab("DIOXIN")

grid.arrange(g1, g2, g3, nrow=2, ncol=2)

# Scatter plots of the measured active variables (investigate different slope for TIME)
g1 <- ggplot(DATA, aes(x=O2COR, y=log(DIOX), colour=TIME)) +
  geom_point() + 
  theme(axis.text=element_text(size=14,face="bold"),
        axis.title=element_text(size=16,face="bold"),
        legend.position="right") + 
  xlab("02COR") + ylab("DIOXIN")

g2 <- ggplot(DATA, aes(x=NEFFEKT, y=log(DIOX), colour=TIME)) +
  geom_point() + 
  theme(axis.text=element_text(size=14,face="bold"),
        axis.title=element_text(size=16,face="bold"),
        legend.position="right") + 
  xlab("NEFFEKT") + ylab("DIOXIN")

g3 <- ggplot(DATA, aes(x=QRAT, y=log(DIOX), colour=TIME)) +
  geom_point() + 
  theme(axis.text=element_text(size=14,face="bold"),
        axis.title=element_text(size=16,face="bold"),
        legend.position="right") + 
  xlab("QRAT") + ylab("DIOXIN")

grid.arrange(g1, g2, g3, nrow=2, ncol=2)


# . ----------------------------------------------------------------------------

# Simple additive model - planned active effects -------------------------------

### Specification of the initial model -----------------------------------------

# We insert the missing values in PRSEK
table(PRSEK, useNA ="ifany")  # Observation 15 and 16 are missing
tapply(QRAT, PRSEK, mean, useNA="ifany")  # Observations seems to have been "L"
DATA$PRSEK[15:16] <- "L"
attach(DATA)


APV1 <- lm(DIOX ~ OXYGEN  + LOAD + PRSEK + PLANT + TIME + LAB)
summary(APV1)


### Model diagnostics ----------------------------------------------------------

stdresid <- rstandard(APV1)

# Model diagnostics plots
par(mfrow=c(2, 2))
plot(APV1, which=1:4) 
par(mfrow=c(1,1))

# BOXCOX transformation
boxcox(APV1)

# # Normality plots
# par(mfrow=c(1, 2))
# hist(stdresid, main="", probability=TRUE, breaks=15)
# curve(dnorm, -3, 3, col="red", lwd=2, add=TRUE)
# plot(APV1, which=2) 
# par(mfrow=c(1, 1))
# 
# 
# # # Standardized residuals vs the three main effects
# # par(mfrow=c(1,3))
# # plot(stdresid ~ OXYGEN, col = rainbow(3)) 
# # plot(stdresid ~ LOAD, col = rainbow(3))
# # plot(stdresid ~ PRESEK, col = rainbow(3))
# # par(mfrow=c(1,1))
# # plot(stdresid ~ Cow, col = rainbow(37))
# 
# 
# # Check for outliers
# temp <- sort(stdresid)
# tail(temp, 20) 
# temp[1:20] 
# length(temp) * (2 * pnorm(-2)) # Expect around 2.5 std.resid numerically > 2.  
# remove(temp)
# 
# 
# # Check influence measures
# infl <- influence.measures(APV1)$infmat
# 
# 
# # Normality tests
# shapiro.test(stdresid)
# lillie.test(stdresid)
# cvm.test(stdresid)
# ad.test(stdresid)


### Specification of the log transformed model ---------------------------------
APV2 <- lm(log(DIOX) ~ OXYGEN  + LOAD + PRSEK + PLANT + TIME + LAB)
summary(APV2)


### Updated model diagnostics --------------------------------------------------

stdresid <- rstandard(APV2)

# Wally plot
# wallyplot(APV2)

# Model diagnostics plots
par(mfrow=c(2, 2))                                                              # Note obs 15 and 16 influential.
plot(APV2, which=1:4) 
par(mfrow=c(1,1))

# Normality plots
par(mfrow=c(1, 2))
hist(stdresid, main="", probability=TRUE, breaks=10)
curve(dnorm, -3, 3, col="red", lwd=2, add=TRUE)
plot(APV2, which=2) 
par(mfrow=c(1, 1))


# # Standardized residuals vs the three main effects
# par(mfrow=c(1,3))
# plot(stdresid ~ OXYGEN, col = rainbow(3)) 
# plot(stdresid ~ LOAD, col = rainbow(3))
# plot(stdresid ~ PRESEK, col = rainbow(3))
# par(mfrow=c(1,1))
# plot(stdresid ~ Cow, col = rainbow(37))


# BOXCOX transformation
boxcox(APV2)

# Check for outliers
temp <- sort(stdresid)
tail(temp, 20)                                                                  # 2 obs with std.resid > 2
temp[1:20]                                                                      # 0 obs with std.resid < 2
length(temp) * (2 * pnorm(-2))                                                  # Expect around 2.5 std.resid numerically > 2.  
remove(temp)


# Check influence measures
infl <- influence.measures(APV1)$infmat


# Normality tests
shapiro.test(stdresid)
lillie.test(stdresid)
cvm.test(stdresid)
ad.test(stdresid)


### Model reduction ------------------------------------------------------------

drop1(APV2, test = 'F')
APV2 <- update(APV2, . ~ . - PRSEK)

drop1(APV2, test = 'F')


### Post Hoc Analysis ----------------------------------------------------------



# . ----------------------------------------------------------------------------
# Simple additive model - measured active effects ------------------------------

### Specification of the initial model -----------------------------------------
AMV1 <- lm(DIOX ~ O2COR  + NEFFEKT + QRAT + PLANT + TIME + LAB)
summary(AMV1)


### Model diagnostics ----------------------------------------------------------

stdresid <- rstandard(AMV1)

# Wally plot
# wallyplot(AMV1)

# Model diagnostics plots
par(mfrow=c(2, 2))
plot(AMV1, which=1:4) 
par(mfrow=c(1,1))

# BOXCOX transformation
boxcox(AMV1)

# # Normality plots
# par(mfrow=c(1, 2))
# hist(stdresid, main="", probability=TRUE, breaks=15)
# curve(dnorm, -3, 3, col="red", lwd=2, add=TRUE)
# plot(AMV1, which=2) 
# par(mfrow=c(1, 1))
# 
# 
# # # Standardized residuals vs the three main effects
# # par(mfrow=c(1,3))
# # plot(stdresid ~ OXYGEN, col = rainbow(3)) 
# # plot(stdresid ~ LOAD, col = rainbow(3))
# # plot(stdresid ~ PRESEK, col = rainbow(3))
# # par(mfrow=c(1,1))
# # plot(stdresid ~ Cow, col = rainbow(37))
# 
# 
# # Check for outliers
# temp <- sort(stdresid)
# tail(temp, 20) 
# temp[1:20] 
# length(temp) * (2 * pnorm(-2)) # Expect around 2.5 std.resid numerically > 2.  
# remove(temp)
# 
# 
# # Check influence measures
# infl <- influence.measures(AMV1)$infmat
# 
# 
# # Normality tests
# shapiro.test(stdresid)
# lillie.test(stdresid)
# cvm.test(stdresid)
# ad.test(stdresid)


### Specification of the log transformed model ---------------------------------
AMV2 <- lm(log(DIOX) ~ O2COR  + NEFFEKT + QRAT + PLANT + TIME + LAB)
summary(AMV2)


### Updated model diagnostics --------------------------------------------------

stdresid <- rstandard(AMV2)

# Wally plot
# wallyplot(AMV2)

# Model diagnostics plots
par(mfrow=c(2, 2))
plot(AMV2, which=1:4) 
par(mfrow=c(1,1))

# Normality plots
par(mfrow=c(1, 2))
hist(stdresid, main="", probability=TRUE, breaks=15)
curve(dnorm, -3, 3, col="red", lwd=2, add=TRUE)
plot(AMV2, which=2) 
par(mfrow=c(1, 1))


# # Standardized residuals vs the three main effects
# par(mfrow=c(1,3))
# plot(stdresid ~ OXYGEN, col = rainbow(3)) 
# plot(stdresid ~ LOAD, col = rainbow(3))
# plot(stdresid ~ PRESEK, col = rainbow(3))
# par(mfrow=c(1,1))
# plot(stdresid ~ Cow, col = rainbow(37))


# BOXCOX transformation
boxcox(AMV2)

# Check for outliers
temp <- sort(stdresid)
tail(temp, 20)                                                                  # 2 obs with std.resid > 2
temp[1:20]                                                                      # 1 obs with std.resid < -2
length(temp) * (2 * pnorm(-2))                                                  # Expect around 2.5 std.resid numerically > 2.  
remove(temp)


# Check influence measures
infl <- influence.measures(AMV2)$infmat


# Normality tests
shapiro.test(stdresid)
lillie.test(stdresid)
cvm.test(stdresid)
ad.test(stdresid)


### Model reduction ------------------------------------------------------------

drop1(AMV2, test = 'F')
AMV2 <- update(AMV2, . ~ . - QRAT)                                              # QRAT is not significant.

drop1(AMV2, test = 'F')
summary(AMV2)


### Interpretation of results --------------------------------------------------

# Question 4:
new_obs <- data.frame(TIME = factor(1), PLANT = "RENO_N", LAB = "KK",
                      O2COR = 0.5, NEFFEKT = -0.01, QRAT = 0.5)

new_pred1 <- predict(AMV2, new_obs, interval = "predict")
new_pred2 <- exp(new_pred1) 
# We get at new prediction of DIOX = 338.55 with the prediction interval 
# (124.31; 922.01).


# Question 5:
summary(AMV2)
# Notice that the operating condition QRAT is insignificant, while the other
# estimates are:
# O2COR = 0.1829
# NEFEKT = 2.9013.

# Remeber however that the response is log-transformed, so we get that
# a unit increase in O2COR yields a relative increase in DIOX of
# exp(0.1829) = 1.201 relative increase in DIOX from one unit increase in O2COR
# exp(2.9013) = 18.198 relative increase in DIOX from one unit increase in NEFFEKT.

# Thus lowering O2COR yiels a small decrease in DIOX and lowering NEFEKT yields
# a drastic decrease in DIOX.


# Question 6:
summary(AMV2)
# Both lab and plant factor are significant, so there is a difference between the
# plants and the labs respectively.
# 
# Looking at the parameters we see that 
# KARA = 6.502
# RENO_N = 5.762
# RENO_S = 4.203
# so we get:
# exp(6.502-5.762) a relative increase in DIOX of 2.096 going from RENO_N to KARA.
# exp(5.762-4.203) a relative increase in DIOX of 4.754 going from RENO_S to RENO_N.
# 
# 
# 
# Looking at the parameters we see that 
# USA = 0.408
# so we get:
# exp(0-0.408) a relative increase in DIOX of 0.665 going from USA to KK.


# Emmeans, Contrasts and CLD.
emmeans(AMV2, ~ PLANT)
emmeans(AMV2, pairwise ~ PLANT, adjust="none")

emmeans(AMV2, ~ LAB)
emmeans(AMV2, pairwise ~ LAB, adjust="none")

em2 <- emmeans(AMV2, pairwise ~ PLANT, adjust="none")
test <- cld(em2[[1]], adjust="bonf", level=0.05, decreasing=TRUE)
# xtable(cld(em2[[1]], adjust="bonf", level=0.05, decreasing=TRUE), 
#        digits=c(1,0,1,1,1,1,1,1))
 
# Boxplot with Compact Letter Display for PLANT
old.par <- par(no.readonly=TRUE)
par(mai=c(1,1,1.25,1))
plot(test, col=2:6)
par(old.par)



### Post Hoc Analysis ----------------------------------------------------------


# . ----------------------------------------------------------------------------
# More advanced model ----------------------------------------------------------


### Imputating missing data ---------------------------------------------------- 

is.na(DATA)    # SO2: Obs 7, 8, 53 and CO: obs 8

SO2_mean <- mean(SO2, na.rm = TRUE)
SO2_sd <- sd(SO2, na.rm = TRUE)
DATA$SO2[7] <- rnorm(1, SO2_mean, SO2_sd)
DATA$SO2[8] <- rnorm(1, SO2_mean, SO2_sd)
DATA$SO2[53] <- rnorm(1, SO2_mean, SO2_sd)

CO_mean <- mean(CO, na.rm = TRUE)
CO_sd <- sd(CO, na.rm = TRUE)
DATA$CO[8] <- rnorm(1, CO_mean, CO_sd)

attach(DATA)


### Specifying initial model ---------------------------------------------------

AM1 <- lm(log(DIOX) ~ poly(O2COR,2)*PLANT + poly(NEFFEKT,2)*PLANT + poly(QRAT,2)*PLANT + 
                      poly(O2COR,2)*LAB + poly(NEFFEKT,2)*LAB + poly(QRAT,2)*LAB +
                      poly(O2COR,2)*TIME + poly(NEFFEKT,2)*TIME + poly(QRAT,2)*TIME +
                      QROEG + TOVN + TROEG + POVN + CO2 + CO + SO2 + HCL + H2O)
  
summary(AM1)


### Model diagnostics ----------------------------------------------------------
stdresid <- rstandard(AM1)

# Model diagnostics plots
par(mfrow=c(2, 2))
plot(AM1, which=1:4) 
par(mfrow=c(1,1))

# BOXCOX transformation
boxcox(AM1)

# Normality plots
par(mfrow=c(1, 2))
hist(stdresid, main="", probability=TRUE, breaks=15)
curve(dnorm, -3, 3, col="red", lwd=2, add=TRUE)
plot(AM1, which=2)
par(mfrow=c(1, 1))

# Check for outliers                                                            # Maybe too many low std.resid.
temp <- sort(stdresid)
tail(temp, 20)                 # 1 obs with std.resid > 2.
temp[1:20]                     # 4 obs with std.resid < -2 
length(temp) * (2 * pnorm(-2)) # Expect around 2.5 std.resid numerically > 2.
remove(temp)


# Check influence measures                                                      # Observation 17 is extremely
infl <- influence.measures(AM1)$infmat                                          # influential but not an outlier.


# Normality tests
shapiro.test(stdresid)
lillie.test(stdresid)
cvm.test(stdresid)
ad.test(stdresid)

# Wally plot
# wallyplot(AM1)


### Model reduction ------------------------------------------------------------

drop1(AM1, test="F")

AM1 <- update(AM1, . ~ . - TROEG)
drop1(AM1, test="F")

AM1 <- update(AM1, . ~ . - poly(QRAT, 2):TIME + poly(QRAT, 1):TIME)
drop1(AM1, test="F")

AM1 <- update(AM1, . ~ . - poly(QRAT, 1):TIME)
drop1(AM1, test="F")

AM1 <- update(AM1, . ~ . - HCL)
drop1(AM1, test="F")

AM1 <- update(AM1, . ~ . - poly(O2COR, 2):LAB + poly(O2COR, 1):LAB)
drop1(AM1, test="F")

AM1 <- update(AM1, . ~ . - poly(O2COR, 1):LAB)
drop1(AM1, test="F")

AM1 <- update(AM1, . ~ . - QROEG)
drop1(AM1, test="F")

AM1 <- update(AM1, . ~ . - poly(O2COR, 2):PLANT + poly(O2COR, 1):PLANT)
drop1(AM1, test="F")

AM1 <- update(AM1, . ~ . - CO2)
drop1(AM1, test="F")

AM1 <- update(AM1, . ~ . - POVN)
drop1(AM1, test="F")

AM1 <- update(AM1, . ~ . - TOVN)
drop1(AM1, test="F")

AM1 <- update(AM1, . ~ . - poly(O2COR, 1):PLANT)
drop1(AM1, test="F")

AM1 <- update(AM1, . ~ . - SO2)
drop1(AM1, test="F")

AM1 <- update(AM1, . ~ . - PLANT:poly(NEFFEKT, 2) + PLANT:poly(NEFFEKT, 1))
drop1(AM1, test="F")

AM1 <- update(AM1, . ~ . - PLANT:poly(QRAT, 2) + PLANT:poly(QRAT, 1))
drop1(AM1, test="F")

AM1 <- update(AM1, . ~ . - PLANT:poly(NEFFEKT, 1))
drop1(AM1, test="F")

AM1 <- update(AM1, . ~ . - poly(O2COR, 2):TIME +  poly(O2COR, 1):TIME)
drop1(AM1, test="F")

AM1 <- update(AM1, . ~ . - poly(O2COR, 1):TIME)
drop1(AM1, test="F")

AM1 <- update(AM1, . ~ . - PLANT:poly(QRAT, 1))
drop1(AM1, test="F")

AM1 <- update(AM1, . ~ . - LAB:poly(QRAT, 2))
drop1(AM1, test="F")

AM1 <- update(AM1, . ~ . - H2O)
drop1(AM1, test="F")

AM1 <- update(AM1, . ~ . - CO)
drop1(AM1, test="F")

AM1 <- update(AM1, . ~ . - poly(QRAT, 2) + poly(QRAT, 1))
drop1(AM1, test="F")

AM1 <- update(AM1, . ~ . - poly(QRAT, 1))
drop1(AM1, test="F")

AM1 <- update(AM1, . ~ . - poly(NEFFEKT, 2):TIME + poly(NEFFEKT, 1):TIME)
drop1(AM1, test="F")


### Post-Hoc analysis ----------------------------------------------------------

# Summary of models
summary(AM1)


# # Parameters with confidence intervals.
# temp <- summary(analysis)$coef
# model_coef <- data.frame(temp[,1:2],lower=temp[,1]-1.96*temp[,2],
#                          upper=temp[,1]+1.96*temp[,2])
# rm(temp)
# model_coef <- round(model_coef,digits=2)
# xtable(model_coef)
# 
# # Estimate of variance parameters
# alpha <- 0.05
# n <- 72
# sigma <- summary(analysis)$sigma
# sigma
# sigma*n/qchisq(1-alpha/2, df = n-2) ; sigma*n/qchisq(alpha/2, df = n-2)
# 
# # Emmeans for parameters
# emmeans(analysis, ~ Var)
# emmeans(analysis, ~ N) 
# emmeans(analysis, ~ Var:N) 
# xtable(emmeans(analysis, ~ Var:N), digits=c(1,1,1,1,0,1,1,1))
# 
# # Contrasts for parameters
# em1 <- emmeans(analysis, pairwise ~ Var, by = "N", adjust="none")
# cld(em1[[1]], adjust="bonf", level=0.05, decreasing=TRUE)
# 
# em2 <- emmeans(analysis, pairwise ~ N, by = "Var", adjust="none")
# cld(em2[[1]], adjust="bonf", level=0.05, decreasing=TRUE)
# xtable(cld(em2[[1]], adjust="bonf", level=0.05, decreasing=TRUE), 
#        digits=c(1,0,1,1,1,1,1,1))
# 
# # Boxplot with Compact Letter Display for ENV
# tuk2 <- glht(analysis, linfct = mcp(Var:N = "Tukey"))
# tuk.cld2 <- cld(tuk2, level=0.05)
# old.par <- par(no.readonly=TRUE)
# par(mai=c(1,1,1.25,1))
# plot(tuk.cld2, col=2:6)
# par(old.par)
# 
# 
# 
# 


# . ----------------------------------------------------------------------------
# Weighing the labs ------------------------------------------------------------


library(nlme)
WL1 <- gls(log(DIOX) ~ poly(O2COR,2)*PLANT + poly(NEFFEKT,2)*PLANT + poly(QRAT,2)*PLANT + 
             poly(O2COR,2)*LAB + poly(NEFFEKT,2)*LAB + poly(QRAT,2)*LAB +
             poly(O2COR,2)*TIME + poly(NEFFEKT,2)*TIME + poly(QRAT,2)*TIME +
             QROEG + TOVN + TROEG + POVN + CO2 + CO + SO2 + HCL + H2O, 
           weights = varIdent(form = ~ 1 | LAB), method = "REML")


WL2 <- gls(log(DIOX) ~ O2COR  + NEFFEKT + PLANT + TIME + LAB,
           weights = varIdent(form = ~ 1 | LAB), method = "ML")

anova(WL2, AMV2)


?logLik

LRT <- -2 * (logLik(AMV2) - logLik(WL2))
p <- 1 - pchisq(LRT, 1)
p
