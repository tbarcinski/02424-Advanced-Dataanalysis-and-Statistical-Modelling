# Preamble #####################################################################

## File Description ############################################################
#   
#   Søren Skjernaa - s223316
#   10/05-2023
#   
#   Advanced Data analysis and Statistical Modelling
#   Assignment 3
#   
#   Notes: 
#         The data file for the analysis is assumed to be found at the relative
#         paths "data/clotingFullAss03.csv".
#
###############################################################################¤

## Clean up ####################################################################
rm(list = ls())
if(!is.null(dev.list())) dev.off()

## Libraries ###################################################################
# Plotting
library(ggplot2)     # Nice plots
library(patchwork)   # Layout of plots
library(gridExtra)

# Handling data frames
library(tidyr)       # Reshape data frames
library(dplyr)       # Rename column in data frames

# Statistics
library(lme4)

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
library(numDeriv)    # Finding Hessians

## Functions ###################################################################


## Data set ####################################################################
dfData <- read.csv("data/clothingFullAss03.csv", header=TRUE, sep=",")
str(dfData)
dfData <- mutate_at(dfData, c("sex", "subjId", "day", "subDay", "time2"), 
                    as.factor)
summary(dfData)
attach(dfData)


## Other #######################################################################


# _ ############################################################################
# Part A: Mixed effects models #################################################


## 1. Plots and summary of data ################################################

# Check for missing data
which(is.na(dfData)) 

# Tabulate data
table(sex)            # not balanced
table(day)            # not balanced
table(time2)          # not balanced
table(subjId)         # not balanced
table(day, time2)     # not balanced
table(subjId, day)    # not balanced
table(subjId, time2)  # not balanced

# Summary of data
tapply(clo, sex, mean)
tapply(clo, sex, sd)

# Boxplot of subjId
p1 <- ggplot(dfData, aes(x=subjId, y=clo, col=sex)) +
  geom_boxplot() +
  labs() +
  theme(plot.title = element_text(hjust=0.5, size=14, face="bold"),
        axis.title = element_text(size=12, face="bold"),
        legend.title = element_text(size=12, face="bold",),
        axis.text.x=element_blank())

# Histograms of clothing insulation level
p2 <- ggplot(dfData, aes(x = clo, fill = sex)) +
  geom_histogram(alpha=0.3, position="identity", bins=20) +
  labs() +
  theme(plot.title = element_text(hjust=0.5, size=14, face="bold"),
        axis.text = element_text(size=12),
        axis.title = element_text(size=12, face="bold"),
        legend.title = element_text(size=12, face="bold"),
        legend.position = "none",
        legend.text = element_text(size=11))

# Plot of clothing insulation levels against temperatures
p3 <- ggplot(dfData, aes(x=tOut, y=clo, col=sex)) +
  geom_point(size = 2) +
  labs(subtitle = "Outdoor temperature") +
  theme(plot.title = element_text(hjust=0.5, size=14, face="bold"),
        axis.title = element_text(size=12, face="bold"),
        legend.title = element_text(size=12, face="bold",),
        legend.position = "none",
        axis.text = element_text(size=12))

p4 <- ggplot(dfData, aes(x=tInOp, y=clo, col=sex)) +
  geom_point(size = 2) +
  labs(subtitle = "Indoor operating temperature") +
  theme(plot.title = element_text(hjust=0.5, size=14, face="bold"),
        axis.title = element_text(size=12, face="bold"),
        legend.title = element_text(size=12, face="bold",),
        legend.position = "none",
        axis.text = element_text(size=12))

p1 + p2 + p3 + p4 + plot_layout(guides = "collect", ncol = 2) +
  plot_annotation(title = "",
                  theme = theme(plot.title = element_text(size = 16,
                                                          face = "bold",
                                                          hjust = 0.5)))

# Boxplot of subDay
ggplot(dfData, aes(x=subDay, y=clo, col=subjId)) +
  geom_boxplot() +
  labs() +
  theme(plot.title = element_text(hjust=0.5, size=14, face="bold"),
        axis.title = element_text(size=12, face="bold"),
        legend.title = element_text(size=12, face="bold",),
        legend.position = "none",
        axis.text.x=element_blank()) + 
  scale_color_manual(values = c(rep(rainbow(6), 7), rainbow(5)))



## 2. subjId model #############################################################

### Definition of models #######################################################

# Model 1 
mSubjId1 <- lmer(sqrt(clo) ~ sex * poly(tInOp, 2) * poly(tOut, 2) + 
                   (1| subjId))
summary(mSubjId1)
AIC(mSubjId1)

# Model 2
mSubjId2 <- lmer(sqrt(clo) ~ sex * poly(tInOp, 2) * poly(tOut, 2) + 
                   (sex + tInOp + tOut| subjId))
summary(mSubjId2)
AIC(mSubjId2)

# Model 3
mSubjId3 <- lmer(sqrt(clo) ~ sex * poly(tInOp, 2) * poly(tOut, 2) + 
                   (sex * tInOp * tOut| subjId))
summary(mSubjId3)
AIC(mSubjId3)

# Model 4
mSubjId4 <- lmer(sqrt(clo) ~ sex * poly(tInOp, 2) * poly(tOut, 2) + 
                   (sex * poly(tInOp, 2) * poly(tOut, 2)| subjId))
summary(mSubjId4)
AIC(mSubjId4)

# Model 5
mSubjId5 <- lmer(sqrt(clo) ~ sex * poly(tInOp, 2) * poly(tOut, 2) + 
                   (sex + tOut| subjId))
summary(mSubjId5)
AIC(mSubjId5)

### Model reduction ############################################################




### Reducing random effects ###

# We find that the interaction of random effects is not necessary
anova(mSubjId4, mSubjId2, refit = FALSE)
LRT <- -2 * (logLik(mSubjId2) - logLik(mSubjId3))
p <- 1 - pchisq(LRT, df = (55 - 29))
c(LRT, p)

# We cannot restrict to the simplest model random effects structure
anova(mSubjId2, mSubjId1, refit = FALSE)
LRT <- -2 * (logLik(mSubjId1) - logLik(mSubjId2))
p <- 1 - pchisq(LRT, df = (29 - 20))
c(LRT, p)

# We cannot restrict the random effects structure by removing tOut
anova(mSubjId2, mSubjId4, refit = FALSE)
LRT <- -2 * (logLik(mSubjId4) - logLik(mSubjId2))
p <- 1 - pchisq(LRT, df = (29 - 25))
c(LRT, p)

# We cannot restrict the random effects structure by removing tInOp
anova(mSubjId2, mSubjId5, refit = FALSE)
LRT <- -2 * (logLik(mSubjId5) - logLik(mSubjId2))
p <- 1 - pchisq(LRT, df = (29 - 20))
c(LRT, p)


### Reducing fixed effects ###

# Remove the three-way interaction
mSubjId2 <- lmer(sqrt(clo) ~ sex * poly(tInOp, 2) * poly(tOut, 2) + 
                   (sex + tInOp + tOut| subjId), REML = FALSE)
summary(mSubjId2)
drop1(mSubjId2, test = "Chisq")
mSub <- update(mSubjId2, . ~ . - sex:poly(tInOp, 2):poly(tOut, 2))
summary(mSub)
AIC(mSub)


# Remove the sex:tInOp interaction
drop1(mSub, test = "Chisq")
mSub <- update(mSub, . ~ . - sex:poly(tInOp, 2))
summary(mSub)
AIC(mSub)


# Further reduction are rejected
drop1(mSub, test = "Chisq")
mSub2 <- lmer(sqrt(clo) ~ sex + poly(tInOp, 2) + poly(tOut, 2) +
                         sex:poly(tOut, 2) +
                         poly(tInOp, 1):poly(tOut, 1) + 
                         (sex + tInOp + tOut| subjId), REML = FALSE)
LRT <- -2 * (logLik(mSub2) - logLik(mSub))
p <- 1 - pchisq(LRT, df = (23 - 21))
c(LRT, p)


### Final model ################################################################
mSubjIdFinal <- update(mSub, REML = TRUE)
summary(mSubjIdFinal)
AIC(mSubjIdFinal)


### Model diagnostics ##########################################################

par(mfrow = c(3,2))
residuals <- resid(mSubjIdFinal)
# infl <- influence(mSubjIdFinal, obs = TRUE)
# leverage <- infl[, ncol(infl)]
cutoff <- 2 * 23 / length(residuals)

p1 <- plot(fitted(mSubjIdFinal), residuals, col = sex, 
           main = "Gaussian with transformation sqrt(clo)",
           xlab = "Fitted values", ylab = "Residuals")
p1 <- p1 + abline(0,0)

qqnorm(residuals, col = sex,
       main = "", xlab = "Theoretical quantiles of N(0,1)-distribution", 
       ylab = "Sample Quantiles of Resid.")

plot(sex, residuals, xlab = "sex", ylab = "clo")

plot(tInOp, residuals, col = sex, xlab = "tInOp", ylab = "clo")
abline(0,0)

plot(tOut, residuals, col = sex, xlab = "tOut", ylab = "clo")
abline(0,0)

# plot(leverage, studDevResid, col = sex, 
#      xlab = "Leverage", ylab = "Stud. Resid")
# abline(h = -2, col = "blue", lwd=2, lty=2)
# abline(h =  2, col = "blue", lwd=2, lty=2)
# abline(v = cutoff, col = "blue", lwd=2, lty=2)
# sum(studDevResid > 2 | studDevResid < -2)


par(mfrow = c(2,2))
qqnorm(ranef(mSubjIdFinal)$subjId$'(Intercept)')
qqnorm(ranef(mSubjIdFinal)$subjId$sexmale)
qqnorm(ranef(mSubjIdFinal)$subjId$tInOp)
qqnorm(ranef(mSubjIdFinal)$subjId$tOut)
par(mfrow = c(1,1))


## 3. subjId + day model #######################################################

### Definition of models #######################################################

# Model 1 
mSubjDay1 <- lmer(sqrt(clo) ~ sex * poly(tInOp, 2) * poly(tOut, 2) + 
                   (1| subjId/subDay))
summary(mSubjDay1)
AIC(mSubjDay1)

# Model 2
mSubjDay2 <- lmer(sqrt(clo) ~ sex * poly(tInOp, 2) * poly(tOut, 2) + 
                   (sex + tInOp + tOut| subjId/subDay))
summary(mSubjDay2)
AIC(mSubjDay2)


# Model 4
mSubjDay4 <- lmer(sqrt(clo) ~ sex * poly(tInOp, 2) * poly(tOut, 2) + 
                   (sex + tInOp| subjId/subDay))
summary(mSubjDay4)
AIC(mSubjDay4)

# Model 5
mSubjDay5 <- lmer(sqrt(clo) ~ sex * poly(tInOp, 2) * poly(tOut, 2) + 
                   (sex + tOut| subjId/subDay))
summary(mSubjDay5)
AIC(mSubjDay5)

### Model reduction ############################################################

### Reducing random effects ###

# We cannot restrict to the simplest model random effects structure
anova(mSubjDay2, mSubjDay1, refit = FALSE)
# LRT <- -2 * (logLik(mSubjDay1) - logLik(mSubjDay2))
# p <- 1 - pchisq(LRT, df = (39 - 21))
# c(LRT, p)


# We cannot restrict the random effects structure by removing tOut
anova(mSubjDay2, mSubjDay4, refit = FALSE)
# LRT <- -2 * (logLik(mSubjDay4) - logLik(mSubjDay2))
# p <- 1 - pchisq(LRT, df = (29 - 25))
# c(LRT, p)

# We cannot restrict the random effects structure by removing tInOp
anova(mSubjDay2, mSubjDay5, refit = FALSE)
# LRT <- -2 * (logLik(mSubjDay5) - logLik(mSubjDay2))
# p <- 1 - pchisq(LRT, df = (29 - 20))
# c(LRT, p)


### Reducing fixed effects ###

# Remove the second order polynomial three way interaction
mSubjDay2 <- lmer(sqrt(clo) ~ sex * poly(tInOp, 2) * poly(tOut, 2) + 
                   (sex + tInOp + tOut| subjId/subDay), REML = FALSE)
summary(mSubjDay2)
drop1(mSubjDay2, test = "Chisq")
mSub <- lmer(sqrt(clo) ~ sex + poly(tInOp, 2) + poly(tOut, 2) + 
               sex:poly(tInOp, 2) + sex:poly(tOut, 2) + 
               poly(tInOp, 1):poly(tOut, 1) + 
               sex:poly(tInOp, 1):poly(tOut, 1) + 
               (sex + tInOp + tOut| subjId/subDay), REML = FALSE)
anova(mSubjDay2, mSub)
summary(mSub)
AIC(mSub)
# drop1(mSub, test = "Chisq")



### Final model ################################################################
mSubjDayFinal <- update(mSub, REML = TRUE)
summary(mSubjDayFinal)
AIC(mSubjDayFinal)

### Model diagnostics ##########################################################

par(mfrow = c(3,2))
residuals <- resid(mSubjDayFinal)
# infl <- influence(mSubjDayFinal, obs = TRUE)
# leverage <- infl[, ncol(infl)]
# cutoff <- 2 * 23 / length(residuals)

p1 <- plot(fitted(mSubjDayFinal), residuals, col = sex, 
           main = "Gaussian with transformation sqrt(clo)",
           xlab = "Fitted values", ylab = "Residuals")
p1 <- p1 + abline(0,0)

qqnorm(residuals, col = sex,
       main = "", xlab = "Theoretical quantiles of N(0,1)-distribution", 
       ylab = "Sample Quantiles of Resid.")

plot(sex, residuals, xlab = "sex", ylab = "clo")

plot(tInOp, residuals, col = sex, xlab = "tInOp", ylab = "clo")
abline(0,0)

plot(tOut, residuals, col = sex, xlab = "tOut", ylab = "clo")
abline(0,0)

# plot(leverage, studDevResid, col = sex, 
#      xlab = "Leverage", ylab = "Stud. Resid")
# abline(h = -2, col = "blue", lwd=2, lty=2)
# abline(h =  2, col = "blue", lwd=2, lty=2)
# abline(v = cutoff, col = "blue", lwd=2, lty=2)
# sum(studDevResid > 2 | studDevResid < -2)


par(mfrow = c(4,2))
qqnorm(ranef(mSubjDayFinal)$subjId$'(Intercept)')
qqnorm(ranef(mSubjDayFinal)$'subDay:subjId'$'(Intercept)')
qqnorm(ranef(mSubjDayFinal)$subjId$sexmale)
qqnorm(ranef(mSubjDayFinal)$'subDay:subjId'$sexmale)
qqnorm(ranef(mSubjDayFinal)$subjId$tInOp)
qqnorm(ranef(mSubjDayFinal)$'subDay:subjId'$tInOp)
qqnorm(ranef(mSubjDayFinal)$subjId$tOut)
qqnorm(ranef(mSubjDayFinal)$'subDay:subjId'$tOut)
par(mfrow = c(1,1))


## 4. RMM subDay ###############################################################

### Definition of models #######################################################
library(nlme)

# Model Compound Symmetry
mNo <- lme(sqrt(clo) ~ sex * poly(tInOp, 2) * poly(tOut, 2),
           random = ~ 1| subDay)
summary(mNo)
c(AIC(mNo), BIC(mNo))

# Model Exp
mExp1 <- lme(sqrt(clo) ~ sex * poly(tInOp, 2) * poly(tOut, 2),
               random = ~ 1 | subDay,
               correlation = corExp(form = ~ time | subDay, nugget = TRUE))
summary(mExp1)
plot(Variogram(mExp1, form = ~ time | subDay))
c(AIC(mExp1), BIC(mExp1))

# Model Gaussian
mGauss1 <- lme(sqrt(clo) ~ sex * poly(tInOp, 2) * poly(tOut, 2),
               random = ~ 1 | subDay,
               correlation = corGaus(form = ~ time | subDay, nugget = TRUE))
summary(mGauss1)
plot(Variogram(mGauss1, form = ~ time | subDay))
c(AIC(mGauss1), BIC(mGauss1))

# Model AR1
mAR1 <- lme(sqrt(clo) ~ sex * poly(tInOp, 2) * poly(tOut, 2),
               random = ~ 1 | subDay,
               correlation = corAR1(form = ~ as.integer(time2) | subDay))
summary(mAR1)
plot(Variogram(mAR1, form = ~ as.integer(time2) | subDay))
c(AIC(mAR1), BIC(mAR1))


### Model reduction ############################################################

### Reducing the random effects structure ###
anova(mAR1, mNo)

### Reducing fixed effects ###
mAR1 <- lme(sqrt(clo) ~ sex * poly(tInOp, 2) * poly(tOut, 2),
             random = ~ 1 | subDay,
             correlation = corAR1(form = ~ as.integer(time2) | subDay),
             method = "ML")

# Remove the three way interaction
drop1(mAR1, test = "Chisq")
mSub <- update(mAR1, . ~ . - sex:poly(tInOp, 2):poly(tOut, 2))
summary(mSub)
c(AIC(mSub), BIC(mSub))

# Remove temperature X temperature interaction
drop1(mSub, test = "Chisq")
mSub <- update(mSub, . ~ . - poly(tInOp, 2):poly(tOut, 2))
summary(mSub)
c(AIC(mSub), BIC(mSub))

# Remove sex X tOut interaction
drop1(mSub, test = "Chisq")
mSub <- update(mSub, . ~ . - sex:poly(tOut, 2))
summary(mSub)
c(AIC(mSub), BIC(mSub))

# Remove tOut
drop1(mSub, test = "Chisq")
mSub <- update(mSub, . ~ . - poly(tOut, 2))
summary(mSub)
c(AIC(mSub), BIC(mSub))

# Remove second order polynomial of tInp and its sex interaction.
drop1(mSub, test = "Chisq")
mSub2 <- lme(sqrt(clo) ~ sex*poly(tInOp, 1),
             random = ~ 1 | subDay,
             correlation = corAR1(form = ~ as.integer(time2) | subDay),
             method = "ML")
anova(mSub, mSub2)

# No further reductions are possible
mSub3 <- lme(sqrt(clo) ~ sex*tInOp,
             random = ~ 1 | subDay,
             correlation = corAR1(form = ~ as.integer(time2) | subDay),
             method = "ML")
drop1(mSub3, test = "Chisq")

### Final model ################################################################
mExpFinal <- lme(sqrt(clo) ~ sex*tInOp,
                 random = ~ 1 | subDay,
                 correlation = corAR1(form = ~ as.integer(time2) | subDay),
                 method = "REML")
summary(mExpFinal)
c(AIC(mExpFinal), BIC(mExpFinal))

### Model diagnostics ##########################################################

par(mfrow = c(2,2))
residuals <- resid(mExpFinal)
# infl <- influence(mSubjDayFinal, obs = TRUE)
# leverage <- infl[, ncol(infl)]
# cutoff <- 2 * 23 / length(residuals)

p1 <- plot(fitted(mExpFinal), residuals, col = sex, 
           main = "Gaussian with transformation sqrt(clo)",
           xlab = "Fitted values", ylab = "Residuals")
p1 <- p1 + abline(0,0)

qqnorm(residuals, col = sex,
       main = "", xlab = "Theoretical quantiles of N(0,1)-distribution", 
       ylab = "Sample Quantiles of Resid.")

plot(sex, residuals, xlab = "sex", ylab = "clo")

plot(tInOp, residuals, col = sex, xlab = "tInOp", ylab = "clo")
abline(0,0)

# plot(leverage, studDevResid, col = sex, 
#      xlab = "Leverage", ylab = "Stud. Resid")
# abline(h = -2, col = "blue", lwd=2, lty=2)
# abline(h =  2, col = "blue", lwd=2, lty=2)
# abline(v = cutoff, col = "blue", lwd=2, lty=2)
# sum(studDevResid > 2 | studDevResid < -2)


par(mfrow = c(1,1))
qqnorm(ranef(mExpFinal)$'(Intercept)')
par(mfrow = c(1,1))



## 5. Presentation of final model ##############################################




## 6. Conclusion ###############################################################










# _ ############################################################################
# Part B: Hierachichal modes ###################################################
library(lme4)

## 1. Model 1 ##################################################################

nll <- function(theta){
  X <- matrix(c(rep(1, 803), 
              as.integer(sex == "male")), 
              byrow = FALSE, ncol = 2)
  
}