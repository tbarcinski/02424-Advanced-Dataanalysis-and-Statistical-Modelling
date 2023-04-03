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
rm(list = ls())
if(!is.null(dev.list())) dev.off()


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
library(Gammareg)

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


# _ ############################################################################
# Part A: Clothing level #######################################################


### Load the data and get initial overview #####################################

# Load data
df_data <- read.csv("data/clothing.csv", header = TRUE) 
str(df_data)
df_data$subjDay <- 2^df_data$subjId * 3^df_data$day
df_data <- mutate_at(df_data, c("sex", "subjId", "day", "subjDay"), as.factor)
summary(df_data)
attach(df_data)

# Tabulate data
table(sex)
table(day, sex)

# Summary statistics of data
tapply(clo, sex, mean)
tapply(clo, sex, sd)

# Check the average standard deviation for each combination subjId and Day
summary(clo)                      # 0.23 - 0.97
sd(clo)                           #      0.16008
mean(tapply(clo, subjDay, sd))    # Only 0.02635
mean(tapply(tInOp, subjDay, sd))
mean(tapply(tOut, subjDay, sd))

### Plot of the data ###########################################################

# Boxplot of subjId
p1 <- ggplot(df_data, aes(x=subjId, y=clo, col=sex)) +
  geom_boxplot() +
  labs() +
  theme(plot.title = element_text(hjust=0.5, size=14, face="bold"),
        axis.title = element_text(size=12, face="bold"),
        legend.title = element_text(size=12, face="bold",),
        axis.text.x=element_blank())

# Histograms of clothing insulation level
p2 <- ggplot(df_data, aes(x = clo, fill = sex)) +
  geom_histogram(alpha=0.3, position="identity", bins=20) +
  labs() +
  theme(plot.title = element_text(hjust=0.5, size=14, face="bold"),
        axis.text = element_text(size=12),
        axis.title = element_text(size=12, face="bold"),
        legend.title = element_text(size=12, face="bold"),
        legend.position = "none",
        legend.text = element_text(size=11))


# Plot of clothing insulation levels against temperatures
# df_data = df_data[df_data$subjId %in% c(11,17,35,49),]
p3 <- ggplot(df_data, aes(x=tOut, y=clo, col=sex)) +
  geom_point(size = 1) +
  labs(subtitle = "Outdoor temperature") +
  theme(plot.title = element_text(hjust=0.5, size=14, face="bold"),
        axis.title = element_text(size=12, face="bold"),
        legend.title = element_text(size=12, face="bold",),
        legend.position = "none",
        axis.text = element_text(size=12))

p4 <- ggplot(df_data, aes(x=tInOp, y=clo, col=sex)) +
  geom_point(size = 1) +
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




### Sex models #################################################################


#### Fitting the four models ###################################################

# Gaussian without transformation of clo
lmId <- lm(clo ~ sex * poly(tInOp, 2) * poly(tOut, 2)) 
summary(lmId)
#boxcox(lmId)


# Gaussian with sqrt transformation
lmSqrt <- lm(sqrt(clo) ~ sex * poly(tInOp, 2) * poly(tOut, 2)) 
summary(lmSqrt)


# Gamma with inverse link
glmGammaInv <- glm(clo ~ sex * poly(tInOp, 2) * poly(tOut, 2),           
                   family = Gamma(link = "inverse"))
summary(glmGammaInv)


# Gamma with logarithmic link
glmGammaLog <- glm(clo ~ sex * poly(tInOp, 2) * poly(tOut, 2),           
                   family = Gamma(link = "log"))
summary(glmGammaLog)



#### Model diagnostics #########################################################

# Gaussian without transformation of clo
par(mfrow = c(3,2))
stdDevResid <- rstandard(lmId)
studDevResid <- studres(lmId)
infl <- influence.measures(lmId)$infmat
leverage <- infl[, ncol(infl)]
cutoff <- 2 * length(coef(lmId)) / length(studDevResid)

p1 <- plot(fitted(lmId), stdDevResid, col = sex, 
           main = "Gaussian without transformation of clo",
           xlab = "Fitted values", ylab = "Residuals")
p1 <- p1 + abline(0,0)

qqnorm(stdDevResid, col = sex,
       main = "", xlab = "Theoretical quantiles of N(0,1)", 
       ylab = "Sample Quantiles of Resid.")
abline(0,1)

plot(sex, stdDevResid, xlab = "sex", ylab = "clo")

plot(tInOp, stdDevResid, col = sex, xlab = "tInOp", ylab = "clo")
abline(0,0)

plot(tOut, stdDevResid, col = sex, xlab = "tOut", ylab = "clo")
abline(0,0)

plot(leverage, studDevResid, col = sex, 
     xlab = "Leverage", ylab = "Stud. Resid")
abline(h = -2, col = "blue", lwd=2, lty=2)
abline(h =  2, col = "blue", lwd=2, lty=2)
abline(v = cutoff, col = "blue", lwd=2, lty=2)
sum(studDevResid > 2 | studDevResid < -2)


# Gaussian with sqrt transformation
par(mfrow = c(3,2))
stdDevResid <- rstandard(lmSqrt)
studDevResid <- studres(lmSqrt)
infl <- influence.measures(lmSqrt)$infmat
leverage <- infl[, ncol(infl)]
cutoff <- 2 * length(coef(lmSqrt)) / length(studDevResid)

p1 <- plot(fitted(lmSqrt), stdDevResid, col = sex, 
           main = "Gaussian with transformation sqrt(clo)",
           xlab = "Fitted values", ylab = "Residuals")
p1 <- p1 + abline(0,0)

qqnorm(stdDevResid, col = sex,
       main = "", xlab = "Theoretical quantiles of N(0,1)", 
       ylab = "Sample Quantiles of Resid.")
abline(0,1)

plot(sex, stdDevResid, xlab = "sex", ylab = "clo")

plot(tInOp, stdDevResid, col = sex, xlab = "tInOp", ylab = "clo")
abline(0,0)

plot(tOut, stdDevResid, col = sex, xlab = "tOut", ylab = "clo")
abline(0,0)

plot(leverage, studDevResid, col = sex, 
     xlab = "Leverage", ylab = "Stud. Resid")
abline(h = -2, col = "blue", lwd=2, lty=2)
abline(h =  2, col = "blue", lwd=2, lty=2)
abline(v = cutoff, col = "blue", lwd=2, lty=2)
sum(studDevResid > 2 | studDevResid < -2)



# Gamma with inverse link
par(mfrow = c(3,2))
stdDevResid <- rstandard(glmGammaInv, type = "deviance")
studDevResid <- rstudent(glmGammaInv, type = "deviance")
infl <- influence.measures(glmGammaInv)$infmat
leverage <- infl[, ncol(infl)]
cutoff <- 2 * length(coef(glmGammaInv)) / length(studDevResid)

p1 <- plot(fitted(glmGammaInv), stdDevResid, col = sex, 
           main = "Gamma with inverse link",
           xlab = "Fitted values", ylab = "Residuals")
p1 <- p1 + abline(0,0)

qqnorm(stdDevResid, col = sex,
       main = "", xlab = "Theoretical quantiles of N(0,1)", 
       ylab = "Sample Quantiles of Resid.")
abline(0,1)

plot(sex, stdDevResid, xlab = "sex", ylab = "clo")

plot(tInOp, stdDevResid, col = sex, xlab = "tInOp", ylab = "clo")
abline(0,0)

plot(tOut, stdDevResid, col = sex, xlab = "tOut", ylab = "clo")
abline(0,0)

plot(leverage, studDevResid, col = sex, 
     xlab = "Leverage", ylab = "Stud. Resid")
abline(h = -2, col = "blue", lwd=2, lty=2)
abline(h =  2, col = "blue", lwd=2, lty=2)
abline(v = cutoff, col = "blue", lwd=2, lty=2)
sum(studDevResid > 2 | studDevResid < -2)


# Gamma with logarithmic link
par(mfrow = c(3,2))
stdDevResid <- rstandard(glmGammaLog, type = "deviance")
studDevResid <- rstudent(glmGammaLog, type = "deviance")
infl <- influence.measures(glmGammaLog)$infmat
leverage <- infl[, ncol(infl)]
cutoff <- 2 * length(coef(glmGammaLog)) / length(studDevResid)

p1 <- plot(fitted(glmGammaLog), stdDevResid, col = sex, 
           main = "Gamma with logarithmic link",
           xlab = "Fitted values", ylab = "Residuals")
p1 <- p1 + abline(0,0)

qqnorm(stdDevResid, col = sex,
       main = "", xlab = "Theoretical quantiles of N(0,1)", 
       ylab = "Sample Quantiles of Resid.")
abline(0,1)

plot(sex, stdDevResid, xlab = "sex", ylab = "clo")

plot(tInOp, stdDevResid, col = sex, xlab = "tInOp", ylab = "clo")
abline(0,0)

plot(tOut, stdDevResid, col = sex, xlab = "tOut", ylab = "clo")
abline(0,0)

plot(leverage, studDevResid, col = sex, 
     xlab = "Leverage", ylab = "Stud. Resid")
abline(h = -2, col = "blue", lwd=2, lty=2)
abline(h =  2, col = "blue", lwd=2, lty=2)
abline(v = cutoff, col = "blue", lwd=2, lty=2)
sum(studDevResid > 2 | studDevResid < -2)


#### Model reduction ###########################################################

# Comparison of the models
AIC(lmId)
AIC(lmSqrt)
AIC(glmGammaInv)
AIC(glmGammaLog)

# Adjust AIC for square root transformation
sqrt_derivative <- function(y_input){
  return(1/(2*sqrt(y_input)))
}
to_subtract <- sum(log(sqrt_derivative(clo)))
AIC(lmSqrt) - 2*to_subtract


# Model reduction
summary(glmGammaInv)

glmGammaInv1 <- glmGammaInv
glmGammaInv2 <- glm(clo ~ sex + poly(tInOp, 2) + poly(tOut, 2) + 
                      sex:poly(tInOp, 2) + sex:poly(tOut, 2) + 
                      poly(tInOp, 1):poly(tOut, 1) + 
                      sex:poly(tInOp, 1):poly(tOut, 1),
                    family = Gamma(link = "inverse"))

anova(glmGammaInv2, glmGammaInv1, test = "LRT")
drop1(glmGammaInv2, test = "LRT")
xtable(drop1(glmGammaInv2, test = "LRT")[, c(1,2,3,4)])



#### Post-Hoc Analysis #########################################################

# Model parameters
parameters <- cbind(summary(glmGammaInv2)$coefficients[, 1:2],
                    confint(glmGammaInv2))
xtable(round(parameters, 2))
summary(glmGammaInv2)$dispersion

# Quantiles of temperatures
tOutLevels <- as.numeric(summary(tOut)[c(2, 3, 5)])
tInOpLevels <- as.numeric(summary(tInOp)[c(2, 3, 5)])

# Model contrasts
observedEmmeans <- emmeans(glmGammaInv2, "tOut", by = c("sex", "tInOp"),
                           at = list(tOut = round(tOutLevels,0),
                                     tInOp = round(tInOpLevels,0)),
                           type = "response")
plot(observedEmmeans)


# Graphical presentation - tInOp
tInOpSeq <- seq(min(tInOp), max(tInOp), 0.1)

df_tInOp <- 
  data.frame(sex = c(rep("female", length(tInOpSeq) * length(tOutLevels)), 
                     rep("male", length(tInOpSeq) * length(tOutLevels))),
             tInOp = rep(tInOpSeq, 2 * length(tOutLevels)),
             tOut = rep(c(rep(tOutLevels[1], length(tInOpSeq)),
                          rep(tOutLevels[2], length(tInOpSeq)),
                          rep(tOutLevels[3], length(tInOpSeq))), 2))
df_tInOp <- cbind(df_tInOp, 
                  "clo" = predict(glmGammaInv2, df_tInOp, type = "response"),
                  "group" = 2^(df_tInOp$tOut) * 
                    3^as.numeric(as.factor(df_tInOp$sex)))

ggplot(df_data, aes(x=tInOp, y=clo, col=sex)) +
  geom_point(size = 1) +
  geom_line(data = df_tInOp, aes(x=tInOp, y=clo, group=group), linewidth = 1)
labs(subtitle = "Outdoor temperature") +
  theme(plot.title = element_text(hjust=0.5, size=14, face="bold"),
        axis.title = element_text(size=12, face="bold"),
        legend.title = element_text(size=12, face="bold",),
        axis.text = element_text(size=12, face="bold"))


# Graphical presentation - tOut
tOutSeq <- seq(min(tOut), max(tOut), 0.1)

df_tOut <- 
  data.frame(sex = c(rep("female", length(tOutSeq) * length(tInOpLevels)), 
                     rep("male", length(tOutSeq) * length(tInOpLevels))),
             tOut = rep(tOutSeq, 2 * length(tInOpLevels)),
             tInOp = rep(c(rep(tInOpLevels[1], length(tOutSeq)),
                           rep(tInOpLevels[2], length(tOutSeq)),
                           rep(tInOpLevels[3], length(tOutSeq))), 2))
df_tOut <- cbind(df_tOut, 
                 "clo" = predict(glmGammaInv2, df_tOut, type = "response"), 
                 "group" = 2^(df_tOut$tInOp) * 
                   3^as.numeric(as.factor(df_tOut$sex)))

ggplot(df_data, aes(x=tOut, y=clo, col=sex)) +
  geom_point(size = 1) +
  geom_line(data = df_tOut, aes(x=tOut, y=clo, group = group), linewidth = 1)
labs(subtitle = "Outdoor temperature") +
  theme(plot.title = element_text(hjust=0.5, size=14, face="bold"),
        axis.title = element_text(size=12, face="bold"),
        legend.title = element_text(size=12, face="bold",),
        axis.text = element_text(size=12, face="bold"))



### Subject model ##############################################################

#### Fit model #################################################################  

glmGammaSubj <- glm(clo ~ subjId * poly(tInOp, 2) * poly(tOut, 2),           
                    family = Gamma(link = "inverse"))
summary(glmGammaSubj)


#### Model Diagnostics #########################################################

par(mfrow = c(3,2))
stdDevResid <- rstandard(glmGammaSubj, type = "deviance")
studDevResid <- rstudent(glmGammaSubj, type = "deviance")
infl <- influence.measures(glmGammaSubj)$infmat
leverage <- infl[, ncol(infl)]
cutoff <- 2 * length(coef(glmGammaSubj)) / length(studDevResid)

p1 <- plot(fitted(glmGammaSubj), stdDevResid, col = sex, 
           main = "Gamma with inverse link",
           xlab = "Fitted values", ylab = "Residuals")
p1 <- p1 + abline(0,0)

qqnorm(stdDevResid, col = sex,
       main = "", xlab = "Theoretical quantiles of N(0,1)", 
       ylab = "Sample Quantiles of Resid.")
abline(0,1)

plot(sex, stdDevResid, xlab = "sex", ylab = "clo")

plot(tInOp, stdDevResid, col = sex, xlab = "tInOp", ylab = "clo")
abline(0,0)

plot(tOut, stdDevResid, col = sex, xlab = "tOut", ylab = "clo")
abline(0,0)

plot(leverage, studDevResid, col = sex, 
     xlab = "Leverage", ylab = "Stud. Resid")
abline(h = -2, col = "blue", lwd=2, lty=2)
abline(h =  2, col = "blue", lwd=2, lty=2)
abline(v = cutoff, col = "blue", lwd=2, lty=2)
sum(studDevResid > 2 | studDevResid < -2)
nrow(summary(glmGammaSubj)$coefficients)


#### Model reduction ###########################################################

# Model reduction
summary(glmGammaSubj)

glmGammaSubj1 <- glmGammaSubj
glmGammaSubj2 <- glm(clo ~ subjId + poly(tInOp, 2) + poly(tOut, 2) + 
                       subjId:poly(tInOp, 2) + subjId:poly(tOut, 2) + 
                       poly(tInOp, 1):poly(tOut, 1) + 
                       subjId:poly(tInOp, 1):poly(tOut, 1),
                     family = Gamma(link = "inverse"))
anova(glmGammaSubj2, glmGammaSubj1, test = "LRT")

drop1(glmGammaSubj, test = "LRT")
xtable(round(drop1(glmGammaSubj, test = "LRT"), 2))

# stepAIC(glmGammaSubj, direction = "backward")


#### Final model ###############################################################

summary(lmSqrt1)
AIC(lmSqrt1)


### Within subj:day autocorrelation ############################################

acf1SubjDay <- numeric(length(unique(subjDay)))
j <- 1
tempDf1 <- cbind(df_data, stdDevResid)

for (i in unique(subjDay)){
  tempDf2 <- filter(tempDf1, subjDay == i)
  if (length(unique(tempDf2$clo)) == 1){
    acf1SubjDay[j] <- 1
  } else{
    acf1SubjDay[j] <- acf(tempDf2$stdDevResid, lag.max = 1, plot = FALSE)$acf[2]
  }
  j <- j + 1
}

par(mfrow = c(1,1))
plot(acf1SubjDay, 
     xlab = "Index for Subject:Day combination",
     ylab = "ACF (lag = 1)")
abline(0,0)
mean(acf1SubjDay)


### Different dispersions ######################################################

#### Optimizing dispersion parameters ##########################################

# Model
glmGammaInv2 <- glm(clo ~ sex + poly(tInOp, 2) + poly(tOut, 2) + 
                      sex:poly(tInOp, 2) + sex:poly(tOut, 2) + 
                      poly(tInOp, 1):poly(tOut, 1) + 
                      sex:poly(tInOp, 1):poly(tOut, 1),
                    family = Gamma(link = "inverse"))


# Negative log likelihood
objective <- function(beta, formula){
  
  Xmu <- model.matrix(formula)
  Xdispersion <- cbind(rep(1, 803), as.integer(df_data$sex == "female"))
  
  eta <- Xmu %*% beta[1:ncol(Xmu)]
  mu <- 1 / eta
  dispersion <- exp(Xdispersion %*% beta[c(ncol(Xmu) + 1, ncol(Xmu) + 2)])
  
  alpha <- 1 / dispersion
  beta <- mu * dispersion
  
  return(-sum(dgamma(clo, shape = alpha, scale = beta, log = TRUE))) 
}

# Model formula
linModel <- clo ~ sex + poly(tInOp, 2) + poly(tOut, 2) +
  sex:poly(tInOp, 2) + sex:poly(tOut, 2) +
  poly(tInOp, 1):poly(tOut, 1) +
  sex:poly(tInOp, 1):poly(tOut, 1)
# linModel <- clo ~ tInOp + tOut

# Initial parameters
coef <- coefficients(glmGammaInv2)
initialBeta <- c(coef, "Overall_dispersion" = 1, "female_weight" = 1)
# initialBeta <- c(rep(1, 5))

# Optimization
opt <- nlminb(initialBeta, objective, formula = linModel,
              control = list(eval.max = 1000, iter.max = 1000))
opt


# Control of results
stdDevResid <- rstandard(glmGammaInv2, type = "deviance")
sd(stdDevResid[sex == "female"])^2 / sd(stdDevResid[sex == "male"])^2
exp(opt$par[length(opt$par)])

### Gammareg ###
# formula.mean = linModel
# Z1 <- as.integer(df_data$sex == "female")
# formula.shape = ~ Z1
# a=Gammareg(formula.mean, formula.shape, meanlink="log")
# summary(a)
# opt$par


#### Profile likelihood ########################################################

# Optimize zeta, for each value fixed tau (see page 34)
profile_objective <- function(tau, formula){
  
  inner_objective <- function(zeta, tau_input){
    objective(c(zeta, tau), formula)
  }
  
  zeta_length <- ncol(model.matrix(formula)) + 1
  initial_zeta <- c(rep(1, zeta_length))
  nlminb(initial_zeta, inner_objective, tau_input = tau)$objective
}

# Calculate profile log likelihoods
tau <- seq(opt$par[14] - 0.5, opt$par[14] + 0.5, 0.01)
profile_logLikelihood <- numeric(length(tau))
for (i in 1:length(tau)){
  profile_logLikelihood[i] <- profile_objective(tau[i], formula = linModel)
  print(c(i, tau[i], profile_logLikelihood[i]))
}

# Finding the profile likelihood confidence interval (see page 36)
profile_likelihood <- exp(-profile_logLikelihood) 
profile_likelihood <- profile_likelihood / max(profile_likelihood)
# profile_likelihood <- exp(profile_likelihood)
L_CI_lower <- min(tau[profile_likelihood > exp(-(1 / 2) *  qchisq(0.95, df=1))])
L_CI_upper <- max(tau[profile_likelihood > exp(-(1 / 2) *  qchisq(0.95, df=1))])


# Finding the quadratic approximation (Wald Confidence intervals) (see page 23)
# observed_hessian <- hessian(objective, opt$par, formula = linModel)
# observed_hessian_tau <- - observed_hessian[14, 14]
observed_hessian_tau <- hessian(profile_objective, opt$par[14], formula = linModel)
quadratic_approx_logLik <- - opt$objective + 
  (1 / 2) * -observed_hessian_tau * (tau - opt$par[14])^2
quadratic_approx_Lik <- exp(quadratic_approx_logLik)
quadratic_approx_Lik <- quadratic_approx_Lik / max(quadratic_approx_Lik)
# plot(tau, quadratic_approx_Lik)
# standard_error = sqrt(diag(solve(observed_hessian)))[14]
standard_error = sqrt(diag(solve(observed_hessian)))[1]


# Plot of profile likelihood and quadratic approximation. 
par(mfrow = c(1,1))
plot(tau, profile_likelihood, type = "l",
     xlab="Female Dispersion Weight", ylab="Profile Likelihood",
     main="The comparison between Wald's and Likelihood based Confidence Intervals")
lines(tau, rep(exp(-(1 / 2) *  qchisq(0.95, df=1)), length(tau)), col = 2)
rug(L_CI_lower, ticksize = 0.2, lwd = 2, col = "red")
rug(L_CI_upper, ticksize = 0.2, lwd = 2, col = "red")
c(L_CI_lower, L_CI_upper)
abline(v = opt$par["female_weight"], lty = 2)
lines(tau, quadratic_approx_Lik, col = "blue")
rug(opt$par["female_weight"] - qnorm(0.975) * standard_error,
    ticksize = 0.2, lwd = 2, col = "blue")
rug(opt$par["female_weight"] + qnorm(0.975) * standard_error,
    ticksize = 0.2, lwd = 2, col = "blue")
opt$par["female_weight"] - qnorm(0.975) * standard_error * c(1, -1)
legend("topright", 95,
       legend=c("Profile likelihood", "Quadratic approximation",
                "95% confidence interval"),
       col=c("black", "blue", "red"), lty = 1:1, cex=0.8,
       inset = 0.02)