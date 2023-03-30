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

# Full model
mf <- glm(infections ~ swimmer*age*sex*location, 
          family = poisson(link = "log"), offset = log(persons),
          data = df_data)
summary(mf)


### Model diagnostics ##########################################################

# Initial Goodnes of fit test
dev <- summary(m1)$deviance
df <- summary(m1)$df.resid
p <- 1 - pchisq(dev, df)
p # p = 0.98 so the data do not reject the initial model.


## Model diagnostic plots ##
par(mfrow=c(2,2))
plot(m1, which = 1:4)
par(mfrow=c(1,1))

# Plot of deviance residuals against the different factors
stddevresid <- rstandard(m1, type = "deviance")
par(mfrow=c(2,2))
plot(swimmer, stddevresid)
abline(h=0, col="red")
plot(location, stddevresid)
abline(h=0, col="red")
plot(age, stddevresid)
abline(h=0, col="red")
plot(sex, stddevresid)
abline(h=0, col="red")
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
