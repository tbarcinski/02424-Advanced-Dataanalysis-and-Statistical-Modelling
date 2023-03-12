rm(list=ls()) 
print(utils::getSrcDirectory(function(){}))
print(utils::getSrcFilename(function(){}, full.names = TRUE))
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()
Sys.setenv(LANG = "en")

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
library(numDeriv)
library(fitdistrplus)
library(lmtest)
library(leaps)
library(plotly)

# Reading data
dat <- read.table("dioxin.csv", sep=',', head=TRUE)
df <- data.frame(dat)
df <- df %>% select(-OBSERV)
summary(df)
is.na(df)
colSums(is.na(df))
# PRSEK, CO, SO2
sum(is.na(df))

cols <- c("PLANT", "TIME", "LAB", "OXYGEN", "LOAD", "PRSEK")
df[cols] <- lapply(df[cols], as.factor)
df <- data.frame(df)
str(df)
attach(df)

########## checking ##############3
AM3 <- lm(log(DIOX) ~ O2COR + PLANT + TIME + LAB + CO2 + TROEG + POVN + 
            poly(log(HCL), 2) + TIME:poly(NEFFEKT, 1), data = data.frame(df))
summary(AM3)
c(AIC(AM3), BIC(AM3))

model_interactions <- lm(log(DIOX) ~ O2COR + PLANT + TIME + LAB + CO2 + TROEG + POVN + 
                           poly(log(HCL), 2) + TIME:poly(NEFFEKT, 1) + O2COR*NEFFEKT*QRAT, data = df)

summary(model_interactions)
drop1(model_interactions, test = "F")


boxcox(AM3)

par(mfrow = c(2, 2))
plot(AM3)

par(mfrow = c(2, 1))
hist(df$CO, breaks = 20)
hist(df$SO2, breaks = 20)


fig <- plot_ly(df, x = ~O2COR, y = ~NEFFEKT, z = ~QRAT,
               marker = list(color = ~(DIOX), colorscale = c('#FFE1A1', '#683531'), showscale = TRUE))
fig <- fig %>% add_markers()
fig <- fig %>% layout(scene = list(xaxis = list(title = 'O2COR'),
                                   yaxis = list(title = 'NEFFEKT'),
                                   zaxis = list(title = 'QRAT')),
                      annotations = list(
                        x = 1.13,
                        y = 1.05,
                        text = 'log(DIOX)',
                        xref = 'paper',
                        yref = 'paper',
                        showarrow = FALSE
                      ))
fig

# changing to factors


## checking for normality
boxcox <- boxcox(DIOX ~ 1)
qqnorm(log(df$DIOX))
qqline(log(df$DIOX), col = "red")
shapiro.test(log(df$DIOX))
lillie.test(log(df$DIOX))
# wallyplot R residuals

# Boxplot of the measured active variables vs the planned active variables
g1 <- ggplot(DATA, aes(x=OXYGEN, y=O2COR)) +
  geom_boxplot() +
  theme(axis.text=element_text(size=14,face="bold"),
        axis.title=element_text(size=16,face="bold"),
        legend.position="none") + 
  xlab("OXYGEN") + ylab("O2COR")

g2 <- ggplot(DATA, aes(x=LOAD, y=NEFFEKT)) +
  geom_boxplot() +
  theme(axis.text=element_text(size=14,face="bold"),
        axis.title=element_text(size=16,face="bold"),
        legend.position="none") + 
  xlab("LOAD") + ylab("NEFFEKT")

g3 <- ggplot(DATA, aes(x=PRSEK, y=QRAT)) + geom_boxplot() +
  theme(axis.text=element_text(size=14,face="bold"),
        axis.title=element_text(size=16,face="bold"),
        legend.position="none") + 
  xlab("PRSEK") + ylab("QRAT")

g1 + g2 + g3 + plot_annotation(
  title = 'Discretized vs continuous active variables',
  theme=theme(plot.title=element_text(hjust=0.5, size = 22))
)

# grid.arrange(g1, g2, g3, nrow=1, ncol=3,
#              top=textGrob("Daily QC: Blue",gp=gpar(fontsize=20,font=3)))
# 

# Boxplots
g1 <- ggplot(DATA, aes(x=OXYGEN, y=DIOX)) +
  geom_boxplot() +
  theme(axis.text=element_text(size=14,face="bold"),
        axis.title=element_text(size=16,face="bold"),
        legend.position="none")

g2 <- ggplot(DATA, aes(x=LOAD, y=DIOX)) +
  geom_boxplot() +
  theme(axis.text=element_text(size=14,face="bold"),
        axis.title=element_text(size=16,face="bold"),
        legend.position="none")

g3 <- ggplot(DATA, aes(x=PRSEK, y=DIOX)) + geom_boxplot() +
  theme(axis.text=element_text(size=14,face="bold"),
        axis.title=element_text(size=16,face="bold"),
        legend.position="none")

g1_log <- ggplot(DATA, aes(x=OXYGEN, y=log(DIOX))) +
  geom_boxplot() +
  theme(axis.text=element_text(size=14,face="bold"),
        axis.title=element_text(size=16,face="bold"),
        legend.position="none")

g2_log <- ggplot(DATA, aes(x=LOAD, y=log(DIOX))) +
  geom_boxplot() +
  theme(axis.text=element_text(size=14,face="bold"),
        axis.title=element_text(size=16,face="bold"),
        legend.position="none")

g3_log <- ggplot(DATA, aes(x=PRSEK, y=log(DIOX))) + geom_boxplot() +
  theme(axis.text=element_text(size=14,face="bold"),
        axis.title=element_text(size=16,face="bold"),
        legend.position="none")

(g1 + g2 + g3) / (g1_log + g2_log + g3_log) + plot_annotation(
  title = 'Dependent variable vs active variables in the original and the log domain',
  theme=theme(plot.title=element_text(hjust=0.5, size = 20))
)


g4 <- ggplot(DATA, aes(x=PLANT, y=DIOX, colour=TIME)) +
  geom_boxplot() +
  theme(axis.text=element_text(size=14,face="bold"),
        axis.title=element_text(size=16,face="bold"),
        legend.position="right")

g5 <- ggplot(DATA, aes(x=LAB, y=DIOX)) +
  geom_boxplot() +
  theme(axis.text=element_text(size=14,face="bold"),
        axis.title=element_text(size=16,face="bold"),
        legend.position="none")

g4_log <- ggplot(DATA, aes(x=PLANT, y=log(DIOX), colour=TIME)) +
  geom_boxplot() +
  theme(axis.text=element_text(size=14,face="bold"),
        axis.title=element_text(size=16,face="bold"),
        legend.position="right")

g5_log <- ggplot(DATA, aes(x=LAB, y=log(DIOX))) +
  geom_boxplot() +
  theme(axis.text=element_text(size=14,face="bold"),
        axis.title=element_text(size=16,face="bold"),
        legend.position="none")

(g4 + g5) / (g4_log + g5_log) + plot_annotation(
  title = 'Dependent variable vs block variables in the original and the log domain',
  theme=theme(plot.title=element_text(hjust=0.5, size = 20))
)

grid.arrange(g1, g2, g3, g4, g5, nrow=2, ncol=3)


### EDA
# 1) target vs active variables + block variables
df_part <- df[,c("DIOX", "O2COR", "NEFFEKT", "QRAT", "LAB")]
df_part["diox_log"] <- log(df_part$DIOX)
df_part <- df_part %>% select(-DIOX)
p_lab <- ggpairs(df_part, aes(color = LAB),
                  columns = c("diox_log", "O2COR", "NEFFEKT", "QRAT"))
p_lab

df_part <- df[,c("DIOX", "O2COR", "NEFFEKT", "QRAT", "PLANT")]
df_part["diox_log"] <- log(df_part$DIOX)
df_part <- df_part %>% select(-DIOX)
p_plant <- ggpairs(df_part, aes(color = PLANT, shape = TIME),
                  columns = c("diox_log", "O2COR", "NEFFEKT", "QRAT"))
p_plant

# Time needs to be only for this one plant that was different
df_part <- df %>% filter(PLANT == "RENO_N") %>%
  select("DIOX", "O2COR", "NEFFEKT", "QRAT", "TIME")
df_part["diox_log"] <- log(df_part$DIOX)
df_part <- df_part %>% select(-DIOX)
p_time <- ggpairs(df_part, aes(color = TIME),
                  columns = c("diox_log", "O2COR", "NEFFEKT", "QRAT"))
p_time

p1 <- ggplot(df, aes(x = seq(1, dim(df)[1], 1), y = DIOX)) + geom_point() +
  xlab("Index")
p2 <- ggplot(df, aes(x = DIOX)) + geom_histogram()
p1_log<- ggplot(df, aes(x = seq(1, dim(df)[1], 1), y = log(DIOX))) + geom_point() +
  xlab("Index")
p2_log <- ggplot(df, aes(x = log(DIOX))) + geom_histogram()
(p1 + p2) / (p1_log + p2_log) + plot_annotation(
  title = 'Distribution of the dependent variable in the original and log domain',
  theme=theme(plot.title=element_text(hjust=0.5, size = 18))
)

# 2) target vs active variables & passive variables + block variables

df_full <- df %>% 
  mutate(co_log = log(CO), hcl_log = log(HCL), diox_log = log(DIOX)) %>% 
  select(-OXYGEN, -LOAD, -PRSEK, -O2, -DIOX, -CO, -HCL)

# plant
df_full_plant <- df_full %>% select(-TIME, -LAB) 
p_full_plant <- ggpairs(df_full_plant, aes(color = PLANT),
                        columns = 2:14)
p_full_plant

# lab
df_full_lab <- df_full %>% select(-TIME, -PLANT) 
p_full_lab <- ggpairs(df_full_lab, aes(color = LAB),
                        columns = 2:14)
p_full_lab

# time
df_full_time <- df_full %>% filter(PLANT == "RENO_N") %>%
  select(-LAB, -PLANT) 
p_full_time  <- ggpairs(df_full_time, aes(color = TIME),
                      columns = 2:14)
p_full_time

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

#### Analysis of time
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

### Analysis of active variables discretized vs continous
p4 <- ggplot(df, aes(x = OXYGEN, y = O2COR)) +
  geom_boxplot()
p5 <- ggplot(df, aes(x = LOAD, y = NEFFEKT)) +
  geom_boxplot()
p6 <- ggplot(df, aes(x = PRSEK, y = QRAT)) +
  geom_boxplot()
p4 + p5 + p6
 
df %>% filter(is.na(PRSEK))
# this is just one observation but sent to different LABs

#### DIFFERENT MODELS ####
model = tree(log(DIOX)~., data = df)
plot(model)
text(model)

# To do
# 1) interactions plots between continous and factor variables
# 2) interactions plots
# 3) box plots of the block variables for the diox response
# 4) Do the plots with factors as in the lecture
# MSW plants
# interactions

######## 1) SIMPLE ADDITIVE MODEL ############

table(df$PRSEK, useNA ="ifany")  # Observation 15 and 16 are missing
tapply(df$QRAT, df$PRSEK, mean, useNA="ifany")  # Observations seems to have been "L"
df$PRSEK[15:16] <- "L"
attach(df)

model_1 <- lm(DIOX ~ OXYGEN + LOAD + PRSEK + PLANT + TIME + LAB)
summary(model_1)
par(mfrow = c(2, 2))
plot(model_1, which = 1:4)

boxcox(model_1)

stdresid <- rstandard(model_1)
par(mfrow=c(2, 2))
plot(model_1, which = 1)
plot(model_1, which = 3)
hist(stdresid, main="", probability=TRUE, breaks=10)
curve(dnorm, -3, 3, col="red", lwd=2, add=TRUE)
plot(model_1, which=2)
par(mfrow=c(1, 1))



model_1_log <- lm(log(DIOX) ~ OXYGEN + LOAD + PRSEK + PLANT + TIME + LAB, data = df)
summary(model_1_log)

stdresid <- rstandard(model_1_log)
par(mfrow=c(2, 2))
plot(model_1_log, which = 1)
plot(model_1_log, which = 3)
hist(stdresid, main="", probability=TRUE, breaks=10)
curve(dnorm, -3, 3, col="red", lwd=2, add=TRUE)
plot(model_1_log, which=2)


AIC(model_1_log)

par(mfrow = c(2, 2))
plot(model_1_log, which = 1:4)

par(mfrow = c(1, 1))
acf(residuals(model_1_log))
rstudent(model_1_log)[which.max(abs(rstudent(model_1_log)))]
# 13
drop1(model_1_log, test = "F")

model_2_log <- update(model_1_log,~. -PRSEK, data = df)
summary(model_2_log)
par(mfrow = c(2, 2))
plot(model_2_log, which = 1:4)





anova(model_2, model_1_log, test="Chisq")

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
summary(model_c1)

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

# 17 observations
# missing values - normal distrubtion

#### Predictions
# can we just exp the residuals

# interpretation
exp(model_c2$coefficients$O2COR)
exp(0.18)
exp(2.9)

#### FULL MODEL WITH INTERACTIONS ########
df_con <- df %>% select(-LOAD, -OXYGEN, -PRSEK, -O2)
m1 <- lm(log(DIOX) ~ . + NEFFEKT:PLANT:TIME +
           TROEG:PLANT:TIME +
           log(CO):PLANT:TIME, data = df_con)
summary(m1)


### DATA IMPUTATION
colSums(is.na(df_con))    # SO2: Obs 7, 8, 53 and CO: obs 8

### CO
par(mfrow = c(1, 1))
descdist(log(df_con$CO[!is.na(df_con$CO)]), boot = 1000)
data_co <- log(df_con$CO[!is.na(df_con$CO)])

lnorm_co <- fitdist(data_co, "lnorm")
norm_co <- fitdist(data_co, "norm")
gamma_co <- fitdist(data_co, "gamma")
plot_legend_ws <- c("lnorm", "norm", "gamma")

ws_fitted <- list(lnorm_co, norm_co, gamma_co)
par(mfrow = c(1, 2))
denscomp(ws_fitted, legendtext = plot_legend_ws)
qqcomp(ws_fitted, legendtext = plot_legend_ws)

summary(lnorm_co)
summary(norm_co)
summary(gamma_co)

lnorm_co$aic
norm_co$aic
gamma_co$aic
### we pick gamma distribution E[gamma]= shape/rate
to_input_co = exp(as.numeric(gamma_co$estimate["shape"] / gamma_co$estimate["rate"]))
df_con$CO[is.na(df_con$CO)] = to_input_co

### SO2
par(mfrow = c(1, 1))
descdist(log(df_con$SO2[!is.na(df_con$SO2)]), boot = 1000)
data_so2 <- log(df_con$SO2[!is.na(df_con$SO2)])

lnorm_co <- fitdist(data_so2, "lnorm")
norm_co <- fitdist(data_so2, "norm")
gamma_co <- fitdist(data_so2, "gamma")
plot_legend_ws <- c("lnorm", "norm", "gamma")

ws_fitted <- list(lnorm_co, norm_co, gamma_co)
par(mfrow = c(1, 2))
denscomp(ws_fitted, legendtext = plot_legend_ws)
qqcomp(ws_fitted, legendtext = plot_legend_ws)

summary(lnorm_co)
summary(norm_co)
summary(gamma_co)

lnorm_co$aic
norm_co$aic
gamma_co$aic

to_input_so2 = exp(as.numeric(norm_co$estimate["mean"]))
df_con$SO2[is.na(df_con$SO2)] = to_input_so2

## in the regular domain: Adjusted R-squared:  0.9931
attach(df_con)
m2 <- lm((DIOX) ~ 
           poly(O2COR,2)*PLANT + poly(NEFFEKT,2)*PLANT + poly(QRAT,2)*PLANT + 
           poly(O2COR,2)*LAB + poly(NEFFEKT,2)*LAB + poly(QRAT,1)*LAB +
           poly(O2COR,2)*TIME + poly(NEFFEKT,2)*TIME + poly(QRAT,1)*TIME +
           QROEG + TOVN + TROEG + POVN + CO2 + CO + SO2 + log(HCL) + H2O + 
           QROEG:PLANT + CO2:PLANT +
           log(HCL):LAB)
summary(m2)
par(mfrow = c(2, 2))
plot(m2)

box_object <- boxcox(m2)

lambda <- box_object$x[which.max(box_object$y)]
diox_box_cox <- (df$DIOX ^ lambda - 1) / lambda

m2_boxcox <- lm(diox_box_cox ~ 
                  poly(O2COR,2)*PLANT + poly(NEFFEKT,2)*PLANT + poly(QRAT,2)*PLANT + 
                  poly(O2COR,2)*LAB + poly(NEFFEKT,2)*LAB + poly(QRAT,1)*LAB +
                  poly(O2COR,2)*TIME + poly(NEFFEKT,2)*TIME + poly(QRAT,1)*TIME +
                  QROEG + TOVN + TROEG + POVN + CO2 + CO + SO2 + log(HCL) + H2O + 
                  QROEG:PLANT + CO2:PLANT +
                  log(HCL):LAB, data = df)
summary(m2_boxcox, cor = T)
par(mfrow = c(2, 2))
plot(m2_boxcox)

# weird stuff happens: some observations have a leverage of 1 ...
plot(rstudent(m2))

shapiro.test(rstudent(m2))
lillie.test(rstudent(m2))
cvm.test(rstudent(m2))
ad.test(rstudent(m2))

### REDUCTION
drop1(m2, test="F")
m2 <- update(m2, . ~ . - TROEG)
drop1(m2, test="F")
m2 <- update(m2, . ~ . - LAB:log(HCL) )
drop1(m2, test="F")
m2 <- update(m2, . ~ . - CO)
drop1(m2, test="F")
m2 <- update(m2, . ~ . - SO2)
drop1(m2, test="F")
m2 <- update(m2, . ~ . - poly(NEFFEKT, 2):LAB )
drop1(m2, test="F")
m2 <- update(m2, . ~ . - poly(QRAT, 1):TIME  )
drop1(m2, test="F")
m2 <- update(m2, . ~ . - poly(O2COR, 2):LAB)
drop1(m2, test="F")
m2 <- update(m2, . ~ . - LAB:poly(QRAT, 1))
drop1(m2, test="F")
m2 <- update(m2, . ~ . - poly(O2COR, 2):TIME)
drop1(m2, test="F")
m2 <- update(m2, . ~ . - LAB:TOVN)
drop1(m2, test="F")
m2 <- update(m2, . ~ . - TOVN)
drop1(m2, test="F")
m2 <- update(m2, . ~ . - poly(NEFFEKT, 2):TIME)
drop1(m2, test="F")
m2 <- update(m2, . ~ . - POVN)
drop1(m2, test="F")
m2 <- update(m2, . ~ . - log(HCL))
drop1(m2, test="F")

summary(m2)
par(mfrow = c(2, 2))
plot(m2)
plot(rstudent(m2))
# the model is not neceserily good, and why so many variables? 
c(AIC(m2), BIC(m2))

shapiro.test(rstudent(m2))
lillie.test(rstudent(m2))
cvm.test(rstudent(m2))
ad.test(rstudent(m2))
# the assumptions are not satisfied


####################### STRANGE MODEL ########################################
m4 <- lm(log(DIOX) ~ 
           poly(O2COR,2)*PLANT + poly(NEFFEKT,2)*PLANT + poly(QRAT,1)*PLANT + 
           poly(O2COR,2)*LAB + poly(NEFFEKT,2)*LAB + poly(QRAT,1)*LAB +
           poly(O2COR,1)*TIME + poly(NEFFEKT,1)*TIME + poly(QRAT,1)*TIME +
           TOVN + TROEG + POVN + CO2 + CO + SO2 + log(HCL) + H2O + 
           QROEG:PLANT + CO2:PLANT + TIME:PLANT)
summary(m4)
par(mfrow = c(2, 2))
plot(m4)
par(mfrow = c(1, 1))
plot(rstudent(m4))
# 15, 16, 17, 18, 19 - do they have leverage 1?

### REDUCTION
drop1(m4, test="F")
m4 <- update(m4, . ~ . - poly(QRAT, 1):TIME)
drop1(m4, test="F")
m4 <- update(m4, . ~ . - POVN)
drop1(m4, test="F")
m4 <- update(m4, . ~ . - poly(QRAT, 1):LAB)
drop1(m4, test="F")
m4 <- update(m4, . ~ . - poly(O2COR, 2):LAB)
drop1(m4, test="F")
m4 <- update(m4, . ~ . - H2O)
drop1(m4, test="F")
m4 <- update(m4, . ~ . - poly(NEFFEKT, 2):LAB)
drop1(m4, test="F")
m4 <- update(m4, . ~ . - PLANT:CO2)
drop1(m4, test="F")
m4 <- update(m4, . ~ . - PLANT:TIME)
drop1(m4, test="F")
m4 <- update(m4, . ~ . - CO2)
drop1(m4, test="F")
m4 <- update(m4, . ~ . - poly(O2COR, 1):TIME)
drop1(m4, test="F")
m4 <- update(m4, . ~ . - TIME:poly(NEFFEKT, 1))
drop1(m4, test="F")
m4 <- update(m4, . ~ . - CO)
drop1(m4, test="F")

summary(m4)
par(mfrow = c(2, 2))
plot(m4)

c(AIC(m4), BIC(m4))



#################### SORENS MODEL #################################
AM1 <- lm(log(DIOX) ~ poly(O2COR,2)*PLANT + poly(NEFFEKT,2)*PLANT + poly(QRAT,2)*PLANT + 
            poly(O2COR,2)*LAB + poly(NEFFEKT,2)*LAB + poly(QRAT,2)*LAB +
            poly(O2COR,2)*TIME + poly(NEFFEKT,2)*TIME + poly(QRAT,2)*TIME +
            QROEG + TOVN + TROEG + POVN + CO2 + CO + SO2 + HCL + H2O, data = df_con)
summary(AM1)

par(mfrow = c(2, 2))
plot(AM1)

AM1 <- update(AM1, . ~ . - TROEG)
drop1(AM1, test="F")

AM1 <- update(AM1, . ~ . - poly(QRAT, 2):TIME)
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
AM1 <- update(AM1, . ~ . - poly(NEFFEKT, 1):TIME)
drop1(AM1, test="F")

summary(AM1)
par(mfrow = c(2, 2))
plot(AM1)
plot(rstudent(AM1))

box_object <- boxcox(m2)

lambda <- box_object$x[which.max(box_object$y)]
diox_box_cox <- (df$DIOX ^ lambda - 1) / lambda

m2_boxcox <- lm(diox_box_cox ~ 
                  poly(O2COR,2)*PLANT + poly(NEFFEKT,2)*PLANT + poly(QRAT,2)*PLANT + 
                  poly(O2COR,2)*LAB + poly(NEFFEKT,2)*LAB + poly(QRAT,1)*LAB +
                  poly(O2COR,2)*TIME + poly(NEFFEKT,2)*TIME + poly(QRAT,1)*TIME +
                  QROEG + TOVN + TROEG + POVN + CO2 + CO + SO2 + log(HCL) + H2O + 
                  QROEG:PLANT + CO2:PLANT +
                  log(HCL):LAB, data = df)

#################### SORENS MODEL UPDATED #################################
m5 <- lm(log(DIOX) ~ poly(O2COR,2)*PLANT + poly(NEFFEKT,2)*PLANT + poly(QRAT,2)*PLANT +
           poly(O2COR,2)*LAB + poly(NEFFEKT,2)*LAB + poly(QRAT,2)*LAB +
           poly(O2COR,2)*TIME + poly(NEFFEKT,2)*TIME + poly(QRAT,2)*TIME +
           QROEG + TOVN + TROEG + POVN + CO2 + CO + SO2 + log(HCL) + H2O + QROEG:PLANT
           ,
         data = df_con)

m5_before <- lm((DIOX) ~ poly(O2COR,2)*PLANT + poly(NEFFEKT,2)*PLANT + poly(QRAT,2)*PLANT +
           poly(O2COR,2)*LAB + poly(NEFFEKT,2)*LAB + poly(QRAT,2)*LAB +
           poly(O2COR,2)*TIME + poly(NEFFEKT,2)*TIME + poly(QRAT,2)*TIME +
           QROEG + TOVN + TROEG + POVN + CO2 + CO + SO2 + log(HCL) + H2O + QROEG:PLANT
         ,
         data = df_con)
summary(m5_before)
par(mfrow = c(2, 2))
plot(m5_before, which = 1:4)

box_object <- boxcox(m5_before)

lambda <- box_object$x[which.max(box_object$y)]
diox_box_cox <- (df$DIOX ^ lambda - 1) / lambda


summary(m5)
par(mfrow = c(2, 2))
plot(m5, which = 1:4)

shapiro.test(rstudent(m5))
lillie.test(rstudent(m5))
cvm.test(rstudent(m5))
ad.test(rstudent(m5))

### REDUCTION
drop1(m5, test="F")
m5 <- update(m5, . ~ . - POVN)
drop1(m5, test="F")
m5 <- update(m5, . ~ . - CO2)
drop1(m5, test="F")
m5 <- update(m5, . ~ . - H2O)
drop1(m5, test="F")
m5 <- update(m5, . ~ . - poly(O2COR, 2):LAB)
drop1(m5, test="F")
m5 <- update(m5, . ~ . - poly(O2COR, 2):PLANT)
drop1(m5, test="F")
m5 <- update(m5, . ~ . - poly(QRAT, 2):TIME)
drop1(m5, test="F")
m5 <- update(m5, . ~ . - SO2)
drop1(m5, test="F")
m5 <- update(m5, . ~ . - TOVN)
drop1(m5, test="F")
m5 <- update(m5, . ~ . - PLANT:poly(NEFFEKT, 2))
drop1(m5, test="F")
m5 <- update(m5, . ~ . - poly(QRAT, 2):LAB)
drop1(m5, test="F")
m5 <- update(m5, . ~ . - CO)
drop1(m5, test="F")
m5 <- update(m5, . ~ . - PLANT:QROEG )
drop1(m5, test="F")
m5 <- update(m5, . ~ . - QROEG)
drop1(m5, test="F")
m5 <- update(m5, . ~ . - poly(NEFFEKT, 2):TIME)
drop1(m5, test="F")
m5 <- update(m5, . ~ . - PLANT:poly(QRAT, 2))
drop1(m5, test="F")
m5 <- update(m5, . ~ . - poly(QRAT, 2))
drop1(m5, test="F")
m5 <- update(m5, . ~ . -poly(O2COR, 2):TIME)

summary(m5) #, cor = TRUE
par(mfrow = c(2, 2))
plot(m5)

######### MODEL VALIDATION ############
# residuals
shapiro.test(rstudent(m5))
lillie.test(rstudent(m5))
cvm.test(rstudent(m5))
ad.test(rstudent(m5))

par(mfrow = c(1, 1))
plot(rstudent(m5))
abline(h = 0, lty = 2, col = "blue")
abline(h = 2, lty = 2, col = "blue")
abline(h = -2, lty = 2, col = "blue")
sort(rstudent(m5))

par(mfrow = c(2, 2))
plot(df_con$O2COR, rstandard(m5))
plot(df_con$NEFFEKT, rstandard(m5))
plot(df_con$TROEG, rstandard(m5))
plot(log(df_con$HCL), rstandard(m5))

par(mfrow = c(2, 2))
plot(m5, which=1:4)
# observation number 17 ...
# bptest(m5)
# bptest(AM1) # hmm

#################### COMPARISON ##################################
c(AIC(AM3), BIC(AM3))
# c(AIC(m5), BIC(m5))
c(AIC(model_2_log), BIC(model_2_log))

# models are nested, the reduction can not be performed
anova(m5, AM1, model_c1)


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

#############################################################
# df_without_missing <- na.omit(df_con)
# attach(df_without_missing)
AM3 <- lm(log(DIOX) ~ O2COR  + NEFFEKT + PLANT + TIME + LAB, data = df)

scope <- ~ . + QROEG + TOVN + TROEG + POVN + CO2 + CO + SO2 + HCL + H2O +       # Passive variables
  poly(O2COR,2) + poly(NEFFEKT,2) +  poly(QRAT,1) + poly(QRAT,2) +              # Active variables higher order terms
  PLANT:poly(O2COR,1)   + PLANT:poly(O2COR,2) +                                 # PLANT interactions
  PLANT:poly(NEFFEKT,1) + PLANT:poly(NEFFEKT,2) + 
  PLANT:poly(QRAT,1)    + PLANT:poly(QRAT,2) + 
  LAB:poly(O2COR,1)     + PLANT:poly(O2COR,2) +                                 # LAB interactions
  LAB:poly(NEFFEKT,1)   + PLANT:poly(NEFFEKT,2) + 
  LAB:poly(QRAT,1)      + PLANT:poly(QRAT,2) + 
  TIME:poly(O2COR,1)    + TIME:poly(O2COR,2) +                                  # TIME interactions
  TIME:poly(NEFFEKT,1)  + TIME:poly(NEFFEKT,2) + 
  TIME:poly(QRAT,1)     + TIME:poly(QRAT,2) + 
  poly(log(HCL),1) +  poly(log(HCL),2) +                                        # Higher order effects of log(HCL) and log(CO)
  PLANT:poly(log(HCL),1) +  PLANT:poly(log(HCL),2) +
  LAB:poly(log(HCL),1) +  LAB:poly(log(HCL),2) + 
  TIME:poly(log(HCL),1) +  TIME:poly(log(HCL),2) + 
  poly(log(CO),1) +  poly(log(CO),2) +  
  PLANT:poly(log(CO),1) +  PLANT:poly(log(CO),2) +
  LAB:poly(log(CO),1) +  LAB:poly(log(CO),2) + 
  TIME:poly(log(CO),1) +  TIME:poly(log(CO),2)  + 
  (O2COR*NEFFEKT*QRAT)^2


add1(AM3, scope, test = "F")

AM3 <- update(AM3, . ~ . + poly(log(HCL), 1))
drop1(AM3, test = "F")
add1(AM3, scope, test = "F")

AM3 <- update(AM3, . ~ . + CO2)
drop1(AM3, test = "F")
add1(AM3, scope, test = "F")

AM3 <- update(AM3, . ~ . + TROEG)
drop1(AM3, test = "F")
add1(AM3, scope, test = "F")

AM3 <- update(AM3, . ~ . + TIME:poly(NEFFEKT, 1))
drop1(AM3, test = "F")
add1(AM3, scope, test = "F")

AM3 <- update(AM3, . ~ . + POVN)
drop1(AM3, test = "F")
add1(AM3, scope, test = "F")

AM3 <- update(AM3, . ~ . + poly(log(HCL), 2))
drop1(AM3, test = "F")
add1(AM3, scope, test = "F")

summary(AM3)
AM3 <- update(AM3, . ~ . - poly(log(HCL), 1) + poly(log(HCL), 2) - 
                NEFFEKT + TIME:poly(NEFFEKT, 1))
summary(AM3)


######## LIKELIHOOD ESTIMATION OF WEIGHTS #################
### testing stuff
# par(mfrow = c(1, 1))
# model_test = lm(Volume ~ Girth + Height, data = trees)
# n = dim(trees)[1]
# p = dim(trees)[2]
# y = trees$Volume
# 
# boxcox(model_test)
# 
# design_matrix = model.matrix(model_test)
# 
# objective <- function(theta){
#   y_hat = design_matrix %*% theta[1:p]
#   sum = 0
#   for(i in 1:n){
#     sum = sum + dnorm(y[i], mean = y_hat[i], sd = sqrt(theta[p+1]),log=TRUE)
#   }
#   return(-sum)
# }
# theta = c("theta_0" = 0, "theta_1" = 0, "theta_3" = 0, "sigma_squared" = 1)
# 
# opt <- nlminb(theta,objective)
# opt$par
# 
# objective_2 <- function(theta){
#   y_hat = design_matrix %*% theta[1:p]
#   result = sum(dmvnorm(y, mean = y_hat, sigma = theta[p+1]*diag(n), log = TRUE))
#   return(-result)
# }
# theta = c("theta_0" = 0, "theta_1" = 0, "theta_3" = 0, "sigma_squared" = 1)
# 
# opt_2 <- nlminb(theta, objective_2)
# opt_2$par

AM3 <- lm(log(DIOX) ~ O2COR + PLANT + TIME + LAB + CO2 + TROEG + POVN + 
            poly(log(HCL), 2) + TIME:poly(NEFFEKT, 1), data = data.frame(df))
summary(AM3)
c(AIC(AM3), BIC(AM3))


#### doing it on our data
model_optimization <- AM3

design_matrix = model.matrix(model_optimization)
n = dim(design_matrix)[1]
p = dim(design_matrix)[2]
y = log(df$DIOX)

objective = function(theta){
  y_hat = design_matrix %*% theta[1:p]
  Sigma_weighted = diag(n)
  diag(Sigma_weighted)[design_matrix[, "LABUSA"] == 0] = 1/theta[p+2]
  # print(diag(Sigma_weighted))
  # diag(Sigma_weighted)[design_matrix[, "LABUSA"] == 1] = (1 - theta[p+2])
  Sigma_weighted = theta[p+1]*Sigma_weighted
  result = sum(dmvnorm(y, mean = y_hat, sigma = Sigma_weighted, log = TRUE))
  return(-result)
}
theat_initial = c(rep(1, length(model_optimization$coefficients)), "sigma_squared" = 1,
                  "weight" = 1)
opt <- nlminb(theat_initial, objective)
opt$par

hessian_matrix = hessian(objective,opt$par)
standard_error = sqrt(diag(solve(hessian_matrix)))

c(opt$par["weight"] - qnorm(0.975)*standard_error[length(standard_error)],
  opt$par["weight"] + qnorm(0.975)*standard_error[length(standard_error)])

### PROFILE LIKELIHOOD CI
profile_objective <- function(weight){
  fun.tmp <- function(theta_inner, weight_input){
    objective(c(theta_inner, weight_input))
  }
  theat_initial_inner_opt = c(rep(1, length(model_optimization$coefficients)),
                              "sigma_squared" = 1)
  nlminb(theat_initial_inner_opt, fun.tmp, weight_input = weight)$objective
}
profile_objective(10)

p1 <- seq(-0.5, opt$par["weight"]*2.7, by = 0.01)
logLp1 <- sapply(p1, profile_objective) ## note sapply!
logLp1 <- logLp1 - min(logLp1) ## normalization
L_CI_lower = min(p1[exp(-logLp1) > exp(-qchisq(0.95,df=1)/2)])
L_CI_upper =  max(p1[exp(-logLp1) > exp(-qchisq(0.95,df=1)/2)])

obs_hessian <- hessian(profile_objective, opt$par["weight"])[1, 1]
quadratic_approximation <-
  exp( -0.5* obs_hessian * (p1 - opt$par["weight"])^2)
quadratic_approximation <- quadratic_approximation / max(quadratic_approximation)

par(mfrow = c(1,1))
plot(p1, exp(-logLp1), type = "l",
     xlab="Weight", ylab="Profile Likelihood",
     main="The comparison between Wald's and Likelihood based Confidence Intervals")
axis(side=1, at=seq(0, 10, by=1))
lines(p1, rep(exp(-qchisq(0.95,df=1)/2), length(p1)), col = 2)
rug(L_CI_lower, ticksize = 0.1, lwd = 2, col = "red")
rug(L_CI_upper, ticksize = 0.1, lwd = 2, col = "red")
c(L_CI_lower, L_CI_upper)
lines(p1, quadratic_approximation, col = "blue")
abline(v = opt$par["weight"], lty = 2)
rug(opt$par["weight"] - qnorm(0.975)*standard_error[length(standard_error)],
    ticksize = 0.1, lwd = 2, col = "blue")
rug(opt$par["weight"] + qnorm(0.975)*standard_error[length(standard_error)],
    ticksize = 0.1, lwd = 2, col = "blue")
legend("topright", 95,
       legend=c("Profile likelihood", "Quadratic approximation",
                "95% confidence interval"),
       col=c("black", "blue", "red"), lty = 1:1, cex=0.8,
       inset = 0.02)

c(opt$par["weight"] - qnorm(0.975)*standard_error[length(standard_error)],
  opt$par["weight"] + qnorm(0.975)*standard_error[length(standard_error)])
c(L_CI_lower, L_CI_upper)

wald_statistic = (opt$par["weight"] - 1)/standard_error[length(standard_error)]

# LRT
ll_full <- -profile_objective(as.numeric(opt$par["weight"]))
ll_test <- -profile_objective(1)
LRT <- -2*(ll_test - (ll_full))
p <- 1 - pchisq(LRT, df = 1)
p         


### GLS

gls_model <- gls(log(DIOX) ~ O2COR + NEFFEKT + PLANT + TIME + LAB,
                 weights = varIdent(form = ~ 1 | LAB), data = df, method = "REML")
gls_lm <- gls(log(DIOX) ~ O2COR + NEFFEKT + PLANT + TIME + LAB, data = df)
anova(gls_model, model_c2)
anova(gls_model, gls_lm)
summary(gls_model)




### alternative parametrization
objective_alternative = function(theta){
  y_hat = design_matrix %*% theta[1:p]
  Sigma_weighted = diag(n)
  diag(Sigma_weighted)[design_matrix[, "LABUSA"] == 0] = theta[p+1]
  diag(Sigma_weighted)[design_matrix[, "LABUSA"] == 1] = theta[p+2]
  
  result = sum(dmvnorm(y, mean = y_hat, sigma = Sigma_weighted, log = TRUE))
  return(-result)
}
# theat_initial_alternative = c(rep(1, length(model_c2$coefficients)),
#                   "sigma_squared_KK" = 1, "sigma_squared_USA" = 1)
theat_initial_alternative = c(m5$coefficients,
                              "sigma_squared_KK" = 1, "sigma_squared_USA" = 1)
opt_alternative <- nlminb(theat_initial_alternative, objective_alternative)
opt_alternative$par

hessian_matrix_alternative = hessian(objective_alternative, opt_alternative$par)
standard_error_alternative = sqrt(diag(solve(hessian_matrix_alternative)))

# ctrl alt t
