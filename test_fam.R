load("Data/GEF14.RData")
library(tidyverse)
library(mgcv)
library(mgcViz)
head(GEF14)

dTrain <- subset(GEF14, year <= 2008)
dStack <- subset(GEF14, year > 2008 & year <= 2010)
dTest  <- subset(GEF14, year == 2011)


# Fit experts ------------------------------------------------------------------------
fit1 <- bam(load ~ s(load24, k = 8, bs = "cr") + dow, discrete = T, data = dTrain)

summary(fit1)

fit2 <- bam(load ~ s(load24, k = 8, bs = "cr") +
              dow + s(tod, k = 12, bs = "cc") +
              s(toy, k = 20, bs = "cc"), discrete = T, data = dTrain)

summary(fit2)


fit3 <- bam(load ~ s(load24, k = 8, bs = "cr") +
              s(tod, k = 12, bs = "cc") +
              s(toy, k = 20, bs = "cc") +
              s(temp, k = 20, bs = "cr") +
              s(temp95, k = 20, bs = "cr") +
              dow, discrete = T, data = subset(dTrain, toy < 0.33))

fit4 <- bam(load ~ s(load24, k = 8, bs = "cr") +
              s(tod, k = 12, bs = "cc") +
              s(toy, k = 20, bs = "cc") +
              s(temp, k = 20, bs = "cr") +
              s(temp95, k = 20, bs = "cr") +
              dow, discrete = T, data = subset(dTrain, toy < 0.66 & toy >= 0.33))

fit5 <- bam(load ~ s(load24, k = 8, bs = "cr") +
              s(tod, k = 12, bs = "cc") +
              s(toy, k = 20, bs = "cc") +
              s(temp, k = 20, bs = "cr") +
              s(temp95, k = 20, bs = "cr") +
              dow, discrete = T, data = subset(dTrain, toy >= 0.66))

p1 <- predict(fit1, newdata = dStack)
p2 <- predict(fit2, newdata = dStack)
p3 <- predict(fit3, newdata = dStack)
p4 <- predict(fit4, newdata = dStack)
p5 <- predict(fit5, newdata = dStack)

den1 <- dnorm(dStack$load, p1, sqrt(fit1$sig2), log = FALSE)
den2 <- dnorm(dStack$load, p2, sqrt(fit2$sig2), log = FALSE)
den3 <- dnorm(dStack$load, p3, sqrt(fit3$sig2), log = FALSE)
den4 <- dnorm(dStack$load, p4, sqrt(fit4$sig2), log = FALSE)
den5 <- dnorm(dStack$load, p5, sqrt(fit5$sig2), log = FALSE)

dens <- cbind(den3, den4, den5)

inner_ordinal <- ordinal(3)
inners <- list(inner_ordinal)
list_of_densities <- list(dens)
pre_fam <- NestedStack(list_of_densities, inners, RidgePen = 1e-04)

fitStack <- gam(list(load ~ s(toy)), data = dStack, family = pre_fam, control = gam.control(trace = TRUE))

summary(fitStack)

preds <- predict(fitStack, type = "response")
fitStack2 <- getViz(fitStack)

library(mgcViz)
plot(ALE(fitStack2, x = "toy", oind = 1, type = "response"))

predict(fitStack)

