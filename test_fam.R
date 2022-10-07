load("Data/GEF14.RData")
library(tidyverse)
library(mgcv)
head(GEF14)

dTrain <- subset(GEF14, year <= 2008)
dStack <- subset(GEF14, year > 2008 & year <= 2010)
dTest  <- subset(GEF14, year == 2011)

dTrain <- dTrain %>% mutate_if(is.numeric,~scale(.) %>% as.vector)
dStack <- dStack %>% mutate_if(is.numeric,~scale(.) %>% as.vector)
dTest <- dTest %>% mutate_if(is.numeric,~scale(.) %>% as.vector)



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
              dow, discrete = T, data = dTrain)

summary(fit3)

dStack <- dStack[1:2000,]

p1 <- predict(fit1, newdata = dStack)
p2 <- predict(fit2, newdata = dStack)
p3 <- predict(fit3, newdata = dStack)

den1 <- dnorm(dStack$load, p1, sqrt(fit1$sig2), log = FALSE)
den2 <- dnorm(dStack$load, p2, sqrt(fit2$sig2), log = FALSE)
den3 <- dnorm(dStack$load, p3, sqrt(fit3$sig2), log = FALSE)

dens <- cbind(den1, den2, den3)
list_of_densities <- list(dens)

#list_of_beta <- list()
#list_of_betaT <- list(c(0))
#list_of_theta <- list(1)

#list_of_X_eta <- list()
#list_of_X_etaT <- (list(matrix(rep(1, 1000), nrow = 1000)))


#fitStack <- gam(list(load ~ 1 + s(temp)), data = dStack, family = NestedStack(list_of_densities, inners), control = list(trace = T))

#G = NestedStack(list_of_densities, inners)
#G$Sl <- Sl.setup(G)
#G$X <- Sl.initial.repara(G$Sl, G$X, both.sides = FALSE)


m1 = mean((predict(fit1, newdata = dStack) - dStack$load)^2)
m2 = mean((predict(fit2, newdata = dStack) - dStack$load)^2)
m3 = mean((predict(fit3, newdata = dStack) - dStack$load)^2)

v1 = sqrt(var((predict(fit1, newdata = dStack) - dStack$load)^2))
v2 = sqrt(var((predict(fit2, newdata = dStack) - dStack$load)^2))
v3 = sqrt(var((predict(fit3, newdata = dStack) - dStack$load)^2))

x <- matrix(c(m1, m2, m3, v1, v2, v3), nrow = 2, byrow = T)
inner_MVN <- MVN_weights(x, 2)
inner_ordinal <- ordinal(3)
inner_id <- id(1, nrow(dens))
inners <- list(inner_ordinal, inner_MVN)

list_of_densities <- list(dens, dens)
pre_fam <- NestedStack(list_of_densities, inners, RidgePen = 1e-04)

fitStack <- gam(list(load ~ s(tod), ~ s(temp), ~1, ~1), data = dStack, family = pre_fam, control = gam.control(trace = TRUE))

summary(fitStack)

preds <- predict(fitStack, type = "response")

preds
