# library(gamstackr)
# library(mgcv)
#
# data("iris")
#
# ex1 <- lm(Sepal.Length ~ Sepal.Width, data = iris)
# ex2 <- lm(Sepal.Length ~ Petal.Length, data = iris)
# ex3 <- lm(Sepal.Length ~ Petal.Width, data = iris)
#
# den1 <- matrix(dnorm(iris$Sepal.Length, predict(ex1), sd(ex1$residuals), log = TRUE))
# den2 <- matrix(dnorm(iris$Sepal.Length, predict(ex2), sd(ex2$residuals), log = TRUE))
# den3 <- matrix(dnorm(iris$Sepal.Length, predict(ex3), sd(ex3$residuals), log = TRUE))
#
# stack_fam = NestedStack(list(den1, den2, den3), list(id(), id(), id()))
# stack_fam = NestedStack(list(cbind(den1, den2, den3)), list(ordinal(3)))
# stack_fam = NestedStack(list(cbind(den1, den2, den3)), list(MVN_weights(x=matrix(c(1,2,3),nrow=1))))
# stack <- gam(list(Sepal.Length ~ Species + Sepal.Width), data = iris, family = stack_fam)
#
#
# stack <- gam(list(Sepal.Length ~ Species + Sepal.Width,  ~Species), data = iris, family = stack_fam)
#
# W <- predict(stack, type = "response")
#
# barplot(t(unique(W[,1:3])), col = unique(iris$Species),
#         names.arg = unique(iris$Species), xlab = "Species", ylab = "Weight")
# legend("topright", legend = c("Expert 1", "Expert 2", "Expert 3"), fill = unique(iris$Species))
#
# list_of_betaT <- list(list(c(1,1)))
# list_of_X_etaT <- list(list(cbind(1, matrix(rnorm(nrow(iris)),ncol=1))))
# list_of_beta <- list()
# list_of_X_eta <- list()
# list_of_theta <- list(c(1))
# dv_calc <- get_derivatives(list_of_beta,
#                 list_of_betaT,
#                 list_of_theta,
#                 list_of_X_eta,
#                 list_of_X_etaT,
#                 list(MVN_weights(x=matrix(c(1,2,3),nrow=1))),
#                 list_of_densities = list(exp(cbind(den1, den2, den3))))
#
# deriv_check <- function(pars) {
#   list_of_betaT <- list(list(pars[1:2]))
#   list_of_theta <- list(pars[3])
#   get_derivatives(list_of_beta,
#                   list_of_betaT,
#                   list_of_theta,
#                   list_of_X_eta,
#                   list_of_X_etaT,
#                   list(MVN_weights(x=matrix(c(1,2,3),nrow=1))),
#                   list_of_densities = list(exp(cbind(den1, den2, den3))))$ll
# }
#
# dv_est <- numDeriv::hessian(deriv_check, c(1,1,1))
#
# dv_calc$llbb/dv_est
#
#
# f <- function(pars) {
#   x <- c(1,1)
#   beta <- pars[1:2]
#   theta <- pars[3]
#   fg <- MVN_weights(x=matrix(c(1,2,3),nrow=1,ncol=3))
#   eta <- sum(beta*x)
#   fg(list(matrix(eta)), theta,deriv=1)$f_eval[1]
# }
#
# numDeriv::hessian(f, c(2,1,1))
# str(fg(list(matrix(c(3))), 1,deriv=2))
# fg(list(matrix(c(1,1),ncol=2)), c(1,1),deriv=1)
#
# fgc <- function(pars) {
#   fg(matrix(pars[1:2],ncol=2), pars[3:4],deriv=1)$f_eval[2]
# }
#
# numDeriv::grad(fgc, c(1,1,1,1))
#
# # Check ordinal
#
# list_of_betaT <- list(list(c(1,1)))
# list_of_X_etaT <- list(list(cbind(1, matrix(rnorm(nrow(iris)),ncol=1))))
# list_of_beta <- list()
# list_of_X_eta <- list()
# list_of_theta <- list(c(1))
# dv_calc <- get_derivatives(list_of_beta,
#                            list_of_betaT,
#                            list_of_theta,
#                            list_of_X_eta,
#                            list_of_X_etaT,
#                            list(ordinal(3)),
#                            list_of_densities = list(exp(cbind(den1, den2, den3))))
#
# deriv_check <- function(pars) {
#   list_of_betaT <- list(list(pars[1:2]))
#   list_of_theta <- list(pars[3])
#   get_derivatives(list_of_beta,
#                   list_of_betaT,
#                   list_of_theta,
#                   list_of_X_eta,
#                   list_of_X_etaT,
#                   list(ordinal(3)),
#                   list_of_densities = list(exp(cbind(den1, den2, den3))))$ll
# }
#
# dv_estg <- numDeriv::grad(deriv_check, c(1,1,1))
# dv_esth <- numDeriv::hessian(deriv_check, c(1,1,1))
# dv_calc$llb/dv_estg
# dv_calc$llbb/dv_esth
#
#
# f <- function(pars) {
#   x <- c(1,1)
#   beta <- pars[1:2]
#   theta <- pars[3]
#   fg <- ordinal(3)
#   eta <- sum(beta*x)
#   fg(list(matrix(eta)), theta,deriv=1)$f_eval[2]
# }
#
# numDeriv::hessian(f, c(2,1,1))
# str(fg(list(matrix(c(3))), 1,deriv=2))
# fg(list(matrix(c(1,1),ncol=2)), c(1,1),deriv=1)
#
# fgc <- function(pars) {
#   fg(matrix(pars[1:2],ncol=2), pars[3:4],deriv=1)$f_eval[2]
# }
#
# numDeriv::grad(fgc, c(1,1,1,1))
#
