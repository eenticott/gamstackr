
#' @title Example Usage of gamstackr Stacking Functions
#' @description
#' This script demonstrates example workflows for using the gamstackr package, including
#' expert model fitting, stacking, and evaluation. All code is commented out for reference.
#'
#' @details
#' The script covers:
#' - Fitting multiple linear models as experts
#' - Creating prediction matrices
#' - Using LossStack and NestedStack families
#' - Evaluating mean squared and absolute errors
#' - Advanced stacking with density matrices
#' - Timing and benchmarking stacking approaches
#' - Derivative checking for custom stacking families
#'
#' To use, uncomment relevant sections and run in an interactive R session.
#'
#' @seealso LossStack, NestedStack, MVN_weights, ordinal, pinball_loss
#'
#' @examples
#' # See code sections below for usage examples

# --- Load libraries and data ---
# library(gamstackr)   # Main package
# library(mgcv)        # For GAM fitting
# data("iris")         # Example dataset

# --- Fit expert models ---
# ex1 <- lm(Sepal.Length ~ Sepal.Width, data = iris)
# ex2 <- lm(Sepal.Length ~ Petal.Length, data = iris)
# ex3 <- lm(Sepal.Length ~ Petal.Width, data = iris)
# ex4 <- lm(Sepal.Length ~ Petal.Width + Sepal.Width, data = iris)
# ex5 <- lm(Sepal.Length ~ Petal.Width + Petal.Length, data = iris)

# --- Generate predictions matrix ---
# pred1 <- predict(ex1)
# pred2 <- predict(ex2)
# pred3 <- predict(ex3)
# pred4 <- predict(ex4)
# pred5 <- predict(ex5)
# preds <- cbind(pred1, pred2, pred3, pred4, pred5)

# --- Example: LossStack with nested weights ---
# pre_fam <- LossStack(preds, loss=pinball_loss(0.5), nested(multinomial(2), list(ordinal(4), id())))
# m1 <- gam(formula = list(Sepal.Length ~ Sepal.Width,~Sepal.Width), family = pre_fam, data = iris)
# # Calculate mean squared and absolute errors for predictions
# MSE1 <- mean((rowSums(predict(m1, type = "response") * preds) - iris$Sepal.Length)^2)
# MAE1 <- mean(abs(rowSums(predict(m1, type = "response") * preds) - iris$Sepal.Length))

# --- Example: LossStack with ordinal weights ---
# pre_fam <- LossStack(preds, loss=pinball_loss(0.5), ordinal(5))
# m2 <- gam(formula = list(Sepal.Length ~ 1), family = pre_fam, data = iris)
# MSE2 <- mean((rowSums(predict(m2, type = "response") * preds) - iris$Sepal.Length)^2)
# MAE2 <- mean(abs(rowSums(predict(m2, type = "response") * preds) - iris$Sepal.Length))
# c(MSE1, MAE1, MSE2, MAE2)

# --- Example: Density stacking ---
# den1 <- matrix(dnorm(iris$Sepal.Length, predict(ex1), sd(ex1$residuals), log = TRUE))
# den2 <- matrix(dnorm(iris$Sepal.Length, predict(ex2), sd(ex2$residuals), log = TRUE))
# den3 <- matrix(dnorm(iris$Sepal.Length, predict(ex3), sd(ex3$residuals), log = TRUE))
# den4 <- matrix(dnorm(iris$Sepal.Length, predict(ex4), sd(ex4$residuals), log = TRUE))
# den5 <- matrix(dnorm(iris$Sepal.Length, predict(ex5), sd(ex5$residuals), log = TRUE))
# stack_fam = NestedStack(list(cbind(den1, den2, den3, den4, den5)), list(ordinal(5)))
# stack <- gam(list(Sepal.Length ~ Species + Sepal.Width), data = iris, family = stack_fam)

# --- Advanced stacking and benchmarking ---
# (See original script for more timing and derivative checking examples)
# # Example for derivative checking (advanced):
# # deriv_check <- function(pars) {
# #   list_of_betaT <- list(list(pars[1:2]))
# #   list_of_theta <- list(pars[3:5])
# #   get_derivatives(list_of_beta,
# #                   list_of_betaT,
# #                   list_of_theta,
# #                   list_of_X_eta,
# #                   list_of_X_etaT,
# #                   list(ordinal(5)),
# #                   list_of_densities = list(exp(cbind(den1, den2, den3, den4, den5))))$ll
# # }
# # dv_estg <- numDeriv::grad(deriv_check, c(1,1,1,2,3))
# # dv_esth <- numDeriv::hessian(deriv_check, c(1,1,1,2,3))
# # dv_calc$llb/dv_estg
# # dv_calc$llbb/dv_esth

# # Example for checking derivatives of ordinal weights:
# # f <- function(pars) {
# #   x <- c(1,1)
# #   beta <- pars[1:2]
# #   theta <- pars[3]
# #   fg <- ordinal(3)
# #   eta <- sum(beta*x)
# #   fg(list(matrix(eta)), theta,deriv=1)$f_eval[2]
# # }
# # numDeriv::hessian(f, c(2,1,1))
# # str(fg(list(matrix(c(3))), 1,deriv=2))
# # fg(list(matrix(c(1,1),ncol=2)), c(1,1),deriv=1)
# # fgc <- function(pars) {
# #   fg(matrix(pars[1:2],ncol=2), pars[3:4],deriv=1)$f_eval[2]
# # }
# # numDeriv::grad(fgc, c(1,1,1,1))
