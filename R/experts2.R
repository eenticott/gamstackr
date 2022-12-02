
fit_func <- function(data) {
  lm(y ~ x1 + x2, data = data)
}

sim_func <- function(fitted_model, data) {
  # returns simulated response from given data
  mu <- predict(fitted_model, newdata = data)
  sd <- sqrt(as.numeric(crossprod(fitted_model$residuals)/fitted_models$df.residuals))
  rnorm(n = nrow(data), mu, sd = sd)
}

dens_func <- function(fitted_model, data) {
  # returns likelihood of each data point given model
  mu <- predict(fitted_model, newdata = data)
  sd <- sqrt(as.numeric(crossprod(fitted_model$residuals)/fitted_models$df.residuals))
  dnorm(n = nrow(data), mu, sd = sd, log = TRUE)
}

create_expert <- function(fit_func, pred_func = NULL, sim_func = NULL, dens_func = NULL) {

  out_list <- list()

  # Initialise from inputs
  out_list$fit_func <- fit_func
  out_list$pred_func <- pred_func
  out_list$sim_func <- sim_func
  out_list$dens_func <- dens_func

  a <- 0
  ex_env = new.env()

  assign("a", a, envir = ex_env)
  out_list$geta<- function() get("a", envir = ex_env)
  out_list$puta <- function(.x) assign("a", .x, envir = ex_env)

  # Fit function
  out_list$fit <- function(expert, data) {
    mod <- fit_func(data)
    expert$fitted_model <- mod
    return(expert)
  }

  # Predict function, if no custom pred_func defined will use default on object
  out_list$predict <- function(expert, newdata) {
    if (is.null(pred_func)) {
      predict(expert$fitted_model, newdata)
    } else {
      pred_func(expert$fitted_model, newdata)
    }
  }

  # Simulate response from given data
  out_list$simulate <- function(expert, newdata) {
    if (is.null(sim_func)) {
      stop("Can only simulate if given sim_func.")
    }
    if (is.null(expert$fitted_model)) {
      stop("Requires fit to be called first")
    }
    sim_func(expert$fitted_model, newdata)
  }

  out_list$density <- function(expert, newdata)
    if (is.null(pred_func)) {
      predict(expert$fitted_model, newdata)
    } else {
      pred_func(expert$fitted_model, newdata)
    }

  return(structure(out_list, class = "expert"))
}

fit <- function(x, ...) {
  UseMethod("fit", x)
}

fit.expert <- function(expert, data) {
  expert$puta(1)
  #mod <- fit_func(data)
  #expert$fitted_model <- mod
}












