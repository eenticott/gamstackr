# Example setup ----------------------------------------------------------------
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

# ------------------------------------------------------------------------------
create_expert <- function(fit_func, pred_func = NULL, sim_func = NULL, dens_func = NULL) {
  out_list <- list()

  # Initialise from inputs
  out_list$fit_func <- fit_func
  out_list$pred_func <- pred_func
  out_list$sim_func <- sim_func
  out_list$dens_func <- dens_func

  # Fit function
  out_list$fit <- function(expert, data) {
    mod <- fit_func(data)
    expert$fitted_model <- mod
    return(expert)
  }

  # Predict function, if no custom pred_func defined will use default on object
  out_list$predict <- function(expert, newdata) {
    if (is.null(pred_func)) {
      predict(expert$fitted_model, newdata, type = "response")
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
      predict(expert$fitted_model, newdata, type = "response")
    } else {
      pred_func(expert$fitted_model, newdata)
    }

  return(structure(out_list, class = "expert"))
}

# ------------------------------------------------------------------------------

expert = create_expert(fit_func = fit_func)
expert$model_fit <- expert$fit(expert, data)


evaluate <- function(expert, windower, data) {
  windowed_dat <- windower(data)
  for (i in 1:length(attr(windowed_dat, "training_windows"))) {
    expert$fit(windowed_dat[attr(windowed_dat, "training_windows")[[i]],])
    expert$predict(windowed_dat[attr(windowed_dat, "testing_windows")[[i]],])
  }
}

inners[[1]] 1:3
inners[[2]] 4:5

# optional inner_idx <- list(c(1,4,5), c(2,3))

# experts taken in order and assigned

starts_fitted


evaluate(stack, window) {

}



create_stack(experts, inners) {
  if (!is.list(experts)) {
    stop("Experts must be supplied in a list.")
  }

  if (!is.list(experts[[0]])) {
    warning("Taking experts in order and assigning based on ")
  }
  stack <- list()
  # How should we define experts, as in dividing into correct outer/inner structure. Need an easy way to define
  if (length(experts) != length(inners)) {
    stop("Number of experts should match number of inner functions.")
  }
  for (ex in experts) {
    expert$fit(train)
    expert$predict(stack)
  }
}




  stack$fit(train, stack)
  stack$predict(newdata, type = "predictions") {
    predict(experts)
    predict(stack, experts)
  }


  create_




