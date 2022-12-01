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

  out_list$density <- function(expert, newdata) {
      dens_func(expert$fitted_model, newdata)
  }

  return(structure(out_list, class = "expert"))
}

# ------------------------------------------------------------------------------

expert = create_expert(fit_func = fit_func)
expert$model_fit <- expert$fit(expert, data)


evaluate_expert <- function(expert, windower, data, type = "predict") {
  list_of_dens <- NULL
  list_of_preds <- NULL

  if (!inherits(expert, "expert")) {
    out_list <- list()
    i <- 1
    for (ex in expert) {
      out_list[[i]] <- evaluate_expert(ex, windower, data, type)
      i <- i + 1
    }
    return(out_list)
  }
  windowed_dat <- windower(data)
  list_of_preds <- list()
  for (i in 1:length(attr(windowed_dat, "training_windows"))) {
    expert <- expert$fit(expert, windowed_dat[attr(windowed_dat, "training_windows")[[i]],])
    if ("density" %in% type) {
      list_of_dens[[i]] <- expert$density(expert, windowed_dat[attr(windowed_dat, "testing_windows")[[i]],])
    }
    if ("predict" %in% type) {
      list_of_preds[[i]] <- expert$predict(expert, windowed_dat[attr(windowed_dat, "testing_windows")[[i]],])
    }
  }
  return(list("preds" = list_of_preds, "dens" = list_of_dens))
}

# force horizon_size = step_size for density evaluation to ensure we can bind properly

evaluate_stack <- function(stacker, windower1, windower2, data) {
  training_dat <- windower1(data)
  stacking_dat <- windower2(data)
  list_of_preds <- list()
  for (i in 1:length(attr(windowed_dat, "training_windows"))) {
    stacker <- stacker$fit_experts(stacker, training_dat[attr(training_dat, "training_windows")[[i]],])
    stacker <- stacker$fit_stack(stacker, stacking_dat[attr(stacking_dat, "training_windows")[[i]],])
    list_of_preds[[i]] <- stacker$predict(stacker, stacking_dat[attr(stacking_dat, "testing_windows")[[i]],])
  }
}



create_stacker <- function(experts, inners) {
  if (!is.list(experts)) {
    stop("Experts must be supplied in a list.")
  }

  if (!is.list(experts[[0]])) {
    warning("Taking experts in order and assigning based on inner numbers.")
  }

  stacker <- list()
  stacker$experts <- experts
  stacker$experts_fitted <- FALSE
  # How should we define experts, as in dividing into correct outer/inner structure. Need an easy way to define
  if (length(experts) != length(inners)) {
    stop("Number of experts should match number of inner functions.")
  }

  K <- length(experts)

  # fit experts to same data
  stacker$fit_experts <- function(stacker, train) {
    stacker$experts_fitted <- TRUE
    experts <- stacker$experts
    for (i in 1:K) {
      for (j in 1:length(experts[[i]])) {
        ex <- experts[[i]][[j]]
        experts[[i]][[j]] <- ex$fit(ex, train)
      }
    }
    stacker$experts <- experts
    return(stacker)
  }

  # fit stack using densities from fitted experts
  stacker$fit_stack <- function(stacker,formula_list, stack) {
    if (!(stacker$experts_fitted)) {
      stop("fit_stack requires experts to be fitted prior")
    }

    list_of_densities <- list()
    experts <- stacker$experts
    list_of_densities <- list()
    for (i in 1:K) {
      dens <- matrix(nrow = nrow(stack), ncol = length(experts[[i]]))
      for (j in 1:length(experts[[i]][[j]])) {
        dens[,j] <- experts[[i]][[j]]$density(experts[[i]][[j]], stack)
      }
      list_of_densities[[i]] <- dens
    }
    stacker$list_of_densities <- list_of_densities

    pre_fam <- NestedStack(list_of_densities, inners, RidgePen = 1e-05)
    fitted_stack <- gam(formula_list, data = stack, family = pre_fam)
    stacker$fitted_stack <- fitted_stack
    return(stacker)
  }

  stacker$predict <- function(stacker, stack_data, list_of_densities, type) {
    if (type == "weights") {
      predict(stacker$fitted_stack, newdata, type = "response")
    }
    if (type == "density") {
      rowSums(predict(stacker$fitted_stack, newdata, type = "response") * do.call("cbind", list_of_densities))
    }
  }
}






