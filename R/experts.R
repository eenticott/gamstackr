
#' @title Create an Expert Object
#' @description
#' Create an expert that can be fitted, used for prediction, simulation, and density evaluation.
#'
#' @param fit_func Function that returns a fitted model given a dataframe.
#' @param pred_func Function that returns predictions given fitted model and data. Optional.
#' @param sim_func Function that simulates from the estimated density given model and data. Optional.
#' @param dens_func Function that returns density under fitted model at given data points. Optional.
#'
#' @return An expert object (S3 class 'expert') for use in stacking workflows.
#' @export
#'
#' @examples
#' fit_func <- function(data) lm(Sepal.Length ~ Sepal.Width, data = data)
#' expert <- create_expert(fit_func)
create_expert <- function(fit_func, pred_func = NULL, sim_func = NULL, dens_func = NULL) {
  out_list <- list()

  # Store user-supplied functions
  out_list$fit_func <- fit_func
  out_list$pred_func <- pred_func
  out_list$sim_func <- sim_func
  out_list$dens_func <- dens_func

  # Fit the expert to data
  out_list$fit <- function(expert, data) {
    mod <- fit_func(data)
    expert$fitted_model <- mod
    return(expert)
  }

  # Predict using the expert
  out_list$predict <- function(expert, newdata) {
    if (is.null(pred_func)) {
      predict(expert$fitted_model, newdata, type = "response")
    } else {
      pred_func(expert$fitted_model, newdata)
    }
  }

  # Simulate from the expert's fitted distribution
  out_list$simulate <- function(expert, newdata) {
    if (is.null(sim_func)) {
      stop("Can only simulate if given sim_func.")
    }
    if (is.null(expert$fitted_model)) {
      stop("Requires fit to be called first")
    }
    sim_func(expert$fitted_model, newdata)
  }

  # Evaluate density at new data
  out_list$density <- function(expert, newdata) {
    dens_func(expert$fitted_model, newdata)
  }

  return(structure(out_list, class = "expert"))
}


#' @title Evaluate an Expert Object Over Rolling Windows
#' @description
#' Evaluates an expert object over a set of rolling or expanding windows, optionally returning predictions, densities, and fitted models.
#'
#' @param expert An expert object created by `create_expert`, or a list of such objects.
#' @param windower A window object created by `create_windower`.
#' @param data A dataframe containing response and features.
#' @param type Character vector: which outputs to return from "predict", "density", and/or "model".
#'
#' @return A list containing the outputs determined by `type`.
#' @export
#'
#' @examples
#' # See package vignette for usage
evaluate_expert <- function(expert, windower, data, type = "predict") {
  list_of_dens <- NULL
  list_of_preds <- NULL

  # If a list of experts, evaluate each recursively
  if (!inherits(expert, "expert")) {
    out_list <- list()
    i <- 1
    for (ex in expert) {
      out_list[[i]] <- evaluate_expert(ex, windower, data, type)
      cat(paste("----- Expert", i, "completed -----", "\n"))
      i <- i + 1
    }
    return(out_list)
  }
  windowed_dat <- windower(data)
  list_of_preds <- list()
  model <- NULL
  pb <- txtProgressBar(min = 0,
                       max = length(attr(windowed_dat, "training_windows")),
                       style = 3,
                       width = 50,
                       char = "=")

  for (i in 1:length(attr(windowed_dat, "training_windows"))) {
    # Fit expert to training window
    expert <- expert$fit(expert, windowed_dat[attr(windowed_dat, "training_windows")[[i]],])
    # Optionally compute density and/or predictions
    if ("density" %in% type) {
      list_of_dens[[i]] <- expert$density(expert, windowed_dat[attr(windowed_dat, "testing_windows")[[i]],])
    }
    if ("predict" %in% type) {
      list_of_preds[[i]] <- expert$predict(expert, windowed_dat[attr(windowed_dat, "testing_windows")[[i]],])
    }
    if ("model" %in% type) {
      model <- expert$fitted_model
    }
    setTxtProgressBar(pb, i)
  }
  close(pb)
  return(list("preds" = list_of_preds, "dens" = list_of_dens, "model" = model))
}


#' @title Evaluate a Set of Experts with Weighting Functions
#' @description
#' Evaluates a stacker object over rolling windows, fitting and predicting weights for each window.
#'
#' @param stacker Stacker object created by `create_stacker`.
#' @param formula Model formula for the weights. For nested structure this should be a list of lists.
#' @param windower A windower object created by `create_windower`.
#' @param stack_data Data for stacking, should contain all columns specified by formula.
#' @param preds Matrix of predictions or densities (rows must match stack_data).
#' @param RidgePen Numeric, ridge penalty for stacking (default 1e-5).
#'
#' @return List of predicted weights from stacking process.
#' @export
#'
#' @examples
#' # See package vignette for usage
evaluate_stack <- function(stacker, formula, windower, stack_data, preds, RidgePen = 1e-5) {
  if (nrow(preds) != nrow(stack_data)) {
    stop("Stack data must correspond to given predictions, found mismatch in nrow.")
  }

  stacking_dat <- windower(stack_data)
  out_list <- list()

  pb <- txtProgressBar(min = 0,
                       max = length(attr(stacking_dat, "training_windows")),
                       style = 3,
                       width = 50,
                       char = "=")

  for (i in 1:length(attr(stacking_dat, "training_windows"))) {
    # Fit stacker to training window
    stack_dat <-  stacking_dat[attr(stacking_dat, "training_windows")[[i]],,drop = FALSE]
    test_dat <- stacking_dat[attr(stacking_dat, "testing_windows")[[i]],, drop = FALSE]
    clod <- preds[attr(stacking_dat, "training_windows")[[i]],,drop = FALSE]
    clotd <- preds[attr(stacking_dat, "testing_windows")[[i]],, drop = FALSE]
    stacker <- stacker$fit_stack(stacker, formula, stack_dat, clod, RidgePen = RidgePen)
    out_list[[i]] <- stacker$predict(stacker, test_dat, "weights")
    setTxtProgressBar(pb, i)
  }
  close(pb)
  return(out_list)
}



#' @title Create a Stacker Object from Experts and Weight Functions
#' @description
#' Creates a stacker object that can fit and predict expert weights using a specified weighting function and stacking type.
#'
#' @param experts List of expert objects created by `create_expert`.
#' @param weight_func Function that determines expert weights.
#' @param type Character: "dens" or "loss". Whether densities or predictions are being stacked.
#' @param loss Loss function to optimise with stacking. Only needed if type = "loss".
#'
#' @return Stacker object that can take data and return predicted weights for experts.
#' @export
#'
#' @examples
#' # See package vignette for usage
create_stacker <- function(experts, weight_func, type, loss = NULL) {
  if (!is.list(experts)) {
    stop("Experts must be supplied in a list.")
  }

  if (!is.list(experts[[1]])) {
    warning("Taking experts in order and assigning based on inner numbers.")
  }

  stacker <- list()
  stacker$experts <- experts
  stacker$weight_func <- weight_func
  stacker$experts_fitted <- FALSE
  stacker$coef <- NULL
  stacker$sp <- NULL
  stacker$scale <- NULL

  # Check number of experts matches number of weights
  if (length(experts) != attr(weight_func, "num_weights")) {
    stop("Number of experts should match number of inner functions.")
  }

  K <- length(experts)

  # Fit all experts to the same data
  stacker$fit_experts <- function(stacker, train) {
    experts <- stacker$experts
    for (i in 1:K) {
      for (j in 1:length(experts[[i]])) {
        ex <- experts[[i]][[j]]
        experts[[i]][[j]] <- ex$fit(ex, train)
      }
    }
    stacker$experts <- experts
    stacker$experts_fitted <- TRUE
    return(stacker)
  }

  # Fit stack using densities or predictions from fitted experts
  stacker$fit_stack <- function(stacker, formula_list, stack, preds = NULL, RidgePen = 1e-5) {
    # preds - matrix of predictions or densities
    if (!(stacker$experts_fitted) & is.null(preds)) {
      stop("fit_stack requires experts to be fitted prior or densities to be directly supplied.")
    }

    if (is.null(preds)) {
      preds <- matrix(nrow = nrow(stack), ncol = K)
      experts <- stacker$experts
      for (k in 1:K) {
        expert <- experts[[k]]
        if (type == "dens") {
          preds[,k] <- expert$density(expert, stack)
        }
        if (type == "loss") {
          preds[,k] <- expert$predict(expert, stack)
        }
      }
    }

    if (type == "dens") {
      pre_fam <- DensStack(preds, weight_func, RidgePen)
    }
    if (type == "loss") {
      pre_fam <- LossStack(preds, loss, weight_func, RidgePen)
    }

    if (!is.null(stacker$sp)) {
      fitted_stack <- gam(formula_list, data = stack, family = pre_fam,
                          in.out=list(sp=stacker$sp, scale=stacker$scale))
    } else {
      fitted_stack <- gam(formula_list, data = stack, family = pre_fam)
    }

    stacker$sp <- fitted_stack$sp
    stacker$scale <- fitted_stack$scale
    stacker$fitted_stack <- fitted_stack
    return(stacker)
  }

  # Predict weights or predictions for new data
  stacker$predict <- function(stacker, newdata, output="weights") {
    if (output == "weights") {
      out <- predict(stacker$fitted_stack, newdata, type = "response")
      return(out)
    }
    if (output == "preds") {
      preds <- matrix(nrow = nrow(stack), ncol = K)
      experts <- stacker$experts
      for (k in 1:K) {
        expert <- experts[[k]]
        if (type == "dens") {
          preds[,k] <- expert$density(expert, newdata)
        }
        if (type == "loss") {
          preds[,k] <- expert$predict(expert, newdata)
        }
      }
      out <- rowSums(predict(stacker$fitted_stack, newdata, type = "response") * preds)
      return(out)
    }
  }
  return(stacker)
}



#' @title Fit a Set of Experts and Stack Them Over an Out-of-Sample Window
#' @description
#' Fits experts and stacker objects over rolling windows for out-of-sample evaluation.
#'
#' @param experts List of expert objects.
#' @param stacker Stacker object created by `create_stacker`.
#' @param formula Model formula for the weights. For nested structure this should be a list of lists.
#' @param expert_windower A windower object created by `create_windower` used for evaluating experts.
#' @param stack_windower A windower object created by `create_windower` used for evaluating stacks.
#' @param dat Data for fitting experts and stacking, should contain all columns specified by formula.
#'
#' @return List of predicted weights from stacking process.
#' @export
#'
#' @examples
#' # See package vignette for usage
evaluate_all <- function(experts, stacker, formula, expert_windower, stack_windower, dat) {
  # First evaluate experts over expert windows.
  expert_dat <- expert_windower(dat)
  stack_dat <- stack_windower(dat[])
  evaluate_expert()
  out_list <- list()
  for (i in 1:length(attr(stacking_dat, "training_windows"))) {
    cat(paste0(i, "/", length(attr(stacking_dat, "training_windows"))), "\n")
    stack_dat <-  stacking_dat[attr(stacking_dat, "training_windows")[[i]],,drop = FALSE]
    test_dat <- stacking_dat[attr(stacking_dat, "testing_windows")[[i]],, drop = FALSE]
    clod <- lapply(list_of_densities, function(x) x[attr(stacking_dat, "training_windows")[[i]],, drop = FALSE])
    clotd <- lapply(list_of_densities, function(x) x[attr(stacking_dat, "testing_windows")[[i]],, drop = FALSE])
    stacker <- stacker$fit_stack(stacker, formula, stack_dat, clod)
    out_list[[i]] <- stacker$predict(stacker, test_dat, clotd, "weights")
  }
  return(out_list)
}
