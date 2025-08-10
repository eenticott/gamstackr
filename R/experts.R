# Core expert/stacking interface
# ------------------------------
# This file defines the user-facing expert wrapper, the evaluation helpers for
# rolling/expanding windows, and the stacker object used to learn windowed
# weights from expert outputs.
#
# Example setup ----------------------------------------------------------------
# fit_func <- function(data) {
#   lm(y ~ x1 + x2, data = data)
# }
#
# sim_func <- function(fitted_model, data) {
#   # returns simulated response from given data
#   mu <- predict(fitted_model, newdata = data)
#   sd <- sqrt(as.numeric(crossprod(fitted_model$residuals)/fitted_model$df.residuals))
#   rnorm(n = nrow(data), mu, sd = sd)
# }
#
# dens_func <- function(fitted_model, data) {
#   # returns likelihood of each data point given model
#   mu <- predict(fitted_model, newdata = data)
#   sd <- sqrt(as.numeric(crossprod(fitted_model$residuals)/fitted_model$df.residuals))
#   dnorm(n = nrow(data), mu, sd = sd, log = TRUE)
# }

# ------------------------------------------------------------------------------
#' Create an expert wrapper for model fitting, prediction, simulation and density evaluation
#'
#' Wraps user-supplied functions into an `expert` object with methods to fit a
#' model, obtain predictions, simulate responses, and evaluate predictive densities.
#'
#' @param fit_func function(data) -> model. Function that fits a model given a
#'   data frame and returns the fitted model object.
#' @param pred_func optional function(model, newdata) -> numeric. Function that
#'   returns predictions given a fitted model and new data. If `NULL`,
#'   `stats::predict(model, newdata, type = "response")` is used.
#' @param sim_func optional function(model, newdata) -> numeric. Function that
#'   simulates responses from the fitted model for given new data. Required if
#'   you plan to call `$simulate` on the returned expert.
#' @param dens_func optional function(model, newdata) -> numeric. Function that
#'   returns (log-)densities under the fitted model at each row of `newdata`.
#'   Required for density stacking (`type = "dens"`). Should typically return
#'   log-densities for numerical stability.
#'
#' @return An object of class `"expert"` (a list) with the following elements/methods:
#'   \describe{
#'     \item{$fit}{function(expert, data) -> expert. Fits the model and stores it in `$fitted_model`.}
#'     \item{$predict}{function(expert, newdata) -> numeric. Returns predictions via `pred_func` or `predict()`.}
#'     \item{$simulate}{function(expert, newdata) -> numeric. Uses `sim_func`; errors if not provided.}
#'     \item{$density}{function(expert, newdata) -> numeric. Uses `dens_func`; required for density stacking.}
#'     \item{$fitted_model}{The last fitted model (populated after `$fit`).}
#'   }
#'
#' @seealso [evaluate_expert()], [create_stacker()], [evaluate_stack()].
#' @export
#'
#' @examples
#' \dontrun{
#' # Simple linear-model expert
#' fit_func <- function(data) lm(y ~ x, data = data)
#' dens_func <- function(mod, newdata) {
#'   mu <- predict(mod, newdata)
#'   sig2 <- as.numeric(crossprod(mod$residuals) / mod$df.residual)
#'   stats::dnorm(newdata$y, mu, sqrt(sig2), log = TRUE)
#' }
#' ex <- create_expert(fit_func = fit_func, dens_func = dens_func)
#' ex <- ex$fit(ex, data.frame(y = rnorm(10), x = rnorm(10)))
#' preds <- ex$predict(ex, data.frame(x = rnorm(5)))
#' }
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

  # Return an object carrying the methods defined above
  return(structure(out_list, class = "expert"))
}

# ------------------------------------------------------------------------------

#' Evaluate an expert (or list of experts) across rolling/expanding windows
#'
#' Applies an `expert` across the windows produced by a `windower` (created by
#' [create_windower()]). When `expert` is a single `expert`, returns lists of
#' per-window outputs. When `expert` is a list of experts, evaluates each and
#' returns a list of outputs (one per expert).
#'
#' @param expert An `expert` created by [create_expert()], or a list of such
#'   experts.
#' @param windower A windowing function created by [create_windower()].
#' @param data A data frame containing response and features used for windowing
#'   and evaluation.
#' @param type Character vector of outputs to return. Any of
#'   `c("predict", "density", "model")`.
#'
#' @return If `expert` is a single expert, a named list with elements:
#'   \describe{
#'     \item{preds}{List of numeric vectors, one per test window (present if requested).}
#'     \item{dens}{List of numeric vectors, one per test window (present if requested).}
#'     \item{model}{The final fitted model object (present if requested).}
#'   }
#'   If `expert` is a list of experts, returns a list whose elements are the above
#'   structure for each expert.
#'
#' @seealso [create_expert()], [create_windower()], [evaluate_stack()].
#' @export
#'
#' @examples
#' \dontrun{
#' windower <- create_windower(initial_size = 100, horizon_size = 20, step_size = 20)
#' out <- evaluate_expert(expert, windower, df, type = c("predict", "density"))
#' }
evaluate_expert <- function(expert, windower, data, type = "predict") {
  list_of_dens <- NULL
  list_of_preds <- NULL

  # If a list of experts is supplied, evaluate each one recursively
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
  # Expand the data into train/test windows using the windower
  windowed_dat <- windower(data)
  list_of_preds <- list()
  model <- NULL
  pb <- txtProgressBar(min = 0,
                       max = length(attr(windowed_dat, "training_windows")),
                       style = 3,
                       width = 50,
                       char = "=")

  # Iterate windows: fit on train indices, compute requested outputs on test
  for (i in seq_along(attr(windowed_dat, "training_windows"))) {
    #cat(paste0(i, "/", length(attr(windowed_dat, "training_windows"))), "\n")
    expert <- expert$fit(expert, windowed_dat[attr(windowed_dat, "training_windows")[[i]],])
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

#' Fit and evaluate a stack over rolling/expanding windows
#'
#' Fits the stacking model on each training window and returns predicted
#' weights on the corresponding test window.
#'
#' @param stacker A stacker object created by [create_stacker()].
#' @param formula Model formula for the weights (passed to [mgcv::gam()]).
#'   Include any covariates in `stack_data` you wish to use to modulate weights.
#' @param windower A windower function created by [create_windower()].
#' @param stack_data Data frame used for stacking; must contain all columns
#'   referenced in `formula`.
#' @param expert_outputs Numeric N x K matrix of expert outputs aligned to
#'   `stack_data` rows (K = number of experts). For `type = "dens"`, supply
#'   predictive (log-)densities; for `type = "loss"`, supply point
#'   predictions. If `NULL`, outputs are computed from `stacker$experts` using
#'   the stacker's `type`.
#' @param RidgePen Non-negative scalar ridge penalty applied within the stacking
#'   families. Default is `1e-5`.
#'
#' @return A list of length equal to the number of windows. Each element is an
#'   `n_test x K` matrix of predicted weights for the test portion of that
#'   window.
#'
#' @seealso [create_stacker()], [create_windower()], [DensStack()], [LossStack()].
#' @export
#'
#' @examples
#' \dontrun{
#' windower <- create_windower(initial_size = 100, horizon_size = 20, step_size = 20)
#' W_list <- evaluate_stack(stacker, y ~ x, windower, stack_data = df, expert_outputs = preds)
#' }
evaluate_stack <- function(stacker, formula, windower, stack_data, expert_outputs, RidgePen = 1e-5) {
  # Optionally derive expert outputs if not provided
  if (is.null(expert_outputs)) {
    # Get predictions from experts
    expert_outputs <- matrix(nrow = nrow(stack_data), ncol = length(stacker$experts))
    for (i in seq_along(stacker$experts)) {
      if (attr(stacker$experts[[i]], "type") == "loss") {
        expert_outputs[,i] <- stacker$experts[[i]]$predict(stacker$experts[[i]], stack_data)
      }
      if (attr(stacker$experts[[i]], "type") == "dens") {
        expert_outputs[,i] <- stacker$experts[[i]]$density(stacker$experts[[i]], stack_data)
      }
    }
  }
  
  if (nrow(expert_outputs) != nrow(stack_data)) {
    stop("Stack data must correspond to given predictions, found mismatch in nrow.")
  }

  # Window the stacking data for train/test splits
  stacking_dat <- windower(stack_data)
  out_list <- list()

  pb <- txtProgressBar(min = 0,
                       max = length(attr(stacking_dat, "training_windows")),
                       style = 3,
                       width = 50,
                       char = "=")

  # For each window: slice rows, fit stacking model on train, predict weights on test
  for (i in seq_along(attr(stacking_dat, "training_windows"))) {
    #cat(paste0(i, "/", length(attr(stacking_dat, "training_windows"))), "\n")
    stack_dat <-  stacking_dat[attr(stacking_dat, "training_windows")[[i]],,drop = FALSE]
    test_dat <- stacking_dat[attr(stacking_dat, "testing_windows")[[i]],, drop = FALSE]
    train_ex_out <- expert_outputs[attr(stacking_dat, "training_windows")[[i]],,drop = FALSE]
    stacker <- stacker$fit_stack(stacker, formula, stack_dat, train_ex_out, RidgePen = RidgePen)
    out_list[[i]] <- stacker$predict(stacker, test_dat, "weights")
    setTxtProgressBar(pb, i)
  }
  close(pb)
  return(out_list)
}


#' Create a stacker object
#'
#' Constructs a stacker encapsulating expert models, a weight-generating
#' function, and the stacking type (density or loss). The returned object
#' provides methods to fit experts, fit the stacking model, and predict weights
#' (or stacked predictions) on new data.
#'
#' @param experts List of experts created by [create_expert()]. Length K equals
#'   the number of experts to be stacked.
#' @param weight_func Inner weight function (e.g. [ordinal()], [multinomial()],
#'   [MVN_weights()]) defining K expert weights; must satisfy
#'   `attr(weight_func, "num_weights") == length(experts)`.
#' @param type Either `"dens"` (stack predictive densities) or `"loss"` (stack
#'   point predictions under a loss).
#' @param loss Loss function used when `type = "loss"` (e.g. [square_loss()],
#'   [pinball_loss()]). Ignored when `type = "dens"`.
#'
#' @return A list with class-less methods/fields:
#'   \describe{
#'     \item{$experts}{The provided list of experts.}
#'     \item{$weight_func}{The inner weight function.}
#'     \item{$fit_experts}{function(stacker, train) -> stacker. Fits all experts to `train`.}
#'     \item{$fit_stack}{function(stacker, formula, stack, preds = NULL, RidgePen = 1e-5) -> stacker. Fits stacking model with [mgcv::gam()].}
#'     \item{$predict}{function(stacker, newdata, output = c("weights","preds")) -> matrix or numeric. Predicts weights (default) or stacked predictions.}
#'     \item{$fitted_stack}{The fitted [mgcv::gam] stacking model (after `$fit_stack`).}
#'     \item{$sp,$scale}{Smoothing/scale parameters stored between fits.}
#'   }
#'
#' @seealso [evaluate_stack()], [evaluate_expert()], [create_windower()],
#'   [DensStack()], [LossStack()].
#' @export
#'
#' @examples
#' \dontrun{
#' stk <- create_stacker(list(ex1, ex2), weight_func = multinomial(2), type = "loss", loss = square_loss)
#' }
create_stacker <- function(experts, weight_func, type = c("dens", "loss"), loss = NULL) {
  if (!is.list(experts)) {
    stop("Experts must be supplied in a list.")
  }

  # Check number of experts match number of weights
  if (length(experts) != attr(weight_func, "num_weights")) {
    stop("Number of experts should match number of inner functions.")
  }

  stacker <- list()
  stacker$experts <- experts
  stacker$weight_func <- weight_func
  stacker$experts_fitted <- FALSE
  stacker$coef <- NULL
  stacker$sp <- NULL
  stacker$scale <- NULL

  K <- length(experts)

  # fit experts to same data
  # Fit each expert in-place on a shared training set
  stacker$fit_experts <- function(stacker, train) {
    experts <- stacker$experts
    for (i in seq_len(K)) {
      for (j in seq_along(experts[[i]])) {
        ex <- experts[[i]][[j]]
        experts[[i]][[j]] <- ex$fit(ex, train)
      }
    }
    stacker$experts <- experts
    stacker$experts_fitted <- TRUE
    return(stacker)
  }

  # Fit the stacking model (GAM) using expert outputs
  stacker$fit_stack <- function(stacker, formula_list, stack, preds = NULL, RidgePen = 1e-5) {
    # preds - matrix of predictions or densities
    if (!(stacker$experts_fitted) && is.null(preds)) {
      stop("fit_stack requires experts to be fitted prior or densities to be directly supplied.")
    }

    # Derive expert outputs if not supplied
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

    # Select the appropriate stacking family
    if (type == "dens") {
      pre_fam <- DensStack(preds, weight_func, RidgePen)
    }
    if (type == "loss") {
      pre_fam <- LossStack(preds, loss, weight_func, RidgePen)
    }

    # Reuse smoothing and scale if available to speed refits
    if (!is.null(stacker$sp)) {
      fitted_stack <- mgcv::gam(formula_list, data = stack, family = pre_fam,
                          in.out=list(sp=stacker$sp, scale=stacker$scale))
    } else {
      fitted_stack <- mgcv::gam(formula_list, data = stack, family = pre_fam)
    }

    stacker$sp <- fitted_stack$sp
    stacker$scale <- fitted_stack$scale
    stacker$fitted_stack <- fitted_stack
    return(stacker)
  }

  # Predict either weights or stacked predictions on new data
  stacker$predict <- function(stacker, newdata, output="weights") {
    if (output == "weights") {
      out <- predict(stacker$fitted_stack, newdata, type = "response")
      return(out)
    }
    if (output == "preds") {
      preds <- matrix(nrow = nrow(newdata), ncol = K)
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