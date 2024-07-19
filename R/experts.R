# TODO: Redesign, interface is confusing. Could put into seperate package so 'gamstackr' focus remains clear.

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
#' Create an expert that can be fitted, retrieve predictions and simulate from fitted distribution.
#'
#' @param fit_func A function that returned a fitted model given a dataframe.
#' @param pred_func A function that returns predictions given fitted model and data.
#' @param sim_func A function that simulated from the estimated density given model and data.
#' @param dens_func A function returns density under fitted model at given data points.
#'
#' @return An expert object that can be used within evaluate_expert and evaluate_stack.
#' @export
#'
#' @examples
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

#' Evaluate an expert object over a sliding or expanding window.
#'
#' @param expert An expert object created by `create_expert`.
#' @param windower A window object created by `create_windower`.
#' @param data A dataframe containing response and features.
#' @param type List the options you want outputted from "predict", "density" and "model".
#'
#' @return A list containing the outputs determined by type.
#' @export
#'
#' @examples
evaluate_expert <- function(expert, windower, data, type = "predict") {
  list_of_dens <- NULL
  list_of_preds <- NULL

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

#' Evaluate a set of experts with weighting functions.
#'
#' @param stacker Stacker object created by `create_stacker`.
#' @param formula Model formula for the weights. For nested stucture this should be a list of lists.
#' @param windower A windower object created by `create_windower`.
#' @param stack_data Data that should be used for stacking, should contain all columns specified by formula.
#' @param list_of_densities For nested stack this is a list of lists. Densities can be retrieved from experts. Rows should match that of `stack_data`.
#'
#' @return list of predicted weights from stacking process.
#' @export
#'
#' @examples
evaluate_stack <- function(stacker, formula, windower, stack_data, list_of_densities) {
  if (nrow(list_of_densities[[1]]) != nrow(stack_data)) {
    stop("Stack data must correspond to given densities, found mismatch in nrow.")
  }

  stacking_dat <- windower(stack_data)
  out_list <- list()

  pb <- txtProgressBar(min = 0,
                       max = length(attr(stacking_dat, "training_windows")),
                       style = 3,
                       width = 50,
                       char = "=")

  for (i in 1:length(attr(stacking_dat, "training_windows"))) {
    #cat(paste0(i, "/", length(attr(stacking_dat, "training_windows"))), "\n")
    stack_dat <-  stacking_dat[attr(stacking_dat, "training_windows")[[i]],,drop = FALSE]
    test_dat <- stacking_dat[attr(stacking_dat, "testing_windows")[[i]],, drop = FALSE]
    clod <- lapply(list_of_densities, function(x) x[attr(stacking_dat, "training_windows")[[i]],, drop = FALSE])
    clotd <- lapply(list_of_densities, function(x) x[attr(stacking_dat, "testing_windows")[[i]],, drop = FALSE])
    stacker <- stacker$fit_stack(stacker, formula, stack_dat, clod)
    out_list[[i]] <- stacker$predict(stacker, test_dat, clotd, "weights")
    setTxtProgressBar(pb, i)
  }
  close(pb)
  return(out_list)
}


#' Create a stacking object from experts and inner functions.
#'
#' @param experts List of expert objects created by `create_expert`.
#' @param inners Function that determines how weights
#'
#' @return Stacking object that can take data and return predicted weights for experts.
#' @export
#'
#' @examples
create_stacker <- function(experts, inners) {
  if (!is.list(experts)) {
    stop("Experts must be supplied in a list.")
  }

  if (!is.list(experts[[1]])) {
    warning("Taking experts in order and assigning based on inner numbers.")
  }

  stacker <- list()
  stacker$experts <- experts
  stacker$experts_fitted <- FALSE
  stacker$coef <- NULL
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
  stacker$fit_stack <- function(stacker,formula_list, stack, list_of_densities = NULL) {
    if (!(stacker$experts_fitted) & is.null(list_of_densities)) {
      stop("fit_stack requires experts to be fitted prior or densities to be directly supplied.")
    }

    if (is.null(list_of_densities)) {
      list_of_densities <- list()
      experts <- stacker$experts
      for (i in 1:K) {
        dens <- matrix(nrow = nrow(stack), ncol = length(experts[[i]]))
        for (j in 1:length(experts[[i]])) {
          dens[,j] <- experts[[i]][[j]]$density(experts[[i]][[j]], stack)
        }
        list_of_densities[[i]] <- dens
      }
      stacker$list_of_densities <- list_of_densities
    } else {
      stacker$list_of_densities <- list_of_densities
    }

    pre_fam <- NestedStack(list_of_densities, inners, RidgePen = 1e-05)
    fitted_stack <- gam(formula_list, data = stack, family = pre_fam, start = stacker$coef)
    stacker$coef <- fitted_stack$coef
    stacker$fitted_stack <- fitted_stack
    return(stacker)
  }

  stacker$predict <- function(stacker, newdata, list_of_densities, type) {
    stacker$fitted_stack$family$putMWF(TRUE)
    if (type == "weights") {
      out <- predict(stacker$fitted_stack, newdata, type = "response")
      stacker$fitted_stack$family$putMWF(FALSE)
      return(out)
    }
    if (type == "density") {
      out <- rowSums(predict(stacker$fitted_stack, newdata, type = "response") * do.call("cbind", list_of_densities))
      stacker$fitted_stack$family$putMWF(FALSE)
      return(out)
    }
  }

  return(stacker)
}


#' Fit a set of experts and stack them over an out-of-sample window.
#'
#' @param experts description
#' @param stacker Stacker object created by `create_stacker`.
#' @param formula Model formula for the weights. For nested stucture this should be a list of lists.
#' @param expert_windower A windower object created by `create_windower` used for evaluating experts.
#' @param stack_windower A windower object created by `create_windower` used for evaluating stacks.
#' @param dat Data that should be used for fitting experts and stacking, should contain all columns specified by formula.
#'
#' @return list of predicted weights from stacking process.
#' @export
#'
#' @examples
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

