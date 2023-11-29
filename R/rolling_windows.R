# TODO: Allow fixed proportional window
# TODO: Don't require all args (have control = named list?)

#' Create train/test rolling windows for model evaluation
#'
#' @param df dataframe of data you want to evaluate
#' @param initial_size Int, How large the starting window is.
#' @param window_size Int, Width of sliding window.
#' @param horizon_size Int,  Number of values to be forecast ahead.
#' @param step_size Int, Number of rows added in each iteration.
#' @param type str, "expanding" or "sliding"
#'
#' @return df with added attributes for training and testing indexes
#' @export
#'
#' @examples
create_windows <- function(df, initial_size, window_size, horizon_size, step_size = NULL, type = "expanding") {
  # Input checks
  if (typeof(type) != "character") {stop("type must be a string, expanding or sliding allowed.")}
  type <- tolower(type)
  if (!(type %in% c("expanding", "sliding"))) {
    stop("Type not recognised.")
  }
  N <- nrow(df)
  if (N < initial_size) {
    stop("Initial window size exceeds given data.")
  }
  if (is.null(step_size)) {
    step_size = horizon_size
  }
  cutoff = initial_size
  train_start_idx = 1:initial_size
  test_start_idx <- (cutoff + 1):(cutoff + horizon_size)
  training_idx <- list()
  testing_idx <- list()
  training_idx[[1]] <- train_start_idx
  testing_idx[[1]] <- test_start_idx
  i <- 2
  cutoff <- cutoff + step_size
  while ((cutoff+horizon_size) <= N) {
    if (type == "expanding") {
      train_idx <- 1:cutoff
      test_idx <- (cutoff + 1):(cutoff + horizon_size)
    } else {
      train_idx <- (cutoff - window_size + 1):(cutoff)
    }
    test_idx <- (cutoff + 1):(cutoff + horizon_size)

    training_idx[[i]] <- train_idx
    testing_idx[[i]] <- test_idx
    i <- i + 1
    cutoff <- cutoff + step_size
  }
  return(structure(df, "training_windows" = training_idx, "testing_windows" = testing_idx))
}

#' Plot windows for easy visualisation.
#'
#' @param windowed_dat Output of create_windows.
#'
#' @return Displays plot showing windowed data
#' @export
#'
#' @examples
plot_windows <- function(windowed_dat) {
  train_idxs <- attr(windowed_dat, "training_windows")
  test_idxs <- attr(windowed_dat, "testing_windows")

  plot(1, type="n", xlab="Index", ylab="Iteration", xlim=c(1, nrow(windowed_dat)), ylim=c(1, length(train_idxs)))

  for (i in 1:length(train_idxs)) {
    lines(x = train_idxs[[i]], y = rep(i, length(train_idxs[[i]])), col = "blue", lwd = 2)
    lines(x = test_idxs[[i]], y = rep(i, length(test_idxs[[i]])), col = "orange", lwd = 2)
  }

  legend("bottomright", legend = c("Train", "Test"), col = c("blue", "orange"), lty = c(1,1))
}

#' Create a windower for use in stacking/experts.
#' @param initial_size Int, How large the starting window is.
#' @param horizon_size Int,  Number of values to be forecast ahead.
#' @param window_size Int, Width of sliding window.
#' @param step_size Int, Number of rows added in each iteration.
#' @param type str, "expanding" or "sliding"
#'
#' @return Windower object.
#' @export
#'
#' @examples
create_windower <- function(initial_size, horizon_size, window_size = NULL, step_size = NULL, type = "expanding") {
  windower <- function(df) {
    create_windows(df, initial_size, window_size, horizon_size, step_size, type)
  }
  return(windower)
}

split_by <- function(df, fact_split) {
   N <- length(unique(df[,fact_split]))
   out_df <- list()
   i <- 1
   for (fs in unique(df[,fact_split])) {
     out_df[[i]] <- df[df[,fact_split] == fs]
   }
   return(out_df)
}






