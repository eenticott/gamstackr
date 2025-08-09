## TODO: Allow fixed proportional window
## TODO: Don't require all args (have control = named list?)

#' Create train/test rolling windows for model evaluation
#'
#' Generates rolling or expanding windows for time series or cross-validation.
#'
#' @param df Data frame of data to evaluate.
#' @param initial_size Integer. Size of the initial training window.
#' @param window_size Integer. Width of the sliding window (ignored for expanding windows).
#' @param horizon_size Integer. Number of values to forecast ahead (test window size).
#' @param step_size Integer. Number of rows added in each iteration (default: horizon_size).
#' @param type Character. Either "expanding" or "sliding".
#'
#' @return Data frame with added attributes: 'training_windows' and 'testing_windows', each a list of row indices.
#' @export
#'
#' @examples
#' df <- data.frame(x = 1:100)
#' win <- create_windows(df, initial_size = 50, window_size = 20, horizon_size = 10)
#' attr(win, "training_windows")
#' attr(win, "testing_windows")
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
  # Add final window
  if (cutoff + horizon_size - step_size < N) {
    remaining = N - (cutoff + horizon_size - step_size)
    if (type == "expanding") {
      train_idx <- 1:cutoff
      test_idx <- (cutoff + 1):(cutoff + remaining)
    } else {
      train_idx <- (cutoff - window_size + 1):(cutoff)
    }
    test_idx <- (cutoff + 1):(cutoff + remaining)

    training_idx[[i]] <- train_idx
    testing_idx[[i]] <- test_idx
  }
  return(structure(df, "training_windows" = training_idx, "testing_windows" = testing_idx))
}

##' Plot windows for easy visualisation.
#'
#' Visualizes the training and testing windows produced by `create_windows`.
#'
#' @param windowed_dat Output of `create_windows`.
#'
#' @return Invisibly returns NULL. Displays a plot showing windowed data.
#' @export
#'
#' @examples
#' df <- data.frame(x = 1:100)
#' win <- create_windows(df, initial_size = 50, window_size = 20, horizon_size = 10)
#' plot_windows(win)
plot_windows <- function(windowed_dat) {
  train_idxs <- attr(windowed_dat, "training_windows")
  test_idxs <- attr(windowed_dat, "testing_windows")

  plot(1, type="n", xlab="Index", ylab="Iteration", xlim=c(1, nrow(windowed_dat)), ylim=c(1, length(train_idxs)))

  for (i in seq_along(train_idxs)) {
    lines(x = train_idxs[[i]], y = rep(i, length(train_idxs[[i]])), col = "blue", lwd = 2)
    lines(x = test_idxs[[i]], y = rep(i, length(test_idxs[[i]])), col = "orange", lwd = 2)
  }

  legend("bottomright", legend = c("Train", "Test"), col = c("blue", "orange"), lty = c(1,1))
}

##' Create a windower for use in stacking/experts.
#'
#' Returns a function that applies `create_windows` with preset parameters.
#'
#' @param initial_size Integer. Size of the initial training window.
#' @param horizon_size Integer. Number of values to forecast ahead (test window size).
#' @param window_size Integer. Width of the sliding window (optional).
#' @param step_size Integer. Number of rows added in each iteration (optional).
#' @param type Character. Either "expanding" or "sliding".
#'
#' @return A function that takes a data frame and returns windowed data with attributes.
#' @export
#'
#' @examples
#' windower <- create_windower(50, 10, 20)
#' win <- windower(data.frame(x = 1:100))
create_windower <- function(initial_size, horizon_size, window_size = NULL, step_size = NULL, type = "expanding") {
  windower <- function(df) {
    create_windows(df, initial_size, window_size, horizon_size, step_size, type)
  }
  return(windower)
}


#' Split a data frame by a factor column
#'
#' Splits a data frame into a list of data frames, one for each unique value of a specified factor column.
#'
#' @param df Data frame to split.
#' @param fact_split Character or integer. Name or index of the column to split by.
#'
#' @return List of data frames, one per unique value in the specified column.
#' @export
#'
#' @examples
#' df <- data.frame(a = rep(1:2, each = 3), b = 1:6)
#' split_by(df, "a")
split_by <- function(df, fact_split) {
  out_df <- list()
  i <- 1
  for (fs in unique(df[, fact_split])) {
    out_df[[i]] <- df[df[, fact_split] == fs, ]
    i <- i + 1
  }
  return(out_df)
}






