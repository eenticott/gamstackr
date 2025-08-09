# Suppress R CMD check notes for use of .data in tidy evaluation
if (getRversion() >= "2.15.1") utils::globalVariables(c(".data"))
# list utils

##' Subtract corresponding elements of two lists
##'
##' @param l1 A list of numeric vectors or matrices
##' @param l2 A list of numeric vectors or matrices (same length as l1)
##' @return A list with each element as l1[[i]] - l2[[i]]
##' @examples
##' list_take_list(list(1:3, 4:6), list(1:3, 1:3))
list_take_list <- function(l1, l2) {
  out <- list()
  if (length(l1) != length(l2)) {
    stop("List length didn't match")
  }
  for (i in seq_along(l1)) {
    out[[i]] <- l1[[i]] - l2[[i]]
  }
  out
}

##' Multiply each element of a list by a vector
##'
##' @param l A list of numeric vectors or matrices
##' @param v A numeric vector (same length as l)
##' @return A list with each element as l[[i]] * v[i]
##' @examples
##' list_by_vector(list(1:2, 3:4), c(2, 3))
list_by_vector <- function(l, v) {
  out <- list()
  for (i in seq_along(v)) {
    out[[i]] <- l[[i]] * v[i]
  }
  return(out)
}

##' Multiply or cross-product corresponding elements of two lists
##'
##' @param list1 A list of numeric vectors or matrices
##' @param list2 A list of numeric vectors or matrices
##' @param elementwise Logical; if TRUE, multiply elementwise, else cross-product
##' @return A list of products or cross-products
##' @examples
##' list_by_list(list(1:2, 3:4), list(2:3, 4:5), elementwise=TRUE)
list_by_list <- function(list1, list2, elementwise = FALSE) {
  if (!is.list(list1)) {
    list1 <- list(list1)
  }
  if (!is.list(list2)) {
    list2 <- list(list2)
  }
  n1 <- length(list1)
  n2 <- length(list2)
  if (n1 != n2) {
    errorCondition("Lists are not of same length")
  }
  out <- list()
  for (i in seq_len(n1)) {
    if (elementwise) {
      out[[i]] <- list1[[i]] * list2[[i]]
    } else {
      out[[i]] <- t(list1[[i]]) %*% list2[[i]]
    }
  }
  return(out)
}

##' Elementwise multiply two lists
##'
##' @param list1 A list of numeric vectors or matrices
##' @param list2 A list of numeric vectors or matrices (same length as list1)
##' @return A list with each element as list1[[i]] * list2[[i]]
##' @examples
##' list_times_list(list(1:2, 3:4), list(2:3, 4:5))
list_times_list <- function(list1, list2) {
  out <- list()
  for (i in seq_along(list1)) {
    out[[i]] <- list1[[i]] * list2[[i]]
  }
  out
}

##' Repeat a matrix vertically and horizontally
##'
##' @param mat A matrix
##' @param v Integer, number of vertical repeats
##' @param h Integer, number of horizontal repeats
##' @return A matrix repeated v times vertically and h times horizontally
##' @examples
##' repmat(matrix(1:4, 2, 2), 2, 2)
repmat <- function(mat, v = 1, h = 1) {
  kronecker(matrix(1, h, v), mat)
}

# Matrix utils
# Found in https://stackoverflow.com/questions/24299171/function-to-split-a-matrix-into-sub-matrices-in-r
##' Split a matrix into submatrices
##'
##' @param M A matrix
##' @param r Number of rows per submatrix
##' @param c Number of columns per submatrix
##' @return An array of submatrices
##' @examples
##' mat_split(matrix(1:16, 4, 4), 2, 2)
mat_split <- function(M, r, c) {
  nr <- ceiling(nrow(M) / r)
  nc <- ceiling(ncol(M) / c)
  newM <- matrix(NA, nr * r, nc * c)
  newM[seq_len(nrow(M)), seq_len(ncol(M))] <- M
  div_k <- kronecker(matrix(seq_len(nr * nc), nr, byrow = TRUE), matrix(1, r, c))
  matlist <- split(newM, div_k)
  N <- length(matlist)
  mats <- unlist(matlist)
  dim(mats) <- c(r, c, N)
  return(mats)
}

##' Cumulative sum of a list of matrices
##'
##' @param list_of_mats A list of matrices
##' @return A list where each element is the sum of matrices up to that index
##' @examples
##' mat_cumsum(list(matrix(1:4,2,2), matrix(5:8,2,2)))
mat_cumsum <- function(list_of_mats) {
  lapply(seq_along(list_of_mats), function(i) {
    Reduce("+", list_of_mats[1:i])
  })
}

##' Extract diagonal elements from a list of lists
##'
##' @param list_of_lists A list of lists (square)
##' @return A list of the diagonal elements
##' @examples
##' get_diag(list(list(1,2), list(3,4)))
get_diag <- function(list_of_lists) {
  out_list <- list()
  for (i in seq_along(list_of_lists[[1]])) {
    out_list[[i]] <- list_of_lists[[i]][[i]]
  }
  return(out_list)
}

##' Convert a matrix to a list of its columns
##'
##' @param X A matrix
##' @return A list where each element is a column of X
##' @examples
##' matrix_to_lov(matrix(1:4,2,2))
matrix_to_lov <- function(X) {
  lapply(seq_len(ncol(X)), function(i) X[,i])
}

##' Bind output from a list of lists by type
##'
##' @param list_of_lists A list of lists
##' @param type Name of element to extract and bind
##' @return A matrix with columns bound by type
##' @examples
##' # bind_output(list(list(a=1:2), list(a=3:4)), 'a')
bind_output <- function(list_of_lists, type) {
  do.call("cbind", lapply(list_of_lists, function(x) do.call("c", x[[type]])))
}

#' Internals
#'
#' @return List of internal functions to be exported
#' @export
.internals <- function() {
  internals <- list()
  internals[["get_list_of_eta"]] <- get_list_of_eta
  internals[["eval_deriv"]] <- eval_deriv
  internals[["eta_to_alpha"]] <- eta_to_alpha
  internals[["list_times_list"]] <- list_times_list
  internals[["matrix_to_lov"]] <- matrix_to_lov
  internals[["bind_output"]] <- bind_output
  return(internals)
}

#' Title
#'
#' @param list_of_lists
#'
#' @return
#' @export
#'
#' @examples
bind_list <- function(list_of_lists) {
  do.call("rbind", lapply(list_of_lists, function(x) do.call("cbind", x)))
}


#' Plot a set of weights as a stacked area.
#'
#' @param W An N by K matrix containing outputted weights.
#'
#' @return A ggplot stacked area plot.
#' @export
plot_weights <- function(W, date_vec = NULL) {
  # Requires: dplyr, tidyr, ggplot2
  if (!requireNamespace("dplyr", quietly = TRUE) ||
      !requireNamespace("tidyr", quietly = TRUE) ||
      !requireNamespace("ggplot2", quietly = TRUE)) {
    stop("plot_weights requires dplyr, tidyr, and ggplot2 packages.")
  }
  W_df <- as.data.frame(W)
  W_df <- dplyr::rename_with(W_df, ~ gsub("^V", "Expert ", .), dplyr::starts_with("V"))
  if (!is.null(date_vec)) {
    W_df$Date <- date_vec
  } else {
    W_df$Date <- seq_len(nrow(W_df))
  }
  W_df <- tidyr::pivot_longer(W_df, cols = -"Date", names_to = "Model", values_to = "Weight")
  out_plot <- ggplot2::ggplot(W_df, ggplot2::aes(x = .data$Date, y = .data$Weight, fill = .data$Model)) +
    ggplot2::geom_area() + ggplot2::theme_classic()
  return(out_plot)
}
