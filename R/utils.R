# list utils

#' Subtract elements of one list from another
#'
#' This function subtracts corresponding elements of list l2 from list l1.
#'
#' @param l1 First list
#' @param l2 Second list (to be subtracted from l1)
#'
#' @return A list containing the differences between corresponding elements
#'
list_take_list <- function(l1, l2) {
  out <- list()
  if (length(l1) != length(l2)) {
    stop("List length didn't match")
  }
  for (i in 1:length(l1)) {
    out[[i]] <- l1[[i]] - l2[[i]]
  }
  out
}

#' Multiply each element of a list by corresponding element of a vector
#'
#' This function multiplies each element of list l by the corresponding element of vector v.
#'
#' @param l A list of elements
#' @param v A vector of multipliers
#'
#' @return A list containing the products
#'
list_by_vector <- function(l, v) {
  out <- list()
  for (i in 1:length(v)) {
    out[[i]] <- l[[i]] * v[i]
  }
  return(out)
}

#' Multiply elements of two lists
#'
#' This function multiplies corresponding elements of two lists, either elementwise
#' or using matrix multiplication.
#'
#' @param list1 First list of elements
#' @param list2 Second list of elements
#' @param elementwise Logical; if TRUE, performs elementwise multiplication instead of matrix multiplication
#'
#' @return A list containing the products
#'
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
  for (i in 1:n1) {
    if (elementwise) {
      out[[i]] <- list1[[i]] * list2[[i]]
    } else {
      out[[i]] <- t(list1[[i]]) %*% list2[[i]]
    }
  }

  return(out)
}

#' Elementwise multiplication of two lists
#'
#' This function performs elementwise multiplication of corresponding elements in two lists.
#'
#' @param list1 First list of elements
#' @param list2 Second list of elements
#'
#' @return A list containing the elementwise products
#'
list_times_list <- function(list1, list2) {
  out <- list()
  for (i in 1:length(list1)) {
    out[[i]] <- list1[[i]] * list2[[i]]
  }
  out
}

#' Replicate a matrix
#'
#' This function replicates a matrix vertically and horizontally.
#'
#' @param mat Input matrix
#' @param v Number of vertical replications
#' @param h Number of horizontal replications
#'
#' @return A replicated matrix
#'
repmat <- function(mat, v = 1, h = 1) {
  kronecker(matrix(1, h, v), mat)
}

# Matrix utils
# Found in https://stackoverflow.com/questions/24299171/function-to-split-a-matrix-into-sub-matrices-in-r
#' Split a matrix into sub-matrices
#'
#' This function splits a matrix into sub-matrices of specified dimensions.
#'
#' @param M Input matrix
#' @param r Number of rows in each sub-matrix
#' @param c Number of columns in each sub-matrix
#'
#' @return An array of sub-matrices
#'
mat_split <- function(M, r, c) {
  nr <- ceiling(nrow(M) / r)
  nc <- ceiling(ncol(M) / c)
  newM <- matrix(NA, nr * r, nc * c)
  newM[1:nrow(M), 1:ncol(M)] <- M

  div_k <- kronecker(matrix(seq_len(nr * nc), nr, byrow = TRUE), matrix(1, r, c))
  matlist <- split(newM, div_k)
  N <- length(matlist)
  mats <- unlist(matlist)
  dim(mats) <- c(r, c, N)
  return(mats)
}

#' Cumulative sum of a list of matrices
#'
#' This function calculates the cumulative sum of a list of matrices.
#'
#' @param list_of_mats A list of matrices
#'
#' @return A list of cumulative sums of matrices
#'
mat_cumsum <- function(list_of_mats) {
  lapply(as.list(1:length(list_of_mats)), function(i) {
    Reduce("+", list_of_mats[1:i])
  })
}

#' Extract diagonal elements from a list of lists
#'
#' This function extracts the diagonal elements from a list of lists.
#'
#' @param list_of_lists A list of lists
#'
#' @return A list containing the diagonal elements
#'
get_diag <- function(list_of_lists) {
  out_list <- list()
  for (i in 1:length(list_of_lists[[1]])) {
    out_list[[i]] <- list_of_lists[[i]][[i]]
  }
  return(out_list)
}

#' Convert a matrix to a list of vectors
#'
#' This function converts a matrix to a list of column vectors.
#'
#' @param X Input matrix
#'
#' @return A list of column vectors
#'
matrix_to_lov <- function(X) {
  lapply(seq_len(ncol(X)), function(i) X[,i])
}

#' Bind outputs from a list of lists
#'
#' This function binds outputs of a specific type from a list of lists.
#'
#' @param list_of_lists A list of lists containing outputs
#' @param type The type of output to extract and bind
#'
#' @return A combined matrix of outputs
#'
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

#' Bind a list of lists into a single matrix
#'
#' This function binds a list of lists into a single matrix by row-binding
#' the column-bound elements of each inner list.
#'
#' @param list_of_lists A list of lists to be bound
#'
#' @return A combined matrix
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
plot_weights <- function(W) {
  W_df <- as.data.frame(W) %>%  rename_with(~ gsub("^V", "Expert ", .), starts_with("V")) %>% mutate(Date = stack_test$Date) %>%
    pivot_longer(cols = -Date,names_to = "Model", values_to = "Weight")

  out_plot <- ggplot(W_df, aes(x = Date, y = Weight, fill = Model)) + geom_area() + theme_classic()
  return(out_plot)
}
