# list utils

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

list_by_vector <- function(l, v) {
  out <- list()
  for (i in 1:length(v)) {
    out[[i]] <- l[[i]] * v[i]
  }
  return(out)
}


list_by_list <- function(list1, list2) {
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
    out[[i]] <- t(list1[[i]]) %*% list2[[i]]
  }

  return(out)
}

list_times_list <- function(list1, list2) {
  out <- list()
  for (i in 1:length(list1)) {
    out[[i]] <- list1[[i]] * list2[[i]]
  }
  out
}

repmat <- function(mat, v = 1, h = 1) {
  kronecker(matrix(1, h, v), mat)
}

# Matrix utils
# Found in https://stackoverflow.com/questions/24299171/function-to-split-a-matrix-into-sub-matrices-in-r
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

mat_cumsum <- function(list_of_mats) {
  lapply(as.list(1:length(list_of_mats)), function(i) {
    Reduce("+", list_of_mats[1:i])
  })
}

get_diag <- function(list_of_lists) {
  out_list <- list()
  for (i in 1:length(list_of_lists[[1]])) {
    out_list[[i]] <- list_of_lists[[i]][[i]]
  }
  return(out_list)
}

matrix_to_lov <- function(X) {
  lapply(seq_len(ncol(X)), function(i) X[,i])
}