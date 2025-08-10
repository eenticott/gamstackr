set.seed(42)

test_that("evaluate_expert runs end-to-end on simple expert", {
  # Synthetic data
  N <- 200
  x <- rnorm(N)
  y <- 2 * x + rnorm(N, sd = 0.5)
  df <- data.frame(y = y, x = x)

  # Simple expert: lm with Gaussian density
  fit_func <- function(data) {
    lm(y ~ x, data = data)
  }
  dens_func <- function(fitted_model, data) {
    mu <- predict(fitted_model, newdata = data)
    # unbiased residual variance estimate
    sigma2 <- as.numeric(crossprod(fitted_model$residuals) / fitted_model$df.residual)
    stats::dnorm(data$y, mean = mu, sd = sqrt(sigma2), log = TRUE)
  }

  expert <- create_expert(fit_func = fit_func, dens_func = dens_func)

  # Windows
  windower <- create_windower(initial_size = 120, horizon_size = 20, step_size = 20, type = "expanding")
  windowed_df <- windower(df)
  n_windows <- length(attr(windowed_df, "training_windows"))

  out <- evaluate_expert(expert, windower, df, type = c("predict", "density", "model"))

  # Structure checks
  expect_type(out, "list")
  expect_true(all(c("preds", "dens", "model") %in% names(out)))
  expect_length(out$preds, n_windows)
  expect_length(out$dens, n_windows)

  # Per-window length checks match test-window sizes
  test_sizes <- vapply(attr(windowed_df, "testing_windows"), length, integer(1))
  pred_sizes <- vapply(out$preds, length, integer(1))
  dens_sizes <- vapply(out$dens, length, integer(1))
  expect_equal(pred_sizes, test_sizes)
  expect_equal(dens_sizes, test_sizes)

  # Model present and is an lm
  expect_s3_class(out$model, "lm")
})



