set.seed(43)

test_that("evaluate_stack runs end-to-end with provided preds", {
  # Synthetic data
  N <- 240
  x <- rnorm(N)
  y <- 0.7 * x + rnorm(N, sd = 0.6)
  df <- data.frame(y = y, x = x)

  # Two simple experts with different fits
  fit1 <- function(data) lm(y ~ x, data = data)
  fit2 <- function(data) lm(y ~ poly(x, 2, raw = TRUE), data = data)

  pred_from <- function(fitted_model, data) stats::predict(fitted_model, newdata = data)

  ex1 <- create_expert(fit_func = fit1, pred_func = pred_from)
  ex2 <- create_expert(fit_func = fit2, pred_func = pred_from)

  # Build predictions matrix for stacking using a rolling window
  windower <- create_windower(initial_size = 120, horizon_size = 24, step_size = 24, type = "expanding")
  windowed_df <- windower(df)
  n_windows <- length(attr(windowed_df, "training_windows"))

  # For stacking we need stack_data and a preds matrix aligned to its rows
  # We'll form preds for the whole dataset to simplify checks, then subset per window inside evaluate_stack
  # Fit both experts on the full data for this synthetic test (we pass preds explicitly, so expert fitting is not required)
  full_pred1 <- pred_from(fit1(df), df)
  full_pred2 <- pred_from(fit2(df), df)
  preds <- cbind(full_pred1, full_pred2)
  colnames(preds) <- c("ex1", "ex2")
  expect_equal(nrow(preds), nrow(df))

  # Stacker with multinomial weights over 2 experts, type = loss using squared loss
  stacker <- create_stacker(experts = list(ex1, ex2), weight_func = multinomial(2), type = "loss", loss = gamstackr:::square_loss)

  # Evaluate stack predicting weights per window; formula can be as simple as ~ x
  out <- evaluate_stack(stacker, formula = y ~ x, windower = windower, stack_data = df, expert_outputs = preds, RidgePen = 1e-5)

  # Structure checks
  expect_type(out, "list")
  expect_length(out, n_windows)

  # Each element is an N_test x K matrix with K=2
  test_sizes <- vapply(attr(windowed_df, "testing_windows"), length, integer(1))
  elem_dims <- lapply(out, dim)
  expect_true(all(vapply(elem_dims, function(d) d[2] == 2, logical(1))))
  expect_equal(vapply(elem_dims, function(d) d[1], integer(1)), test_sizes)

  # Weights should be finite and between 0 and 1 for multinomial
  all_W <- do.call(rbind, out)
  expect_true(all(is.finite(all_W)))
  expect_true(all(all_W >= 0 - 1e-8))
  expect_true(all(all_W <= 1 + 1e-8))
})


