library(ggplot2)
library(mgcv)
library(lubridate)
library(dplyr)
library(zoo)
library(gamstackr)
library(xgboost)

metrics <- function(y_true, y_pred) {
  RMSE <- sqrt(mean((y_true - y_pred)**2))
  print(paste("RMSE:", RMSE))
  MAE <- mean(abs(y_true - y_pred))
  print(paste("MAE:", MAE))
}

# Define experts ===============================================================
fit_func1 <- function(data) {
  data <- tail(data, 14) # take last two weeks
  gam(SpotPrice ~ s(LagMatrix, by = LoadMatrix), data = data)
}

fit_func2 <- function(data) {
  data <- tail(data, 12*7)
  gam(SpotPrice ~ ti(LagMatrix, LoadMatrix, mc = c(FALSE, TRUE)) + s(LagMatrix, by = SpotMatrixL2) + s(LagMatrix, by = SpotMatrixL7) +
        s(Lag_J1_GazPrice, k = 3) + s(Nuclear_availability, k = 3) + s(yday), data = data,
      family = gaussian())
}

fit_func3 <- function(data) {
  data <- tail(data, 52*7)
  gam(SpotPrice ~ ti(LagMatrix, LoadMatrix, mc = c(FALSE, TRUE)) +
        ti(LagMatrix, SpotMatrixL2, mc = c(FALSE, TRUE)) +
        ti(LagMatrix, SpotMatrixL7, mc = c(FALSE, TRUE)) +
        s(Lag_J1_GazPrice, k = 10) + s(Nuclear_availability, k = 3) + s(yday, k = 10), data = data,
      family = gaussian())
}

fit_func4 <- function(data) {
  data <- data %>% select(-c("Time", "hour_pred", "split_type", "Date"))
  X <- data %>% select(-c("SpotPrice"))
  y <- data %>% select("SpotPrice")
  N <- nrow(X)
  train_idx <- 1:floor(0.95*N)
  test_idx <- -train_idx
  dTrain <- xgb.DMatrix(data = as.matrix(X[train_idx,]), label = as.matrix(y[train_idx,]))
  dTest <- xgb.DMatrix(data = as.matrix(X[test_idx,]), label = as.matrix(y[test_idx,]))
  watchlist = list(train = dTrain, test = dTest)
  mod <- xgb.train(data = dTrain, verbose = FALSE, max.depth = 5, eta = 0.1, objective = "reg:squarederror", nrounds = 200, watchlist = watchlist, early_stopping_rounds = 10)
  mod$fitted_data <- dTrain
  mod$sig2 <- sd(predict(mod, dTest) - y[test_idx,])**2
  return(mod)
}

fit_func5 <- function(data) {
  data <- tail(data, 7*4)
  gam(SpotPrice ~ s(LagMatrix, by = LoadMatrix) + Lag_J1_GazPrice + Nuclear_availability + M1_Oil + M1_Coal, data = data)
}

fit_func6 <- function(data) {
  data <- tail(data, 7*12)
  gam(SpotPrice ~ s(LagMatrix, by = LoadMatrix) + Lag_J1_GazPrice + Nuclear_availability + M1_Oil + M1_Coal, data = data)
}

fit_func7 <- function(data) {
  gam(SpotPrice ~ s(LagMatrix, by = LoadMatrix) + Lag_J1_GazPrice + Nuclear_availability + M1_Oil + M1_Coal + s(yday, k = 5), data = data)
}

# Dens func -----------------------------------

dens_func <- function(fitted_model, data) {
  mu <- predict(fitted_model, data)
  sd <- sqrt(fitted_model$sig2)
  dnorm(data[,"SpotPrice"], mu, sd, log = TRUE)
}

dens_func4 <- function(fitted_model, data) {
  mu <- pred_func4(fitted_model, data)
  sd <- sqrt(fitted_model$sig2)
  dnorm(data[,"SpotPrice"], mu, sd, log = TRUE)
}

# Preds func --------------------------------

pred_func4 <- function(fitted_model, data) {
  data <- data %>% select(-c("Time", "hour_pred", "split_type", "Date"))
  X <- data %>% select(-c("SpotPrice"))
  y <- data %>% select("SpotPrice")
  dTest <- xgb.DMatrix(data = as.matrix(X), label = as.matrix(y))
  predict(fitted_model, dTest)
}

# Create experts ---------------------------------------

ex1 <- create_expert(fit_func = fit_func1, dens_func = dens_func)
ex2 <- create_expert(fit_func = fit_func2, dens_func = dens_func)
ex3 <- create_expert(fit_func = fit_func3, dens_func = dens_func)

ex4 <- create_expert(fit_func = fit_func4, dens_func = dens_func4, pred_func = pred_func4)

ex5 <- create_expert(fit_func = fit_func5, dens_func = dens_func)
ex6 <- create_expert(fit_func = fit_func6, dens_func = dens_func)
ex7 <- create_expert(fit_func = fit_func7, dens_func = dens_func)

# ==============================================================================

price_data <- read.csv(file = "Data/central_gas_price_1005_R.csv")

price_data <- price_data %>% mutate(Date = as_datetime(ymd(Date)), yday = yday(Date))
price_data <- price_data[order(price_data$Date),]

price_data$LoadMatrix <- price_data %>% dplyr::select(starts_with("Prev_Residual_Load_H")) %>% as.matrix
price_data$LagMatrix <- matrix(1:24, nrow(price_data), 24, byrow = TRUE)

price_data$SpotMatrixL2 <- price_data %>% dplyr::select(starts_with("Lag_J2_SpotPrice_H")) %>% as.matrix
price_data$SpotMatrixL7 <- price_data %>% dplyr::select(starts_with("Lag_J7_SpotPrice_H")) %>% as.matrix
price_data$Trend <- 1:nrow(price_data)

price_data_H <- price_data %>% filter(Hour == 0)

# ==============================================================================

# Fit experts
ex_windower <- create_windower(52*7, horizon_size = 7, window_size = 52*7, step_size = 7, type = "sliding")

out <- evaluate_expert(list(ex1, ex2, ex3, ex4, ex5, ex6, ex7), ex_windower, price_data_H, type = c("density", "predict"))
preds <- bind_output(out, "preds")

for (i in 1:3) {
  metrics(price_data_H[rownames(preds), "SpotPrice"], preds[,i])
}
dens <- bind_output(out, "dens")
list_of_densities <- list(dens[,1:3], dens[,5:7], dens[,4, drop = FALSE])
predictions <- bind_output(out, "preds")

add_past_performance <- function(stack_data, dens,  over = c(1, 3, 7)) {
  out <- list()
  for (i in 1:length(over)) {
    pp <- lag(rollmean(dens, over[i], align = "right", fill = NA))
    out[[i]] <- (exp(pp)/rowSums(exp(pp)))[,-1]
  }
  out <- do.call("cbind", out)
  colnames(out) <- paste0("PP", apply(expand.grid(2:ncol(dens), over), 1, paste, collapse="."))
  out <- cbind(stack_data, out)
  return(out)
}

stack_data <- price_data_H[rownames(predictions),]
stack_data <- add_past_performance(stack_data, dens)

# get rid of NA
stack_data <- tail(stack_data, -7)
list_of_densities <- lapply(list_of_densities, function(x) tail(x, -7))


stack_windower <- create_windower(52*7, horizon_size = 7, window_size = NULL, step_size = 7, type = "expanding")
stacker <- create_stacker(list(list(ex1, ex2, ex3), list(ex5, ex6, ex7), list(ex4)), inners = list(ordinal(3), ordinal(3), id()))
sout <- evaluate_stack(stacker, list(SpotPrice ~ s(PP5.1) + s(PP5.7),~s(PP4.1) + s(PP4.1), ~ s(PP2.1) + s(PP3.1) + s(PP2.7) + s(PP3.7) + Lag_J2_SpotPrice_H0, ~ s(PP5.1) + s(PP5.7) ), stack_windower, stack_data, list_of_densities = list_of_densities)
sout <- evaluate_stack(stacker, list(SpotPrice ~ 1,~1, ~ 1, ~ 1), stack_windower, stack_data, list_of_densities = list_of_densities)

W <- do.call("rbind", sout)
W_dat <- data.frame(Date = price_data_H[rownames(W), "Date"])
W_dat <- cbind(W_dat, W)
W_dat <- reshape::melt(W_dat, id.vars = c("Date"))
ggplot(data=W_dat, aes(x=Date, y = value, fill=variable, group = variable)) +
  geom_area() + scale_fill_brewer(palette="Set2")

metrics(price_data_H[rownames(W), "SpotPrice"], rowSums(W * predictions[rownames(W),]))
metrics(price_data_H[rownames(W), "SpotPrice"], rowMeans(predictions[rownames(W),]))

ggplot() + geom_line(aes(price_data_H[rownames(W), "Date"],price_data_H[rownames(W), "SpotPrice"])) +
  geom_line(aes(price_data_H[rownames(W), "Date"], rowSums(W * predictions[rownames(W),])), col = "red")

## Output of evaluate_stack
# most_recent_stack_model
# weights out of sample
# density as matrix


