library(ggplot2)
library(mgcv)
library(lubridate)
library(dplyr)
library(zoo)
library(gamstackr)

# TODO: Add 20% cv into xgboost fit
# TODO: Create slides

metrics <- function(y_true, y_pred) {
  RMSE <- sqrt(mean((y_true - y_pred)**2))
  print(paste("RMSE:", RMSE))
}

# Define experts ===============================================================
fit_func1 <- function(data) {
  data <- tail(data, 14) # take last two weeks
  gam(SpotPrice ~ s(LagMatrix, by = LoadMatrix), data = data)
}

fit_func2 <- function(data) {
  data <- tail(data, 12*7)
  gam(SpotPrice ~ ti(LagMatrix, LoadMatrix, mc = c(FALSE, TRUE)) + s(LagMatrix, by = SpotMatrixL2) + s(LagMatrix, by = SpotMatrixL7) +
        s(Lag_J1_GazPrice, k = 3) + s(Nuclear_availability, k = 3) + s(toy_sin), data = data,
      family = gaussian())
}

fit_func3 <- function(data) {
  data <- tail(data, 52*7)
  gam(SpotPrice ~ ti(LagMatrix, LoadMatrix, mc = c(FALSE, TRUE)) +
        ti(LagMatrix, SpotMatrixL2, mc = c(FALSE, TRUE)) +
        ti(LagMatrix, SpotMatrixL7, mc = c(FALSE, TRUE)) +
        s(Lag_J1_GazPrice, k = 10) + s(Nuclear_availability, k = 3) + s(toy_sin, k = 10), data = data,
      family = gaussian())
}

fit_func4 <- function(data) {
  data <- data %>% select(-c("Time", "hour_pred", "split_type", "Date"))
  X <- data %>% select(-c("SpotPrice"))
  y <- data %>% select("SpotPrice")
  dTrain <- xgb.DMatrix(data = as.matrix(X), label = as.matrix(y))
  mod <- xgb.train(data = dTrain, max.depth = 5, eta = 0.1, objective = "reg:squarederror", nrounds = 100)
  return(mod)
}

pred_func4 <- function(fitted_model, data) {
  data <- data %>% select(-c("Time", "hour_pred", "split_type", "Date"))
  X <- data %>% select(-c("SpotPrice"))
  y <- data %>% select("SpotPrice")
  dTest <- xgb.DMatrix(data = as.matrix(X), label = as.matrix(y))
  predict(fitted_model, dTest)
}

dens_func <- function(fitted_model, data) {
  mu <- predict(fitted_model, data)
  sd <- sqrt(fitted_model$sig2)
  dnorm(data[,"SpotPrice"], mu, sd, log = TRUE)
}

ex1 <- create_expert(fit_func = fit_func1, dens_func = dens_func)
ex2 <- create_expert(fit_func = fit_func2, dens_func = dens_func)
ex3 <- create_expert(fit_func = fit_func3, dens_func = dens_func)

ex4 <- create_expert(fit_func = fit_func4, dens_func = dens_func, pred_func = pred_func4)
# ==============================================================================

price_data <- read.csv(file = "Data/central_gas_price_1005_R.csv")

price_data <- price_data %>% mutate(Date = as_datetime(ymd(Date)))
price_data <- price_data[order(price_data$Date),]

price_data$LoadMatrix <- price_data %>% dplyr::select(starts_with("Prev_Residual_Load_H")) %>% as.matrix
price_data$LagMatrix <- matrix(1:24, nrow(price_data), 24, byrow = TRUE)

price_data$SpotMatrixL2 <- price_data %>% dplyr::select(starts_with("Lag_J2_SpotPrice_H")) %>% as.matrix
price_data$SpotMatrixL7 <- price_data %>% dplyr::select(starts_with("Lag_J7_SpotPrice_H")) %>% as.matrix
price_data$Trend <- 1:nrow(price_data)

price_data_H <- price_data %>% filter(Hour == 0)

# Fit experts
ex_windower <- create_windower(52*7, horizon_size = 7, window_size = 52*7, step_size = 7, type = "sliding")
out <- evaluate_expert(list(ex1, ex2, ex3), ex_windower, price_data_H, type = c("density", "predict"))

ex_windower <- create_windower(52*7*24, horizon_size = 24*7, window_size = 52*7*24, step_size = 24*7, type = "sliding")
out4 <- evaluate_expert(ex4, ex_windower, price_data, type = "predict")

list_of_densities <- list(bind_output(out, "dens"))
dens <- bind_output(out, "dens")
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

plot_dat <- cbind(price_data_H[rownames(predictions),c("Date", "SpotPrice")], predictions)
colnames(plot_dat) <- c("Date", "SpotPrice", "Ex1", "Ex2", "Ex3")
ggplot(plot_dat) + geom_line(aes(Date, SpotPrice))
ggplot(plot_dat) + geom_line(aes(Date, SpotPrice)) + geom_line(aes(Date, Ex1), col = "red") + geom_line(aes(Date, Ex2), col = "blue") + geom_line(aes(Date, Ex3), col = "green")
metrics(plot_dat$SpotPrice, plot_dat$Ex1)
metrics(plot_dat$SpotPrice, plot_dat$Ex2)
metrics(plot_dat$SpotPrice, plot_dat$Ex3)


stack_windower <- create_windower(20*7, horizon_size = 7, window_size = NULL, step_size = 7, type = "expanding")
stacker <- create_stacker(list(list(ex1, ex2, ex3)), inners = list(ordinal(3)))
sout <- evaluate_stack(stacker, list(SpotPrice ~ s(PP2.1) + s(PP3.1) + s(PP2.7) + s(PP3.7)), stack_windower, stack_data, list_of_densities)

W <- do.call("rbind", sout)

metrics(price_data_H[rownames(W), "SpotPrice"], rowSums(W * predictions[rownames(W),]))
## Output of evaluate_stack
# most_recent_stack_model
# weights out of sample
# density as matrix


