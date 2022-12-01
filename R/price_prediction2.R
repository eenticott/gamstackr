library(ggplot2)
library(mgcv)
library(lubridate)
library(dplyr)

price_data <- read.csv(file = "Data/central_gas_price_1005_R.csv")

price_data_H0 <- price_data %>% filter(Hour == 0)

price_data_H0 <- price_data_H0 %>% mutate(Date = as_datetime(ymd(Date)))

price_data_H0$LoadMatrix <- price_data_H0 %>% dplyr::select(starts_with("Prev_Residual_Load_H")) %>% as.matrix
price_data_H0$LagMatrix <- matrix(1:24, nrow(price_data_H0), 24, byrow = TRUE)

price_data_H0$SpotMatrixL2 <- price_data_H0 %>% dplyr::select(starts_with("Lag_J2_SpotPrice_H")) %>% as.matrix
price_data_H0$SpotMatrixL7 <- price_data_H0 %>% dplyr::select(starts_with("Lag_J7_SpotPrice_H")) %>% as.matrix

fit_func1 <- function(data) {
  data <- tail(data, 14) # take last two weeks
  gam(SpotPrice ~ s(LagMatrix, by = LoadMatrix), data = data)
}

fit_func2 <- function(data) {
  data <- tail(data, 12*7)
  gam(SpotPrice ~ ti(LagMatrix, LoadMatrix, mc = c(FALSE, TRUE)) +
        ti(LagMatrix, SpotMatrixL2, mc = c(FALSE, TRUE)) +
        ti(LagMatrix, SpotMatrixL7, mc = c(FALSE, TRUE)) +
        s(Lag_J1_GazPrice, k = 3) + s(Nuclear_availability, k = 3) + s(toy_sin), data = data,
      family = gaussian())
}

fit_func3 <- function(data) {
  data <- tail(data, 52*7)
  gam(SpotPrice ~ ti(LagMatrix, LoadMatrix, mc = c(FALSE, TRUE)) +
        ti(LagMatrix, SpotMatrixL2, mc = c(FALSE, TRUE)) +
        ti(LagMatrix, SpotMatrixL7, mc = c(FALSE, TRUE)) +
        s(Lag_J1_GazPrice, k = 3) + s(Nuclear_availability, k = 3) + s(toy_sin), data = data,
      family = gaussian())
}

dens_func <- function(fitted_model, data) {
  mu <- predict(fitted_model, data)
  sd <- sqrt(fitted_model$sigma2)
  dnorm(data[,"SpotPrice"], mu, sd, log = TRUE)
}


ex1 <- create_expert(fit_func = fit_func1, dens_func = dens_func)
ex2 <- create_expert(fit_func = fit_func2, dens_func = dens_func)
ex3 <- create_expert(fit_func = fit_func3, dens_func = dens_func)


bind_df <- function(list_of_lists, type) {
  do.call("cbind", lapply(list_of_lists, function(x) do.call("c", x[[type]])))
}
bind_df(out, "dens")

windower <- create_windower(52*7, horizon_size = 7, window_size = 52*7, step_size = 7, type = "sliding")

out <- evaluate_expert(list(ex1, ex2, ex3), windower, price_data_H0, type = c("density", "predict"))
out
stacker <- create_stack(list(ex1, ex2, ex3), inners = list(ordinal(3)))

stack_out <- evaluate_stack(stacker, data = cbind(stack, past_performance), windower = windower)

## Output of evaluate_stack
# most_recent_stack_model
# weights out of sample
# density as matrix
#
