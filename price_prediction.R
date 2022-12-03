library(ggplot2)
library(mgcv)
library(lubridate)
library(dplyr)
library(mgcViz)
library(zoo)

# cha dayto scat dist
# add day of week to experts
# fit to each hour
# fit MVN stacking

my_data <- read.csv(file = "Data/central_gas_price_1005_R.csv")

dat_save <- my_data %>% filter(Hour == 0)

dat_h <- my_data %>% filter(Hour == 0)
dat_h <- dat_h %>% mutate(Date = as_datetime(ymd(Date)))
ggplot(data = my_data, aes(Date, SpotPrice)) + geom_line()
plot1 <- ggplot(data = dat_h, aes(Date, SpotPrice)) + geom_line()
print(plot1)
plot(dat_h$SpotPrice)

dat_h$LoadMatrix <- dat_h %>% dplyr::select(starts_with("Prev_Residual_Load_H")) %>% as.matrix
dat_h$LagMatrix <- matrix(1:24, nrow(dat_h), 24, byrow = TRUE)

dat_h$SpotMatrixL2 <- dat_h %>% dplyr::select(starts_with("Lag_J2_SpotPrice_H")) %>% as.matrix
dat_h$SpotMatrixL7 <- dat_h %>% dplyr::select(starts_with("Lag_J7_SpotPrice_H")) %>% as.matrix

start_date <- ymd("2019-01-01")

my_gam1 <- function(date_idx) {
  gam(SpotPrice ~ s(LagMatrix, by = LoadMatrix), data = dat_h[date_idx, ])
}

my_gam <- function(date_idx) {
  gam(SpotPrice ~ ti(LagMatrix, LoadMatrix, mc = c(FALSE, TRUE)) +
        ti(LagMatrix, SpotMatrixL2, mc = c(FALSE, TRUE)) +
        ti(LagMatrix, SpotMatrixL7, mc = c(FALSE, TRUE)) +
        s(Lag_J1_GazPrice, k = 3) + s(Nuclear_availability, k = 3) + s(toy_sin), data = dat_h[date_idx, ],
      family = gaussian())
}

K <- 3
dens_matrix <- matrix(nrow = nrow(dat_h), ncol = K)
predict_matrix <- matrix(nrow = nrow(dat_h), ncol = K)
i <- 1
cur_date <- start_date
while (cur_date < ymd("2022-01-01")) {
  print(i)
  #w1_lag_idx <- dat_h$Date %within% interval(cur_date - days(3), cur_date - days(1))
  #w1_lag_mod <- my_gam0(w1_lag_idx)

  w2_lag_idx <- dat_h$Date %within% interval(cur_date - weeks(2), cur_date - days(1))
  w2_lag_mod <- my_gam1(w2_lag_idx)

  w4_lag_idx <- dat_h$Date %within% interval(cur_date - weeks(12), cur_date - days(1))
  w4_lag_mod <- my_gam1(w4_lag_idx)

  #w10_lag_idx <- dat_h$Date %within% interval(cur_date - weeks(26), cur_date - days(1))
  #w10_lag_mod <- my_gam(w10_lag_idx)

  w20_lag_idx <- dat_h$Date %within% interval(cur_date - weeks(52), cur_date - days(1))
  w20_lag_mod <- my_gam(w20_lag_idx)

  # Predict week ahead
  new_week_idx <- dat_h$Date %within% interval(cur_date, cur_date + days(6))
  #p1 <- predict(w1_lag_mod, newdata = dat_h[new_week_idx, ])
  p2 <- predict(w2_lag_mod, newdata = dat_h[new_week_idx, ])
  p3 <- predict(w4_lag_mod, newdata = dat_h[new_week_idx, ])
  #p4 <- predict(w10_lag_mod, newdata = dat_h[new_week_idx, ])
  p5 <- predict(w20_lag_mod, newdata = dat_h[new_week_idx, ])

  #den1 <- dnorm(dat_h[new_week_idx, ]$SpotPrice, p1, summary(w1_lag_mod)$sigma, log = FALSE)
  den2 <- dnorm(dat_h[new_week_idx, ]$SpotPrice, p2, sqrt(w2_lag_mod$sig2), log = FALSE)
  den3 <- dnorm(dat_h[new_week_idx, ]$SpotPrice, p3, sqrt(w4_lag_mod$sig2), log = FALSE)
  #den4 <- dnorm(dat_h[new_week_idx, ]$SpotPrice, p4, sqrt(w10_lag_mod$sig2), log = FALSE)
  den5 <- dnorm(dat_h[new_week_idx, ]$SpotPrice, p5, sqrt(w20_lag_mod$sig2), log = FALSE)

  predict_matrix[new_week_idx, 1:K] <- cbind(p2, p3, p5)
  dens_matrix[new_week_idx, 1:K] <- cbind(den2, den3, den5)

  cur_date <- cur_date + weeks(1)
  i <- i + 1
}
dat_h[-1, "pp1"] <- rollmean((predict_matrix - dat_h$SpotPrice)**2, k = 1)[-nrow(rollmean((predict_matrix - dat_h$SpotPrice)**2, k = 1)), 1]
dat_h[-1, "pp2"] <- rollmean((predict_matrix - dat_h$SpotPrice)**2, k = 1)[-nrow(rollmean((predict_matrix - dat_h$SpotPrice)**2, k = 1)), 2]
dat_h[-(1:3), "pp3"] <- rollmean((predict_matrix - dat_h$SpotPrice)**2, k = 3)[-nrow(rollmean((predict_matrix - dat_h$SpotPrice)**2, k = 3)), 1]
dat_h[-(1:3), "pp4"] <- rollmean((predict_matrix - dat_h$SpotPrice)**2, k = 3)[-nrow(rollmean((predict_matrix - dat_h$SpotPrice)**2, k = 3)), 2]
# dat_h[-1, "pp3"] <- ((predict_matrix - dat_h$SpotPrice)**2)[-nrow(dens_matrix), 3]

cur_date <- start_date + months(12) + days(1)

weight_matrix <- matrix(nrow = nrow(dat_h), ncol = K)


for (i in 1:104) {
  print(i)
  inners <- list(ordinal(K))
  list_of_densities <- list(as.matrix(na.omit(as.data.frame(dens_matrix[(dat_h$Date < cur_date), ])))[-c(1:3), ])
  pre_fam <- NestedStack(list_of_densities, inners, RidgePen = 1e-05)
  fit_stack <- gam(list(SpotPrice ~ s(pp1, pp2) + s(pp3, pp4) + s(toy_sin)),
                   data = na.omit(dat_h[dat_h$Date < cur_date, c("SpotPrice", "pp1", "pp2","pp3","pp4", "toy_sin")]),
                   family = pre_fam,
                   control = gam.control(trace = FALSE))

  new_week_idx <- dat_h$Date %within% interval(cur_date, cur_date + days(6))

  weight_matrix[new_week_idx, ] <- ordinal(K)(
    predict(fit_stack, newdata = dat_h[new_week_idx, ]),
    fit_stack$coefficients[(length(fit_stack$coefficients)):length(fit_stack$coefficients)],
    deriv = 0)$f_eval

  cur_date <- cur_date + weeks(1)
}

idx <- dat_h$Date %within% interval("2019-01-01", "2022-01-01")

m_preds <- predict_matrix[idx,] * sd(dat_save$SpotPrice) + mean(dat_save$SpotPrice)
preds <- (rowSums(predict_matrix * weight_matrix)[idx] * sd(dat_save$SpotPrice)) + mean(dat_save$SpotPrice)
preds2 <- rowSums(m_preds * weight_matrix[idx,])
date <- dat_h$Date[idx]
y <- (dat_h$SpotPrice[idx] * sd(dat_save$SpotPrice)) + mean(dat_save$SpotPrice)

test_dat <- data.frame(y_true = y, y_pred = preds, date = date)

ggplot(data = test_dat) + geom_line(aes(date, y)) + geom_line(aes(date, y_pred), col = "red") + ylab("Spot Price") + xlab("Date")

ggplot(data = test_dat) + geom_line(aes(date, y_pred - y_true), col = "red")

sqrt(mean((m_preds[!is.na(preds2),3] - y[!is.na(preds2)])**2))
sqrt(mean((preds2[!is.na(preds2)] - y[!is.na(preds2)])**2, na.rm = F))

weight_dat <- data.frame(weight_matrix)
weight_dat$date <- dat_h$Date
wide_weight <- melt(setDT(weight_dat), id.vars = "date")
levels(wide_weight$variable) <- c("2 weeks", "3 months", "1 year")
wide_weight <- rename(wide_weight, "Expert" = variable, "Weight" = value)

wplot <- ggplot(data = wide_weight) + geom_area(aes(x = date, y = Weight, fill = Expert))
wplot <- wplot + guides(shape = guide_legend(override.aes = list(size = 0.3)))
wplot <- wplot + guides(color = guide_legend(override.aes = list(size = 0.3)))
wplot <- wplot + theme(legend.title = element_text(size = 4),
                       legend.text = element_text(size = 4)) + xlab("Date")

ggsave("weight_plot.png", wplot)

pplot <- ggplot(data = test_dat) + geom_line(aes(date, y)) + geom_line(aes(date, y_pred), col = "red", alpha = 0.7) + ylab("Spot Price") + xlab("Date")
ggsave("pred_plot.png", pplot)
ggplot(data = weight_dat[idx, ]) +
  geom_line(aes(date, X1), col = "#1eccaf") +
  geom_line(aes(date, X2), col = "blue") +
  geom_line(aes(date, X3), col = "green") +
  geom_line(aes(date, X4), col = "purple")


