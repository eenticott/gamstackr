    ## Loading required package: nlme

    ## This is mgcv 1.8-41. For overview type 'help("mgcv-package")'.

    ## Loading required package: timechange

    ## 
    ## Attaching package: 'lubridate'

    ## The following objects are masked from 'package:base':
    ## 
    ##     date, intersect, setdiff, union

    ## 
    ## Attaching package: 'dplyr'

    ## The following object is masked from 'package:nlme':
    ## 
    ##     collapse

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

    ## 
    ## Attaching package: 'zoo'

    ## The following objects are masked from 'package:base':
    ## 
    ##     as.Date, as.Date.numeric

    ## 
    ## Attaching package: 'gamstackr'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     id

    ## 
    ## Attaching package: 'xgboost'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     slice

We will demonstrate how to use the `gamstackr` package through a short
example applied to day-ahead electricity spot price forecasting. We look
at hourly data from 2018 - 2022. Price data can be found at
`https://www.entsoe.eu/`, exogenous data was provided by EDF, some of
which is available at part of the eco2mix data set found here:
`https://odre.opendatasoft.com/explore/dataset/eco2mix-regional-tr`/.

![](forecasting_notebook_files/figure-markdown_strict/unnamed-chunk-1-1.png)

## *gamstackr* workflow

The basic workflow for fitting models includes four key steps:

1.  **Define experts**

2.  **Fit experts**

3.  **Fit the stack**

4.  **Forecast expert densities and weights for desired period**

We will go through each of these steps in more detail later.

## Time series workflow

Often when dealing with time series data you will want to perform these
steps iteratively.

The time series workflow would be

1.  **Define experts**

2.  **Define windows**

3.  *For each window*:

    1.  **Fit experts**
    2.  **Fit stack**
    3.  **Forecast expert densities and weights on chosen horizon**

There are tools included in `gamstackr` for handling windows and
automatically fitting/predicting across a given window.

# Window tools

There are two main types of window that we consider, *sliding window*
and *expanding window*. You first need to create your window function
with given parameters, and then you can call this on any data frame that
you wish to split into windows.

    slide_window <- create_windower(12*7, horizon_size = 7, window_size = 12*7, step_size = 7, type = "sliding")
    df_slide <- slide_window(price_data_H)
    plot_windows(df_slide)

![](forecasting_notebook_files/figure-markdown_strict/unnamed-chunk-2-1.png)

Sliding windows keep the length of data the same as you move through
time. This is useful when older data becomes less relevant as you move
through time.

    expand_window <- create_windower(12*7, horizon_size = 7, window_size = 12*7, step_size = 7, type = "expanding")
    df_expand <- expand_window(price_data_H)
    plot_windows(df_expand)

![](forecasting_notebook_files/figure-markdown_strict/unnamed-chunk-3-1.png)

Expanding windows keep adding more data to your training set, this is
useful when you don’t expect the trend to change.

# Creating experts

The next step is to define the experts that you want to use in your
stacking, we include a flexible framework that allows for almost any
model to be compatable at the cost of some extra work when defining
models.

An expert object is created through the function `create_expert` which
requires the user to provide two functions, a fitting function
`fit_func` and a density function `dens_func`.

The `fit_func` should take as input a data frame containing both the
response and the features.

The `dens_func` should take as input a fitted model, and the data for
which you want the densities, it should return a vector of the density
at each row of the data given.

You can optionally specify a custom `pred_func` if your predictions
aren’t directly obtained by calling `predict(fitted_model)`.

Once you have defined these two functions you can create your expert.

    fit_func <- function(data) {
      data <- tail(data, 7*4) # trim to last 4 weeks of data
      gam(SpotPrice ~  s(LagMatrix, by = LoadMatrix) + Lag_J1_GazPrice + Nuclear_availability + M1_Oil + M1_Coal, data = data)
    }


    dens_func <- function(fitted_model, data) {
      mu <- predict(fitted_model, data)
      sd <- sqrt(fitted_model$sig2)
      dnorm(data[,"SpotPrice"], mu, sd, log = TRUE)
    }


    expert1 <- create_expert(fit_func, dens_func = dens_func)

An `expert` is an object with the following methods:

-   `fit`
-   `predict`
-   `simulate`
-   `density`

For example to fit it to the first window of data we can call,

    expert1 <- expert1$fit(expert1, df_slide[attr(df_slide, "training_windows")[[1]],])

The package is better set up for fitting and predicting an expert across
an entire window at once, we can do this with `evaluate_expert`. This
takes an expert object, window object and a data frame and then
sequentially fits to each window. By default this will only return the
test predictions from each window, however we can also ask it to return
each fitted model object and/or the evaluated density.

    exp1_out <- evaluate_expert(expert1, slide_window, price_data_H, type = c("model", "predict", "density"))

    ##   |                                                          |                                                  |   0%  |                                                          |=                                                 |   2%  |                                                          |==                                                |   5%  |                                                          |====                                              |   8%  |                                                          |=====                                             |  10%  |                                                          |======                                            |  12%  |                                                          |========                                          |  15%  |                                                          |=========                                         |  18%  |                                                          |==========                                        |  20%  |                                                          |===========                                       |  22%  |                                                          |============                                      |  25%  |                                                          |==============                                    |  28%  |                                                          |===============                                   |  30%  |                                                          |================                                  |  32%  |                                                          |==================                                |  35%  |                                                          |===================                               |  38%  |                                                          |====================                              |  40%  |                                                          |=====================                             |  42%  |                                                          |======================                            |  45%  |                                                          |========================                          |  48%  |                                                          |=========================                         |  50%  |                                                          |==========================                        |  52%  |                                                          |============================                      |  55%  |                                                          |=============================                     |  58%  |                                                          |==============================                    |  60%  |                                                          |===============================                   |  62%  |                                                          |================================                  |  65%  |                                                          |==================================                |  68%  |                                                          |===================================               |  70%  |                                                          |====================================              |  72%  |                                                          |======================================            |  75%  |                                                          |=======================================           |  78%  |                                                          |========================================          |  80%  |                                                          |=========================================         |  82%  |                                                          |==========================================        |  85%  |                                                          |============================================      |  88%  |                                                          |=============================================     |  90%  |                                                          |==============================================    |  92%  |                                                          |================================================  |  95%  |                                                          |================================================= |  98%  |                                                          |==================================================| 100%

Then we can retrieve relevant information from our list of outputs.

    exp1_out$preds[1:3] # predictions for first 3 test windows

    ## [[1]]
    ##       85       86       87       88       89       90       91 
    ## 27.29038 35.93240 34.26679 38.55249 25.73162 21.59209 20.81119 
    ## 
    ## [[2]]
    ##       92       93       94       95       96       97       98 
    ## 32.58839 31.75952 30.03517 32.10625 28.71492 20.24273 14.80221 
    ## 
    ## [[3]]
    ##       99      100      101      102      103      104      105 
    ## 41.25832 34.98615 26.87001 31.35441 31.53950 17.39344 10.69683

    exp1_out$dens[1:3] # evaluated density on first 3 test windows

    ## [[1]]
    ## [1] -3.141503 -2.747355 -2.808424 -3.477677 -2.977127 -2.770492 -4.455639
    ## 
    ## [[2]]
    ## [1] -3.540320 -5.147748 -4.848701 -6.889294 -4.230990 -3.613639 -2.808886
    ## 
    ## [[3]]
    ## [1] -3.065358 -2.810018 -3.720443 -2.957988 -2.877830 -2.965067 -3.657297

    exp1_out$model # final fitted model

    ## 
    ## Family: gaussian 
    ## Link function: identity 
    ## 
    ## Formula:
    ## SpotPrice ~ s(LagMatrix, by = LoadMatrix) + Lag_J1_GazPrice + 
    ##     Nuclear_availability + M1_Oil + M1_Coal
    ## 
    ## Estimated degrees of freedom:
    ## 3.28  total = 7.28 
    ## 
    ## GCV score: 17.84823     rank: 14/15

Now we need to create a stack, before we can create a stack we will need
some more experts, we will create an expert with the same structure but
increase the window size. When adding multiple experts, it is
recommended that they all take the same data frame as input, and to
transform inside `fit_func` if required.

    fit_func2 <- function(data) {
      data <- tail(data, 8*4) # trim to last 8 weeks of data
      gam(SpotPrice ~  s(LagMatrix, by = LoadMatrix) + Lag_J1_GazPrice + Nuclear_availability + M1_Oil + M1_Coal, data = data)
    }

    fit_func3 <- function(data) {
      data <- tail(data, 12*4) # trim to last 12 weeks of data
      gam(SpotPrice ~  s(LagMatrix, by = LoadMatrix) + Lag_J1_GazPrice + Nuclear_availability + M1_Oil + M1_Coal, data = data)
    }

    expert2 <- create_expert(fit_func2, dens_func = dens_func)
    expert3 <- create_expert(fit_func3, dens_func = dens_func)

We can fit multiple experts at once with `evaluate_expert`, simply enter
them as a list.

    experts_out <- evaluate_expert(list(expert1, expert2, expert3), slide_window, price_data_H, type = c("predict", "density"))

    ##   |                                                          |                                                  |   0%  |                                                          |=                                                 |   2%  |                                                          |==                                                |   5%  |                                                          |====                                              |   8%  |                                                          |=====                                             |  10%  |                                                          |======                                            |  12%  |                                                          |========                                          |  15%  |                                                          |=========                                         |  18%  |                                                          |==========                                        |  20%  |                                                          |===========                                       |  22%  |                                                          |============                                      |  25%  |                                                          |==============                                    |  28%  |                                                          |===============                                   |  30%  |                                                          |================                                  |  32%  |                                                          |==================                                |  35%  |                                                          |===================                               |  38%  |                                                          |====================                              |  40%  |                                                          |=====================                             |  42%  |                                                          |======================                            |  45%  |                                                          |========================                          |  48%  |                                                          |=========================                         |  50%  |                                                          |==========================                        |  52%  |                                                          |============================                      |  55%  |                                                          |=============================                     |  58%  |                                                          |==============================                    |  60%  |                                                          |===============================                   |  62%  |                                                          |================================                  |  65%  |                                                          |==================================                |  68%  |                                                          |===================================               |  70%  |                                                          |====================================              |  72%  |                                                          |======================================            |  75%  |                                                          |=======================================           |  78%  |                                                          |========================================          |  80%  |                                                          |=========================================         |  82%  |                                                          |==========================================        |  85%  |                                                          |============================================      |  88%  |                                                          |=============================================     |  90%  |                                                          |==============================================    |  92%  |                                                          |================================================  |  95%  |                                                          |================================================= |  98%  |                                                          |==================================================| 100%
    ## ----- Expert 1 completed ----- 
    ##   |                                                          |                                                  |   0%  |                                                          |=                                                 |   2%  |                                                          |==                                                |   5%  |                                                          |====                                              |   8%  |                                                          |=====                                             |  10%  |                                                          |======                                            |  12%  |                                                          |========                                          |  15%  |                                                          |=========                                         |  18%  |                                                          |==========                                        |  20%  |                                                          |===========                                       |  22%  |                                                          |============                                      |  25%  |                                                          |==============                                    |  28%  |                                                          |===============                                   |  30%  |                                                          |================                                  |  32%  |                                                          |==================                                |  35%  |                                                          |===================                               |  38%  |                                                          |====================                              |  40%  |                                                          |=====================                             |  42%  |                                                          |======================                            |  45%  |                                                          |========================                          |  48%  |                                                          |=========================                         |  50%  |                                                          |==========================                        |  52%  |                                                          |============================                      |  55%  |                                                          |=============================                     |  58%  |                                                          |==============================                    |  60%  |                                                          |===============================                   |  62%  |                                                          |================================                  |  65%  |                                                          |==================================                |  68%  |                                                          |===================================               |  70%  |                                                          |====================================              |  72%  |                                                          |======================================            |  75%  |                                                          |=======================================           |  78%  |                                                          |========================================          |  80%  |                                                          |=========================================         |  82%  |                                                          |==========================================        |  85%  |                                                          |============================================      |  88%  |                                                          |=============================================     |  90%  |                                                          |==============================================    |  92%  |                                                          |================================================  |  95%  |                                                          |================================================= |  98%  |                                                          |==================================================| 100%
    ## ----- Expert 2 completed ----- 
    ##   |                                                          |                                                  |   0%  |                                                          |=                                                 |   2%  |                                                          |==                                                |   5%  |                                                          |====                                              |   8%  |                                                          |=====                                             |  10%  |                                                          |======                                            |  12%  |                                                          |========                                          |  15%  |                                                          |=========                                         |  18%  |                                                          |==========                                        |  20%  |                                                          |===========                                       |  22%  |                                                          |============                                      |  25%  |                                                          |==============                                    |  28%  |                                                          |===============                                   |  30%  |                                                          |================                                  |  32%  |                                                          |==================                                |  35%  |                                                          |===================                               |  38%  |                                                          |====================                              |  40%  |                                                          |=====================                             |  42%  |                                                          |======================                            |  45%  |                                                          |========================                          |  48%  |                                                          |=========================                         |  50%  |                                                          |==========================                        |  52%  |                                                          |============================                      |  55%  |                                                          |=============================                     |  58%  |                                                          |==============================                    |  60%  |                                                          |===============================                   |  62%  |                                                          |================================                  |  65%  |                                                          |==================================                |  68%  |                                                          |===================================               |  70%  |                                                          |====================================              |  72%  |                                                          |======================================            |  75%  |                                                          |=======================================           |  78%  |                                                          |========================================          |  80%  |                                                          |=========================================         |  82%  |                                                          |==========================================        |  85%  |                                                          |============================================      |  88%  |                                                          |=============================================     |  90%  |                                                          |==============================================    |  92%  |                                                          |================================================  |  95%  |                                                          |================================================= |  98%  |                                                          |==================================================| 100%
    ## ----- Expert 3 completed -----

`experts_out` is a list of lists, where the top level contains each of
the evaluated experts and the inner level contains the
predictions/density etc. We can turn these into data frames using the
`bind_output` function, this is the format we need to input to the
stacking.

    preds <- bind_output(experts_out, "preds")
    dens <- bind_output(experts_out, "dens")

    head(preds)

    ##        [,1]     [,2]     [,3]
    ## 85 27.29038 27.71385 28.86378
    ## 86 35.93240 34.32104 35.65636
    ## 87 34.26679 32.31312 32.36643
    ## 88 38.55249 36.69834 37.91584
    ## 89 25.73162 26.23389 30.17281
    ## 90 21.59209 20.73251 26.67370

    head(dens)

    ##           [,1]      [,2]      [,3]
    ## [1,] -3.141503 -3.213065 -3.382221
    ## [2,] -2.747355 -2.805436 -2.849941
    ## [3,] -2.808424 -2.800873 -2.847972
    ## [4,] -3.477677 -3.816736 -3.547068
    ## [5,] -2.977127 -2.964273 -2.846361
    ## [6,] -2.770492 -2.805496 -3.317935

# Stacking

Similiarly to the expert we need to create a stacking object.

This takes a list of lists experts (designed for nested stacking) and a
list of *inner\_functions* these are the functions used to obtain a set
of inner weights, currently we have 3 methods available: - `ordinal(n)`:
Assumes an ordered set of experts, `n` is number of experts included -
`MVN_weights(x)`: Requires experts to be placed in some coordinate
space, `x` is given coordinates - `id()` which is used when you want to
include models in the stack with fixed weights.

    stacker <- create_stacker(list(list(expert1, expert2, expert3)), inners = list(ordinal(3)))

The stacker object can be used to fit experts to new data or to fit a
weighted stack. We can fit it to the experts we defined earlier. We want
to fit to out of sample density to avoid overfitting.

    stack_dat <- price_data_H[rownames(preds),]

We need to give the stacker a formula, a list of out of sample densities
and

    stack <- evaluate_stack(stacker,
                   formula = list(SpotPrice ~ VMA),
                   windower = expand_window, 
                   stack_data = stack_dat, 
                   list_of_densities = list(dens))

    ##   |                                                          |                                                  |   0%  |                                                          |==                                                |   4%  |                                                          |====                                              |   7%  |                                                          |=====                                             |  11%  |                                                          |=======                                           |  14%  |                                                          |=========                                         |  18%  |                                                          |===========                                       |  21%  |                                                          |============                                      |  25%  |                                                          |==============                                    |  29%  |                                                          |================                                  |  32%  |                                                          |==================                                |  36%  |                                                          |====================                              |  39%  |                                                          |=====================                             |  43%  |                                                          |=======================                           |  46%  |                                                          |=========================                         |  50%  |                                                          |===========================                       |  54%  |                                                          |=============================                     |  57%  |                                                          |==============================                    |  61%  |                                                          |================================                  |  64%  |                                                          |==================================                |  68%  |                                                          |====================================              |  71%  |                                                          |======================================            |  75%  |                                                          |=======================================           |  79%  |                                                          |=========================================         |  82%  |                                                          |===========================================       |  86%  |                                                          |=============================================     |  89%  |                                                          |==============================================    |  93%  |                                                          |================================================  |  96%  |                                                          |==================================================| 100%

Our fitted stack contains all the evaluated weights.

    stack_weights <- do.call("rbind", stack)
    stack_idx <- rownames(stack_weights)
    plot(price_data_H[stack_idx, "VMA"], stack_weights[,1])

![](forecasting_notebook_files/figure-markdown_strict/unnamed-chunk-15-1.png)

    metrics( price_data_H[stack_idx, "SpotPrice"], rowSums(stack_weights * preds[stack_idx,]))

    ## [1] "RMSE: 8.28239727693582"
    ## [1] "MAE: 6.27734475072377"

    metrics( price_data_H[stack_idx, "SpotPrice"], preds[stack_idx,3])

    ## [1] "RMSE: 8.47467798427476"
    ## [1] "MAE: 6.34155051940999"

Only a small improvement but this was a very simplified example. We can
include other models alongside this ordered stack.

    fit_func4 <- function(data) {
      data <- data %>% select(-c("Time", "hour_pred", "split_type", "Date", "clock", "Trend"))
      X <- data %>% select(-c("SpotPrice"))
      y <- data %>% select("SpotPrice")
      N <- nrow(X)
      train_idx <- 1:floor(0.9*N)
      test_idx <- -train_idx
      dTrain <- xgb.DMatrix(data = as.matrix(X[train_idx,]), label = as.matrix(y[train_idx,]))
      dTest <- xgb.DMatrix(data = as.matrix(X[test_idx,]), label = as.matrix(y[test_idx,]))
      watchlist = list(train = dTrain, test = dTest)
      mod <- xgb.train(data = dTrain, verbose = FALSE, eta = 0.05, max_depth = 5, objective = "reg:squarederror", nrounds = 200, watchlist = watchlist, early_stopping_rounds = 10,
      )
      mod$fitted_data <- dTrain
      mod$sig2 <- sd(predict(mod, dTest) - y[test_idx,])**2
      return(mod)
    }

    dens_func4 <- function(fitted_model, data) {
      mu <- pred_func4(fitted_model, data)
      sd <- sqrt(fitted_model$sig2)
      dnorm(data[,"SpotPrice"], mu, sd, log = TRUE)
    }

    pred_func4 <- function(fitted_model, data) {
      data <- data %>% select(-c("Time", "hour_pred", "split_type", "Date", "clock", "Trend"))
      X <- data %>% select(-c("SpotPrice"))
      y <- data %>% select("SpotPrice")
      dTest <- xgb.DMatrix(data = as.matrix(X), label = as.matrix(y))
      predict(fitted_model, dTest)
    }

As we are using xgboost we need a custom predict function.

    expert4 <- create_expert(fit_func = fit_func4, dens_func = dens_func4, pred_func = pred_func4)

Evaluate all the experts.

    out <- evaluate_expert(list(expert1, expert2, expert3, expert4), expand_window, price_data_H, type = c("density", "predict", "model"))

    ##   |                                                          |                                                  |   0%  |                                                          |=                                                 |   2%  |                                                          |==                                                |   5%  |                                                          |====                                              |   8%  |                                                          |=====                                             |  10%  |                                                          |======                                            |  12%  |                                                          |========                                          |  15%  |                                                          |=========                                         |  18%  |                                                          |==========                                        |  20%  |                                                          |===========                                       |  22%  |                                                          |============                                      |  25%  |                                                          |==============                                    |  28%  |                                                          |===============                                   |  30%  |                                                          |================                                  |  32%  |                                                          |==================                                |  35%  |                                                          |===================                               |  38%  |                                                          |====================                              |  40%  |                                                          |=====================                             |  42%  |                                                          |======================                            |  45%  |                                                          |========================                          |  48%  |                                                          |=========================                         |  50%  |                                                          |==========================                        |  52%  |                                                          |============================                      |  55%  |                                                          |=============================                     |  58%  |                                                          |==============================                    |  60%  |                                                          |===============================                   |  62%  |                                                          |================================                  |  65%  |                                                          |==================================                |  68%  |                                                          |===================================               |  70%  |                                                          |====================================              |  72%  |                                                          |======================================            |  75%  |                                                          |=======================================           |  78%  |                                                          |========================================          |  80%  |                                                          |=========================================         |  82%  |                                                          |==========================================        |  85%  |                                                          |============================================      |  88%  |                                                          |=============================================     |  90%  |                                                          |==============================================    |  92%  |                                                          |================================================  |  95%  |                                                          |================================================= |  98%  |                                                          |==================================================| 100%
    ## ----- Expert 1 completed ----- 
    ##   |                                                          |                                                  |   0%  |                                                          |=                                                 |   2%  |                                                          |==                                                |   5%  |                                                          |====                                              |   8%  |                                                          |=====                                             |  10%  |                                                          |======                                            |  12%  |                                                          |========                                          |  15%  |                                                          |=========                                         |  18%  |                                                          |==========                                        |  20%  |                                                          |===========                                       |  22%  |                                                          |============                                      |  25%  |                                                          |==============                                    |  28%  |                                                          |===============                                   |  30%  |                                                          |================                                  |  32%  |                                                          |==================                                |  35%  |                                                          |===================                               |  38%  |                                                          |====================                              |  40%  |                                                          |=====================                             |  42%  |                                                          |======================                            |  45%  |                                                          |========================                          |  48%  |                                                          |=========================                         |  50%  |                                                          |==========================                        |  52%  |                                                          |============================                      |  55%  |                                                          |=============================                     |  58%  |                                                          |==============================                    |  60%  |                                                          |===============================                   |  62%  |                                                          |================================                  |  65%  |                                                          |==================================                |  68%  |                                                          |===================================               |  70%  |                                                          |====================================              |  72%  |                                                          |======================================            |  75%  |                                                          |=======================================           |  78%  |                                                          |========================================          |  80%  |                                                          |=========================================         |  82%  |                                                          |==========================================        |  85%  |                                                          |============================================      |  88%  |                                                          |=============================================     |  90%  |                                                          |==============================================    |  92%  |                                                          |================================================  |  95%  |                                                          |================================================= |  98%  |                                                          |==================================================| 100%
    ## ----- Expert 2 completed ----- 
    ##   |                                                          |                                                  |   0%  |                                                          |=                                                 |   2%  |                                                          |==                                                |   5%  |                                                          |====                                              |   8%  |                                                          |=====                                             |  10%  |                                                          |======                                            |  12%  |                                                          |========                                          |  15%  |                                                          |=========                                         |  18%  |                                                          |==========                                        |  20%  |                                                          |===========                                       |  22%  |                                                          |============                                      |  25%  |                                                          |==============                                    |  28%  |                                                          |===============                                   |  30%  |                                                          |================                                  |  32%  |                                                          |==================                                |  35%  |                                                          |===================                               |  38%  |                                                          |====================                              |  40%  |                                                          |=====================                             |  42%  |                                                          |======================                            |  45%  |                                                          |========================                          |  48%  |                                                          |=========================                         |  50%  |                                                          |==========================                        |  52%  |                                                          |============================                      |  55%  |                                                          |=============================                     |  58%  |                                                          |==============================                    |  60%  |                                                          |===============================                   |  62%  |                                                          |================================                  |  65%  |                                                          |==================================                |  68%  |                                                          |===================================               |  70%  |                                                          |====================================              |  72%  |                                                          |======================================            |  75%  |                                                          |=======================================           |  78%  |                                                          |========================================          |  80%  |                                                          |=========================================         |  82%  |                                                          |==========================================        |  85%  |                                                          |============================================      |  88%  |                                                          |=============================================     |  90%  |                                                          |==============================================    |  92%  |                                                          |================================================  |  95%  |                                                          |================================================= |  98%  |                                                          |==================================================| 100%
    ## ----- Expert 3 completed ----- 
    ##   |                                                          |                                                  |   0%  |                                                          |=                                                 |   2%  |                                                          |==                                                |   5%  |                                                          |====                                              |   8%  |                                                          |=====                                             |  10%  |                                                          |======                                            |  12%  |                                                          |========                                          |  15%  |                                                          |=========                                         |  18%  |                                                          |==========                                        |  20%  |                                                          |===========                                       |  22%  |                                                          |============                                      |  25%  |                                                          |==============                                    |  28%  |                                                          |===============                                   |  30%  |                                                          |================                                  |  32%  |                                                          |==================                                |  35%  |                                                          |===================                               |  38%  |                                                          |====================                              |  40%  |                                                          |=====================                             |  42%  |                                                          |======================                            |  45%  |                                                          |========================                          |  48%  |                                                          |=========================                         |  50%  |                                                          |==========================                        |  52%  |                                                          |============================                      |  55%  |                                                          |=============================                     |  58%  |                                                          |==============================                    |  60%  |                                                          |===============================                   |  62%  |                                                          |================================                  |  65%  |                                                          |==================================                |  68%  |                                                          |===================================               |  70%  |                                                          |====================================              |  72%  |                                                          |======================================            |  75%  |                                                          |=======================================           |  78%  |                                                          |========================================          |  80%  |                                                          |=========================================         |  82%  |                                                          |==========================================        |  85%  |                                                          |============================================      |  88%  |                                                          |=============================================     |  90%  |                                                          |==============================================    |  92%  |                                                          |================================================  |  95%  |                                                          |================================================= |  98%  |                                                          |==================================================| 100%
    ## ----- Expert 4 completed -----

Split the densities into data frames corresponding to groups of experts,
in this case experts 1 - 3 belong to the ordered group and expert 4
stands alone.

    dens <- bind_output(out, "dens")
    preds <- bind_output(out, "preds")

    list_of_densities <- list(dens[,1:3], dens[,4,drop = FALSE])

    stack_dat <- price_data_H[rownames(preds),]

The list of lists for the experts now makes more sense, you also need a
list of the inner functions as before. `id` simply assigns a weight of
1.

    stacker <- create_stacker(list(list(expert1, expert2, expert3), list(expert4)), inners = list(ordinal(3), id()))

    st_out <- evaluate_stack(stacker, list(SpotPrice ~ 1, ~ VMA), windower = expand_window, stack_data = stack_dat, list_of_densities = list_of_densities)

    ##   |                                                          |                                                  |   0%  |                                                          |==                                                |   4%  |                                                          |====                                              |   7%  |                                                          |=====                                             |  11%  |                                                          |=======                                           |  14%  |                                                          |=========                                         |  18%  |                                                          |===========                                       |  21%  |                                                          |============                                      |  25%  |                                                          |==============                                    |  29%  |                                                          |================                                  |  32%  |                                                          |==================                                |  36%  |                                                          |====================                              |  39%  |                                                          |=====================                             |  43%  |                                                          |=======================                           |  46%  |                                                          |=========================                         |  50%  |                                                          |===========================                       |  54%  |                                                          |=============================                     |  57%  |                                                          |==============================                    |  61%  |                                                          |================================                  |  64%  |                                                          |==================================                |  68%  |                                                          |====================================              |  71%  |                                                          |======================================            |  75%  |                                                          |=======================================           |  79%  |                                                          |=========================================         |  82%  |                                                          |===========================================       |  86%  |                                                          |=============================================     |  89%  |                                                          |==============================================    |  93%  |                                                          |================================================  |  96%  |                                                          |==================================================| 100%

    stack_weights <- do.call("rbind", st_out)

    head(stack_weights)

    ##           [,1]         [,2]      [,3]      [,4]
    ## 169 0.03385924 1.678391e-07 0.8062688 0.1598718
    ## 170 0.03343485 1.658226e-07 0.8066932 0.1598718
    ## 171 0.03648869 1.802833e-07 0.8036393 0.1598718
    ## 172 0.03553032 1.757575e-07 0.8045977 0.1598718
    ## 173 0.03507731 1.736143e-07 0.8050507 0.1598718
    ## 174 0.03737048 1.844374e-07 0.8027575 0.1598718

    metrics( price_data_H[stack_idx, "SpotPrice"], rowSums(stack_weights * preds[stack_idx,]))

    ## [1] "RMSE: 7.78850369238737"
    ## [1] "MAE: 5.83607844905244"

    metrics( price_data_H[stack_idx, "SpotPrice"], preds[stack_idx,4])

    ## [1] "RMSE: 10.5808976672181"
    ## [1] "MAE: 7.88642634138769"

    with(price_data_H[stack_idx, ], plot(Date, SpotPrice, type = "l"))
    lines(price_data_H[stack_idx, "Date"], rowSums(stack_weights * preds[stack_idx,]), col = "red")

![](forecasting_notebook_files/figure-markdown_strict/unnamed-chunk-23-1.png)

    W <- do.call("rbind", st_out)
    W_dat <- data.frame(Date = price_data_H[rownames(W), "Date"])
    W_dat <- cbind(W_dat, W)
    W_dat <- reshape::melt(W_dat, id.vars = c("Date"))
    ggplot(data=W_dat, aes(x=Date, y = value, fill=variable, group = variable)) +
      geom_area() + scale_fill_brewer(palette="Set2")

![](forecasting_notebook_files/figure-markdown_strict/unnamed-chunk-24-1.png)
