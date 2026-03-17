library(readxl)
library(forecast)
library(tseries)
library(stats)
library(ggplot2)
library(Metrics)
library(patchwork)

# Consider only the price column in the data
data <- read_excel("arima-data.xlsx")
View(data)
series <- ts(data$PRICE, start = c(2015,1), frequency = 12)
series

plot(series)

par(mfrow = c(2, 2))
acf(series, lag.max = 48, ci.type = "white", main = "White Noise CI") 
acf(series, lag.max = 48,  ci.type = "ma", main = "Bartlett's CI")
pacf(series, lag.max = 48)

series_diff = diff(series, 1, 1)
plot(series_diff)
plot(series_diff)
acf(series_diff, lag.max = 48, ci.type = "white", main = "White Noise CI") 
acf(series_diff, lag.max = 48,  ci.type = "ma", main = "Bartlett's CI")
pacf(series_diff, lag.max = 48)
par(mfrow = c(1, 1))


train <- window(series, end = c(2022, 12))
test  <- window(series, start = c(2023, 1))

arima_fit <- auto.arima(train,
                        d=1,
                        ic='bic',
                        lambda=NULL,
                        stepwise = FALSE,
                        trace = TRUE)


arima_forecast <- forecast(arima_fit, h = length(test))

actual_fit <- as.numeric(train)
predicted <- as.numeric(arima_fit$fitted)
actual_test <- as.numeric(test)
forecasted <- as.numeric(arima_forecast$mean)

results <- data.frame(
  Set  = c("Training", "Test"),
  MAE  = c(mae(actual_fit, predicted),   mae(actual_test, forecasted)),
  RMSE = c(rmse(actual_fit, predicted),  rmse(actual_test, forecasted)),
  MAPE = c(mape(actual_fit, predicted),  mape(actual_test, forecasted)) * 100
)

print(results)

plot(forecast(arima_fit, h = 12))

checkresiduals(arima_fit)


# Use this section only for manual implementation
arima.fit <- Arima(train,
                   order = c(2, 1, 0),
                   method = 'CSS-ML')

arima_forecast <- forecast(arima.fit, h = length(test))

actual_fit <- as.numeric(train)
predicted <- as.numeric(arima.fit$fitted)
actual_test <- as.numeric(test)
forecasted <- as.numeric(arima_forecast$mean)

results <- data.frame(
  Set  = c("Training", "Test"),
  MAE  = c(mae(actual_fit, predicted),   mae(actual_test, forecasted)),
  RMSE = c(rmse(actual_fit, predicted),  rmse(actual_test, forecasted)),
  MAPE = c(mape(actual_fit, predicted),  mape(actual_test, forecasted)) * 100
)

plot(forecast(arima.fit, h = 12))

checkresiduals(arima.fit)
