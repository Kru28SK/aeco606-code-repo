# ── Libraries ──────────────────────────────────────────────────────────────────
library(forecast)
library(Metrics)

# ── 1. Load & Create Time Series ───────────────────────────────────────────────
data     <- read.csv("pundanuru_market_price.CSV")
price.ts <- ts(data$price, start = c(2008, 1), frequency = 12)

train <- window(price.ts, end   = c(2022, 12))
test  <- window(price.ts, start = c(2023,  1))

# ── 2. Visualise ───────────────────────────────────────────────────────────────
plot(price.ts)

par(mfrow = c(1, 2))
acf (price.ts, lag.max = 36)
pacf(price.ts, lag.max = 36)
par(mfrow = c(1, 1))

# ── Helper: print RMSE & MAE ───────────────────────────────────────────────────
acc <- function(label, actual, predicted) {
  cat("\n──", label, "──\n")
  cat("  RMSE:", round(rmse(actual, predicted), 4), "\n")
  cat("  MAE :", round(mae (actual, predicted), 4), "\n")
}

# ── 3. SES ─────────────────────────────────────────────────────────────────────
ses.model <- ses(train, h = length(test))

cat("\n── SES Coefficients ──\n"); print(ses.model$model$par)

autoplot(train, series = "Actual") + autolayer(fitted(ses.model), series = "Fitted")
autoplot(ses.model) + autolayer(test, series = "Actual")
checkresiduals(ses.model)

acc("SES – Train", as.numeric(train), as.numeric(fitted(ses.model)))
acc("SES – Test",  as.numeric(test),  as.numeric(ses.model$mean))

# ── 4. Holt's Method ───────────────────────────────────────────────────────────

# 4a. Undamped
holt.model <- holt(train, h = length(test))
cat("\n── Holt (undamped) Coefficients ──\n"); print(round(holt.model$model$par, 4))

autoplot(train, series = "Actual") + autolayer(fitted(holt.model), series = "Fitted")
autoplot(holt.model) + autolayer(test, series = "Actual")
checkresiduals(holt.model)

acc("Holt Undamped – Train", as.numeric(train), as.numeric(fitted(holt.model)))
acc("Holt Undamped – Test",  as.numeric(test),  as.numeric(holt.model$mean))

# 4b. Damped
holt.d <- holt(train, damped = TRUE, h = length(test))
cat("\n── Holt (damped) Coefficients ──\n"); print(round(holt.d$model$par, 4))

autoplot(train, series = "Actual") + autolayer(fitted(holt.d), series = "Fitted")
autoplot(holt.d) + autolayer(test, series = "Actual")
checkresiduals(holt.d)

acc("Holt Damped – Train", as.numeric(train), as.numeric(fitted(holt.d)))
acc("Holt Damped – Test",  as.numeric(test),  as.numeric(holt.d$mean))

# ── 5. Future Forecast (full series, 24 months) ────────────────────────────────
autoplot(holt(price.ts, h = 24))