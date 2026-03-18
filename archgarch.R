# --- Libraries ---
library(readxl)
library(tseries)
library(FinTS)
library(forecast)
library(rugarch)

# --- Data ---
data  <- read_excel(file.choose())
x     <- diff(data$price)
y     <- ts(x, frequency = 12, start = c(2015, 2))
plot(y, main = "Differenced Price Series")

# ============================================================
# STEP 1: ARIMA on differenced data
# ============================================================
model_arima  <- auto.arima(y, d = 0, ic = "bic", stepwise = FALSE, trace = TRUE)
summary(model_arima)

# Extract residuals
arima_resid <- resid(model_arima)

# ============================================================
# STEP 2: Residual Diagnostics
# ============================================================
par(mfrow = c(2, 2))
plot(arima_resid,       main = "Residuals")
acf(arima_resid,        main = "ACF of Residuals")
acf(arima_resid^2,      main = "ACF of Squared Residuals")
pacf(arima_resid^2,     main = "PACF of Squared Residuals")
par(mfrow = c(1, 1))

# ARCH LM Test
ArchTest(arima_resid, lags = 12, demean = TRUE)

# ============================================================
# STEP 3: Auto ARCH Model (GARCH with beta fixed to 0)
# ============================================================
arch_spec <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1, 0)), # q=1, p=0 → pure ARCH(1)
  mean.model     = list(armaOrder = c(1, 1), include.mean = TRUE),
  distribution.model = "norm"
)
arch_fit <- ugarchfit(spec = arch_spec, data = y)
show(arch_fit)

# ============================================================
# STEP 4: Auto GARCH Model
# ============================================================
garch_spec <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1, 1)), # GARCH(1,1)
  mean.model     = list(armaOrder = c(1, 1), include.mean = TRUE),
  distribution.model = "norm"
)
garch_fit <- ugarchfit(spec = garch_spec, data = y)
show(garch_fit)

# ============================================================
# STEP 5: Volatility Plot
# ============================================================
par(mfrow = c(2, 1))
plot(sigma(arch_fit),  main = "ARCH(1,0) - Conditional Volatility",  ylab = "Volatility")
plot(sigma(garch_fit), main = "GARCH(1,1) - Conditional Volatility", ylab = "Volatility")
par(mfrow = c(1, 1))

# ============================================================
# STEP 6: Residual Diagnostics on GARCH
# ============================================================
garch_resid <- residuals(garch_fit, standardize = TRUE)

par(mfrow = c(2, 2))
plot(garch_resid,        main = "Standardized Residuals")
acf(garch_resid,         main = "ACF of Std. Residuals")
acf(garch_resid^2,       main = "ACF of Squared Std. Residuals")
pacf(garch_resid^2,      main = "PACF of Squared Std. Residuals")
par(mfrow = c(1, 1))

# Remaining ARCH effect?
ArchTest(garch_resid, lags = 12, demean = TRUE)
Box.test(garch_resid, lag = 12, type = "Ljung-Box")