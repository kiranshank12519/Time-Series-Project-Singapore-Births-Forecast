# Analysis.R - Model fitting and forecast evaluation for TFR and TLB
# generates forecasts, saves figures and results to disk.

library(tseries)

if (!dir.exists("figures")) dir.create("figures")


## Load and clean data 

raw <- read.csv("BirthsAndFertilityRates.csv",
                na.strings = "na",
                stringsAsFactors = FALSE,
                check.names = FALSE)
labels_raw <- trimws(raw[, 1])

# extract row by exact name, reverse year order, drop 2025
extract_row <- function(name) {
  matches <- labels_raw == name
  values <- as.numeric(unlist(raw[matches, -1]))
  rev(values)[1:65]
}

fertility <- data.frame(
  Year = 1960:2024,
  TFR  = extract_row("Total Fertility Rate (TFR) (Per Female)"),
  TLB  = extract_row("Total Live-Births (Number)")
)

TFR <- ts(fertility$TFR, start = 1960, frequency = 1)
TLB <- ts(fertility$TLB, start = 1960, frequency = 1)

# 1960-2012 train, 2013-2024 test
TFR_train <- window(TFR, end = 2012)
TLB_train <- window(TLB, end = 2012)
TFR_test  <- window(TFR, start = 2013)
TLB_test  <- window(TLB, start = 2013)


# Instead of repeating the same code over and over for each model, I used
# ChatGPT to help write a helper function that fits each model, checks the
# residuals, produces forecasts and calculates RMSE.

# fits on log scale, forecasts back-transformed with exp()
fit_and_evaluate <- function(train, test, order, seasonal = NULL, label) {
  
  log_train <- log(train)
  
  if (is.null(seasonal)) {
    fit <- tryCatch(arima(log_train, order = order),
                    error = function(e) NULL)
  } else {
    fit <- tryCatch(arima(log_train, order = order, seasonal = seasonal),
                    error = function(e) NULL)
  }
  
  if (is.null(fit)) {
    return(list(label = label, fit = NULL, viable = FALSE))
  }
  
  n_params <- sum(order[c(1, 3)]) + 
    ifelse(is.null(seasonal), 0, sum(seasonal$order[c(1, 3)]))
  n <- length(log_train)
  
  # white noise check
  lb <- Box.test(residuals(fit), lag = 24, type = "Ljung-Box",
                 fitdf = n_params)
  
  # small-sample correction
  k <- n_params + 1
  aicc <- AIC(fit) + (2 * k * (k + 1)) / (n - k - 1)
  
  fcast <- exp(predict(fit, n.ahead = length(test))$pred)
  err <- as.numeric(test) - as.numeric(fcast)
  
  list(label    = label,
       fit      = fit,
       viable   = lb$p.value > 0.05,
       aic      = AIC(fit),
       aicc     = aicc,
       lb_p     = lb$p.value,
       forecast = fcast,
       rmse     = sqrt(mean(err^2)))
}


# I also used ChatGPT to help write this function so the results from
# each candidate model could be put into the same comparison table.
# This made it easier to compare AIC, AICc, Ljung-Box p-values and RMSE.

print_results <- function(results, name) {
  cat("\n", strrep("=", 80), "\n", sep = "")
  cat("Results for ", name, "\n", sep = "")
  cat(strrep("=", 80), "\n\n", sep = "")
  
  tbl <- do.call(rbind, lapply(results, function(r) {
    if (is.null(r$fit)) {
      data.frame(Model = r$label, AIC = NA, AICc = NA, LjungBox_p = NA,
                 Viable = "Fit failed", RMSE = NA)
    } else {
      data.frame(
        Model      = r$label,
        AIC        = round(r$aic, 2),
        AICc       = round(r$aicc, 2),
        LjungBox_p = round(r$lb_p, 4),
        Viable     = ifelse(r$viable, "Yes", "No"),
        RMSE       = signif(r$rmse, 4)
      )
    }
  }))
  print(tbl, row.names = FALSE)
  
  viable <- results[sapply(results, function(r) isTRUE(r$viable))]
  if (length(viable) > 0) {
    best <- viable[[which.min(sapply(viable, function(r) r$rmse))]]
    cat("\nBest viable model: ", best$label, "\n", sep = "")
  } else {
    cat("\nNo viable models.\n")
  }
  
  return(tbl)
}


## TFR candidate models 

# three families: non-seasonal d=1, non-seasonal d=2 (Route 2),
# SARIMA with seasonal MA, SARIMA with seasonal AR
tfr_results <- list(
  
  # non-seasonal, d=1
  fit_and_evaluate(TFR_train, TFR_test, order = c(0,1,13),
                   label = "ARIMA(0,1,13)"),
  fit_and_evaluate(TFR_train, TFR_test, order = c(13,1,0),
                   label = "ARIMA(13,1,0)"),
  fit_and_evaluate(TFR_train, TFR_test, order = c(0,1,11),
                   label = "ARIMA(0,1,11)"),
  fit_and_evaluate(TFR_train, TFR_test, order = c(12,1,1),
                   label = "ARIMA(12,1,1)"),
  
  # non-seasonal, d=2
  fit_and_evaluate(TFR_train, TFR_test, order = c(13,2,2),
                   label = "ARIMA(13,2,2)"),
  fit_and_evaluate(TFR_train, TFR_test, order = c(12,2,3),
                   label = "ARIMA(12,2,3)"),
  
  # SARIMA with seasonal MA
  fit_and_evaluate(TFR_train, TFR_test,
                   order = c(1,1,0),
                   seasonal = list(order = c(0,1,1), period = 12),
                   label = "SARIMA(1,1,0)(0,1,1)[12]"),
  fit_and_evaluate(TFR_train, TFR_test,
                   order = c(0,1,1),
                   seasonal = list(order = c(0,1,1), period = 12),
                   label = "SARIMA(0,1,1)(0,1,1)[12]"),
  fit_and_evaluate(TFR_train, TFR_test,
                   order = c(1,1,1),
                   seasonal = list(order = c(0,1,1), period = 12),
                   label = "SARIMA(1,1,1)(0,1,1)[12]"),
  
  # SARIMA with seasonal AR
  fit_and_evaluate(TFR_train, TFR_test,
                   order = c(1,1,0),
                   seasonal = list(order = c(1,1,0), period = 12),
                   label = "SARIMA(1,1,0)(1,1,0)[12]"),
  fit_and_evaluate(TFR_train, TFR_test,
                   order = c(0,1,1),
                   seasonal = list(order = c(1,1,0), period = 12),
                   label = "SARIMA(0,1,1)(1,1,0)[12]"),
  fit_and_evaluate(TFR_train, TFR_test,
                   order = c(1,1,1),
                   seasonal = list(order = c(1,1,0), period = 12),
                   label = "SARIMA(1,1,1)(1,1,0)[12]"),
  
  # full seasonal
  fit_and_evaluate(TFR_train, TFR_test,
                   order = c(1,1,1),
                   seasonal = list(order = c(1,1,1), period = 12),
                   label = "SARIMA(1,1,1)(1,1,1)[12]")
)

tfr_table <- print_results(tfr_results, "TFR")

# For TFR, nine of the thirteen models were viable after the Ljung-Box check.
# The best model was ARIMA(13,2,2), so the non-seasonal family performed best.
# The best SARIMA model was close behind, with RMSE 0.09650 compared with
# 0.08623. The d = 2 models were competitive, with ARIMA(13,2,2) giving the
# lowest RMSE overall.


##TLB candidate models 

tlb_results <- list(
  
  # non-seasonal, d=1
  fit_and_evaluate(TLB_train, TLB_test, order = c(0,1,13),
                   label = "ARIMA(0,1,13)"),
  fit_and_evaluate(TLB_train, TLB_test, order = c(13,1,0),
                   label = "ARIMA(13,1,0)"),
  fit_and_evaluate(TLB_train, TLB_test, order = c(0,1,11),
                   label = "ARIMA(0,1,11)"),
  fit_and_evaluate(TLB_train, TLB_test, order = c(12,1,1),
                   label = "ARIMA(12,1,1)"),
  
  # non-seasonal, d=2
  fit_and_evaluate(TLB_train, TLB_test, order = c(13,2,2),
                   label = "ARIMA(13,2,2)"),
  fit_and_evaluate(TLB_train, TLB_test, order = c(12,2,3),
                   label = "ARIMA(12,2,3)"),
  
  # SARIMA with seasonal MA
  fit_and_evaluate(TLB_train, TLB_test,
                   order = c(1,1,0),
                   seasonal = list(order = c(0,1,1), period = 12),
                   label = "SARIMA(1,1,0)(0,1,1)[12]"),
  fit_and_evaluate(TLB_train, TLB_test,
                   order = c(0,1,1),
                   seasonal = list(order = c(0,1,1), period = 12),
                   label = "SARIMA(0,1,1)(0,1,1)[12]"),
  fit_and_evaluate(TLB_train, TLB_test,
                   order = c(1,1,1),
                   seasonal = list(order = c(0,1,1), period = 12),
                   label = "SARIMA(1,1,1)(0,1,1)[12]"),
  
  # SARIMA with seasonal AR
  fit_and_evaluate(TLB_train, TLB_test,
                   order = c(1,1,0),
                   seasonal = list(order = c(1,1,0), period = 12),
                   label = "SARIMA(1,1,0)(1,1,0)[12]"),
  fit_and_evaluate(TLB_train, TLB_test,
                   order = c(0,1,1),
                   seasonal = list(order = c(1,1,0), period = 12),
                   label = "SARIMA(0,1,1)(1,1,0)[12]"),
  fit_and_evaluate(TLB_train, TLB_test,
                   order = c(1,1,1),
                   seasonal = list(order = c(1,1,0), period = 12),
                   label = "SARIMA(1,1,1)(1,1,0)[12]"),
  
  # full seasonal
  fit_and_evaluate(TLB_train, TLB_test,
                   order = c(1,1,1),
                   seasonal = list(order = c(1,1,1), period = 12),
                   label = "SARIMA(1,1,1)(1,1,1)[12]")
)

tlb_table <- print_results(tlb_results, "TLB")

# For TLB, seven of the thirteen models were viable after the Ljung-Box check.
# The best model was SARIMA(1,1,1)(0,1,1)_{12}, so the seasonal family
# performed best. The best viable non-seasonal model was ARIMA(0,1,13), with
# RMSE 6,386 compared with 2,730 for the SARIMA model. This shows that the
# 12-year seasonal structure gave a much better forecast for TLB.

# Picking the best viable model per series 

best_viable <- function(results) {
  viable <- results[sapply(results, function(r) isTRUE(r$viable))]
  if (length(viable) == 0) return(NULL)
  viable[[which.min(sapply(viable, function(r) r$rmse))]]
}

best_tfr <- best_viable(tfr_results)
best_tlb <- best_viable(tlb_results)


##residual diagnostics for best models I didnt check both ACF and PACF because
#the PACF more or less gave the same results as the ACF. To make the figure nice
#only ACF and Q-Q plot is created.

plot_diagnostics <- function(m, filename) {
  res <- residuals(m$fit)
  
  # save to file
  png(filename, width = 1800, height = 900, res = 200)
  par(mfrow = c(1,2), mar = c(4,4,3,1))
  acf(res, lag.max = 24, main = paste0("Residual ACF: ", m$label))
  qqnorm(res, main = paste0("Q-Q Plot: ", m$label))
  qqline(res, col = "red")
  par(mfrow = c(1,1))
  dev.off()
  
  # also draw to plots pane
  par(mfrow = c(1,2), mar = c(4,4,3,1))
  acf(res, lag.max = 24, main = paste0("Residual ACF: ", m$label))
  qqnorm(res, main = paste0("Q-Q Plot: ", m$label))
  qqline(res, col = "red")
  par(mfrow = c(1,1))
}

if (!is.null(best_tfr)) plot_diagnostics(best_tfr, 
                                         "figures/fig5_tfr_residuals.png")

# Looking at the ACF plot for the TFR residuals, all spikes are within the
# significance bounds. This suggests there is no remaining significant
# autocorrelation, so the residuals behave like white noise. The Q-Q plot also
# follows the reference line closely, with only small deviations at the tails.
# Overall, the diagnostics suggest that ARIMA(13,2,2) is a suitable model for
# the TFR series.

if (!is.null(best_tlb)) plot_diagnostics(best_tlb,
                                         "figures/fig6_tlb_residuals.png")

# Looking at the ACF plot for the TLB residuals, most spikes are within the
# significance bounds, although there is one borderline negative spike around
# lag 4. This suggests that the SARIMA(1,1,1)(0,1,1)_{12} model has captured
# most of the autocorrelation in the TLB series, but there may still be some
# small residual structure left. The Q-Q plot mostly follows the reference line
# in the middle, but there are noticeable deviations in the tails. Overall, the
# diagnostics are acceptable, but not as strong as the TFR model.

#Forecasting plots with prediction intervals

plot_forecast <- function(train, test, fit, label, name, filename) {
  
  fc_log <- predict(fit, n.ahead = length(test))
  
  # 80% PI uses z=1.282, 95% uses z=1.96
  lo80_log <- fc_log$pred - 1.282 * fc_log$se
  hi80_log <- fc_log$pred + 1.282 * fc_log$se
  lo95_log <- fc_log$pred - 1.960 * fc_log$se
  hi95_log <- fc_log$pred + 1.960 * fc_log$se
  
  fc   <- as.numeric(exp(fc_log$pred))
  lo80 <- as.numeric(exp(lo80_log))
  hi80 <- as.numeric(exp(hi80_log))
  lo95 <- as.numeric(exp(lo95_log))
  hi95 <- as.numeric(exp(hi95_log))
  
  train_yr <- as.numeric(time(train))
  test_yr  <- as.numeric(time(test))
  tr <- as.numeric(train)
  te <- as.numeric(test)
  
  # prepend last training value so lines connect at 2012
  last_yr  <- tail(train_yr, 1)
  last_val <- tail(tr, 1)
  test_yr_ext <- c(last_yr, test_yr)
  fc_ext      <- c(last_val, fc)
  te_ext      <- c(last_val, te)
  lo80_ext    <- c(last_val, lo80)
  hi80_ext    <- c(last_val, hi80)
  lo95_ext    <- c(last_val, lo95)
  hi95_ext    <- c(last_val, hi95)
  
  y_range <- range(c(tr, te, lo95, hi95))
  
  draw <- function() {
    plot(range(c(train_yr, test_yr)), y_range,
         type = "n", xlab = "Year", ylab = name,
         main = paste0(name, " Forecast: ", label))
    
    polygon(c(test_yr_ext, rev(test_yr_ext)),
            c(lo95_ext, rev(hi95_ext)),
            col = adjustcolor("red", alpha.f = 0.15), border = NA)
    polygon(c(test_yr_ext, rev(test_yr_ext)),
            c(lo80_ext, rev(hi80_ext)),
            col = adjustcolor("red", alpha.f = 0.25), border = NA)
    
    lines(train_yr,    tr,     col = "black",  lwd = 1.5)
    lines(test_yr_ext, te_ext, col = "grey40", lwd = 1.5)
    lines(test_yr_ext, fc_ext, col = "red",    lwd = 2, lty = 2)
    
    abline(v = 2012, lty = 3, col = "grey60")
    
    legend("topright",
           legend = c("Training (1960-2012)", "Test actual (2013-2024)",
                      "Point forecast", "80% PI", "95% PI"),
           col = c("black", "grey40", "red",
                   adjustcolor("red", alpha.f = 0.4),
                   adjustcolor("red", alpha.f = 0.2)),
           lty = c(1, 1, 2, 1, 1),
           lwd = c(1.5, 1.5, 2, 8, 8),
           bty = "n", cex = 0.8)
  }
  
  png(filename, width = 1800, height = 1000, res = 200)
  draw()
  dev.off()
  
  draw()
}

if (!is.null(best_tfr)) plot_forecast(TFR_train, TFR_test, best_tfr$fit,
                                      best_tfr$label, "TFR",
                                      "figures/fig7_tfr_forecast.png")

# TFR forecast looks pretty good.
# Red forecast line stays close to the grey actual line.
# Actual values look like they stay inside the 95% band.
# By 2024, the model still follows the downward trend.
# Overall, ARIMA(13,2,2) seems to work well for TFR.

if (!is.null(best_tlb)) plot_forecast(TLB_train, TLB_test, best_tlb$fit,
                                      best_tlb$label, "TLB",
                                      "figures/fig8_tlb_forecast.png")

# TLB forecast is a bit messier than TFR.
# Red forecast line is below the actual values for most years.
# Model expects a jump near 2024 because of the Dragon year cycle.
# Actual TLB drops instead, so the point forecast misses the direction.
# Actual values still look like they sit inside the wider 95% band.

#Saveing results for FinalReport.Rmd 

saveRDS(tfr_results, "tfr_results.rds")
saveRDS(tlb_results, "tlb_results.rds")
saveRDS(tfr_table,   "tfr_table.rds")
saveRDS(tlb_table,   "tlb_table.rds")
saveRDS(best_tfr,    "best_tfr.rds")
saveRDS(best_tlb,    "best_tlb.rds")

