# Author:  Kiran Shankaran
# ============================================================================

# Updated EDA. Feed back from the first EDA submission was taken into
# consideration and the EDA was redone

# tseries provides kpss.test();
pacman::p_load(tseries)

# Create figures folder if it doesn't exist
if (!dir.exists("figures")) dir.create("figures")



# LOAD DATA
#----------------------------------------------------------------------------

# Data including nationality, age and etc was downloaded this time
raw <- read.csv("BirthsAndFertilityRates.csv",
                na.strings = "na",
                stringsAsFactors = FALSE,
                check.names = FALSE)

# Trim whitespace from variable labels 
labels_raw <- trimws(raw[, 1])

# Helper: extract a named row by exact match, convert to numeric, reverse 
# year order (CSV runs 2025..1960), drop incomplete 2025 column.
extract_row <- function(name) {
  matches <- labels_raw == name
  if (sum(matches) != 1) {
    stop(sprintf("Label '%s' matched %d rows (expected 1)",
                 name, sum(matches)))
  }
  row_vals <- raw[matches, -1]
  values <- as.numeric(unlist(row_vals))
  values <- rev(values)
  values[1:65]
}

# Build the dataframe. Match the full row labels exactly as they appear 
# in the CSV 
fertility <- data.frame(
  Year         = 1960:2024,
  TFR          = extract_row("Total Fertility Rate (TFR) (Per Female)"),
  TLB          = extract_row("Total Live-Births (Number)"),
  TFR_Chinese  = extract_row("Chinese (Per Female)"),
  TFR_Malay    = extract_row("Malays (Per Female)"),
  TFR_Indian   = extract_row("Indians (Per Female)")
)

# Quick sanity check 
cat("TLB 1960: ", fertility$TLB[1],  " (expected 61775)\n", sep = "")
cat("TLB 2024: ", fertility$TLB[65], " (expected 33703)\n", sep = "")

# From the suggested reading from the feedback, it appears that there is some
# relationship between the chinese zodiac years and the birth rate in singapore.
# Literature shows The birth rate is significantly higher during the Year of the
# Dragon compared to the Year of the Tiger because of deep-rooted cultural
# superstitions and astrological preferences in many Asian societies. So the
# zodiac signs are considered in this analysis, as previously the stop-at

# CHINESE ZODIAC YEAR LABELS
# ----------------------------------------------------------------------------

# The Chinese zodiac cycles every 12 years, with Rat at 1900. 
# (Year - 1900) mod 12 gives the position in the cycle.
zodiac_names <- c("Rat", "Ox", "Tiger", "Rabbit",
                  "Dragon", "Snake", "Horse", "Goat",
                  "Monkey", "Rooster", "Dog", "Pig")
fertility$Zodiac <- zodiac_names[((fertility$Year - 1900) %% 12) + 1]

# Convenience vectors for use in plotting
dragon_years <- fertility$Year[fertility$Zodiac == "Dragon"]
tiger_years  <- fertility$Year[fertility$Zodiac == "Tiger"]

cat("Dragon years (expected boost in Chinese births): ", dragon_years, "\n")
cat("Tiger years  (expected dip in Chinese births):   ", tiger_years,  "\n\n")



# TIME SERIES OBJECTS AND TRAIN/TEST SPLIT
# ----------------------------------------------------------------------------

# Base R's arima() function requires ts objects.
TLB <- ts(fertility$TLB, start = 1960, frequency = 1)
TFR <- ts(fertility$TFR, start = 1960, frequency = 1)

# Training period 1960-2012 (53 years); test 2013-2024 (12 years)
TLB_train <- window(TLB, end = 2012)
TFR_train <- window(TFR, end = 2012)
TLB_test  <- window(TLB, start = 2013)
TFR_test  <- window(TFR, start = 2013)



# 5. TIME PLOTS OF RAW SERIES WITH DRAGON / TIGER MARKERS
# ----------------------------------------------------------------------------

# Figure 1: Raw TFR. Dragon years as red circles, Tiger years as blue 
# triangles, train/test boundary as dashed vertical line at 2012.
plot_tfr_raw <- function() {
  plot(fertility$Year, fertility$TFR, type = "l", lwd = 1.5,
       xlab = "Year", ylab = "TFR",
       main = "Singapore Total Fertility Rate, 1960-2024")
  points(dragon_years,
         fertility$TFR[fertility$Year %in% dragon_years],
         pch = 19, col = "red", cex = 1.4)
  points(tiger_years,
         fertility$TFR[fertility$Year %in% tiger_years],
         pch = 17, col = "blue", cex = 1.4)
  abline(v = 2012, lty = 2, col = "grey60")
  legend("topright",
         legend = c("Dragon years", "Tiger years", "Train/test boundary"),
         pch = c(19, 17, NA),
         lty = c(NA, NA, 2),
         col = c("red", "blue", "grey60"),
         bty = "n", cex = 0.9)
}

png("figures/fig1_tfr_raw.png", width = 1800, height = 1000, res = 200)
plot_tfr_raw()
dev.off()
plot_tfr_raw()

# From the graph it seems like there are spikes during the year of the dragon in
# the TFR plot which reflects what was found in the literature

# Figure 2: Raw TLB with same markers
plot_tlb_raw <- function() {
  plot(fertility$Year, fertility$TLB, type = "l", lwd = 1.5,
       xlab = "Year", ylab = "TLB",
       main = "Singapore Total Live Births, 1960-2024")
  points(dragon_years,
         fertility$TLB[fertility$Year %in% dragon_years],
         pch = 19, col = "red", cex = 1.4)
  points(tiger_years,
         fertility$TLB[fertility$Year %in% tiger_years],
         pch = 17, col = "blue", cex = 1.4)
  abline(v = 2012, lty = 2, col = "grey60")
  legend("topright",
         legend = c("Dragon years", "Tiger years", "Train/test boundary"),
         pch = c(19, 17, NA),
         lty = c(NA, NA, 2),
         col = c("red", "blue", "grey60"),
         bty = "n", cex = 0.9)
}

png("figures/fig2_tlb_raw.png", width = 1800, height = 1000, res = 200)
plot_tlb_raw()
dev.off()
plot_tlb_raw()

# Same trend can be seen where the dragon year has a spike in births and there
# is a bit of a dip during the year of the tiger.

# The raw TFR and TLB plots show a strong long-term decline, especially from
# 1960 to the early 1980s. Dragon years appear to coincide with local increases
# in fertility and births, particularly in TLB. The effect is most visually
# noticeable around 1988 and 2000. Tiger years appear to show some weaker dips,
# but this pattern is less consistent.

# TRANSFORMATIONS
# ----------------------------------------------------------------------------

# Log transform to stabilise variance
TFR_log <- log(TFR_train)
TLB_log <- log(TLB_train)

# Seasonal differencing at lag 12 to remove the Dragon Year cycle
TFR_log_sdiff <- diff(TFR_log, lag = 12)
TLB_log_sdiff <- diff(TLB_log, lag = 12)

# First-order regular differencing on top to remove the trend
TFR_log_sdiff_rdiff <- diff(TFR_log_sdiff, lag = 1)
TLB_log_sdiff_rdiff <- diff(TLB_log_sdiff, lag = 1)


# TIME PLOTS OF TRANSFORMATION STAGES
# ----------------------------------------------------------------------------

# Figure 5: TFR transformation stages
plot_tfr_transformations <- function() {
  par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
  
  plot(TFR_train, type = "l", lwd = 1.5,
       main = "Raw TFR (training)",
       xlab = "Year", ylab = "TFR")
  
  plot(TFR_log, type = "l", lwd = 1.5,
       main = "log(TFR)",
       xlab = "Year", ylab = "log(TFR)")
  
  plot(TFR_log_sdiff, type = "l", lwd = 1.5,
       main = "log(TFR) with seasonal diff (lag 12)",
       xlab = "Year", ylab = "")
  abline(h = 0, lty = 3, col = "grey50")
  
  plot(TFR_log_sdiff_rdiff, type = "l", lwd = 1.5,
       main = "log(TFR) with seasonal + regular diff",
       xlab = "Year", ylab = "")
  abline(h = 0, lty = 3, col = "grey50")
  
  par(mfrow = c(1, 1))
}

png("figures/fig5_tfr_transformations.png", width = 1800, height = 1400, res = 200)
plot_tfr_transformations()
dev.off()
plot_tfr_transformations()


# Figure 6: TLB transformation stages
plot_tlb_transformations <- function() {
  par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
  
  plot(TLB_train, type = "l", lwd = 1.5,
       main = "Raw TLB (training)",
       xlab = "Year", ylab = "TLB")
  
  plot(TLB_log, type = "l", lwd = 1.5,
       main = "log(TLB)",
       xlab = "Year", ylab = "log(TLB)")
  
  plot(TLB_log_sdiff, type = "l", lwd = 1.5,
       main = "log(TLB) with seasonal diff (lag 12)",
       xlab = "Year", ylab = "")
  abline(h = 0, lty = 3, col = "grey50")
  
  plot(TLB_log_sdiff_rdiff, type = "l", lwd = 1.5,
       main = "log(TLB) with seasonal + regular diff",
       xlab = "Year", ylab = "")
  abline(h = 0, lty = 3, col = "grey50")
  
  par(mfrow = c(1, 1))
}

png("figures/fig6_tlb_transformations.png", width = 1800, height = 1400, res = 200)
plot_tlb_transformations()
dev.off()
plot_tlb_transformations()


# KPSS STATIONARITY TESTS AT EACH STAGE
# ----------------------------------------------------------------------------

# Helper: run KPSS, suppress the truncation warning, return named vector.
run_kpss <- function(series) {
  out <- suppressWarnings(kpss.test(series, null = "Level"))
  c(KPSS_statistic = unname(round(out$statistic, 4)),
    p_value        = round(out$p.value, 4))
}

# Run KPSS at four stages for each series: raw, log, log + seasonal diff, 
# log + seasonal + regular diff. Small p-value -> reject stationarity.
kpss_table <- rbind(
  "TFR raw"                            = run_kpss(TFR_train),
  "log(TFR)"                           = run_kpss(TFR_log),
  "log(TFR) seasonal diff(12)"         = run_kpss(TFR_log_sdiff),
  "log(TFR) seasonal + regular diff"   = run_kpss(TFR_log_sdiff_rdiff),
  "TLB raw"                            = run_kpss(TLB_train),
  "log(TLB)"                           = run_kpss(TLB_log),
  "log(TLB) seasonal diff(12)"         = run_kpss(TLB_log_sdiff),
  "log(TLB) seasonal + regular diff"   = run_kpss(TLB_log_sdiff_rdiff)
)


print(kpss_table)


# The seasonal + regular differencing establishes one route to 
# stationarity, supporting SARIMA models. We also test whether second-order 
# regular differencing alone achieves stationarity, which would support 
# non-seasonal ARIMA models with d=2.

TFR_log_d2 <- diff(diff(log(TFR_train)))
TLB_log_d2 <- diff(diff(log(TLB_train)))

kpss_d2 <- rbind(
  "log(TFR) twice-differenced (d=2)" = run_kpss(TFR_log_d2),
  "log(TLB) twice-differenced (d=2)" = run_kpss(TLB_log_d2)
)

cat("\n=== KPSS for alternative d=2 route ===\n")
print(kpss_d2)

# KPSS confirms that the raw and log-transformed series are non-stationary. For
# TFR, lag-12 seasonal differencing alone is not enough, so both seasonal and
# regular differencing are needed. For TLB, lag-12 seasonal differencing appears
# sufficient, so regular differencing may not be required. This suggests using
# d=1,D=1 for TFR, but considering d=0,D=1 for TLB.

# KPSS also shows that a d=2 route makes both logged series stationary.
# However, this is a more aggressive non-seasonal approach and does not directly
# model the 12-year Dragon-cycle structure, so it is treated as a comparison.


# ACF AND PACF OF STATIONARY SERIES (lag.max = 36)
# ----------------------------------------------------------------------------

# After log + seasonal + regular differencing, both series should be 
# stationary. ACF and PACF identify candidate ARIMA orders. lag.max = 36 
# covers three full 12-year cycles so seasonal harmonics at 12, 24, 36 
# are visible.

# Figure 3: ACF and PACF for stationary log(TFR)
plot_tfr_acf_pacf <- function() {
  par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
  acf(TFR_log_sdiff_rdiff, lag.max = 36,
      main = "ACF: log(TFR), seasonal + regular diff")
  pacf(TFR_log_sdiff_rdiff, lag.max = 36,
       main = "PACF: log(TFR), seasonal + regular diff")
  par(mfrow = c(1, 1))
}

png("figures/fig3_tfr_acf_pacf.png", width = 1800, height = 900, res = 200)
plot_tfr_acf_pacf()
dev.off()
plot_tfr_acf_pacf()

# Figure 4: ACF and PACF for stationary log(TLB)
plot_tlb_acf_pacf <- function() {
  par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
  acf(TLB_log_sdiff_rdiff, lag.max = 36,
      main = "ACF: log(TLB), seasonal + regular diff")
  pacf(TLB_log_sdiff_rdiff, lag.max = 36,
       main = "PACF: log(TLB), seasonal + regular diff")
  par(mfrow = c(1, 1))
}

png("figures/fig4_tlb_acf_pacf.png", width = 1800, height = 900, res = 200)
plot_tlb_acf_pacf()
dev.off()
plot_tlb_acf_pacf()

# ACF/PACF of LOG SERIES WITH FIRST DIFFERENCE ONLY
# ----------------------------------------------------------------------------

# The lecturer's feedback referred specifically to "log or 1st order 
# difference" — a simpler transformation than the seasonal + regular 
# differencing applied above. ACF/PACF of this simpler transformation 
# should show the most prominent lag 11, 12, 13 spikes, which directly 
# motivate the inclusion of a seasonal component.

TFR_log_rdiff_only <- diff(log(TFR_train), lag = 1)
TLB_log_rdiff_only <- diff(log(TLB_train), lag = 1)

# Figure 7: ACF/PACF of first-differenced log TFR
plot_tfr_acf_pacf_first <- function() {
  par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
  acf (TFR_log_rdiff_only, lag.max = 36,
       main = "ACF: log(TFR), first diff only")
  pacf(TFR_log_rdiff_only, lag.max = 36,
       main = "PACF: log(TFR), first diff only")
  par(mfrow = c(1, 1))
}

png("figures/fig7_tfr_acf_pacf_first_diff.png", width = 1800, height = 900, res = 200)
plot_tfr_acf_pacf_first()
dev.off()
plot_tfr_acf_pacf_first()

# Figure 8: ACF/PACF of first-differenced log TLB
plot_tlb_acf_pacf_first <- function() {
  par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
  acf (TLB_log_rdiff_only, lag.max = 36,
       main = "ACF: log(TLB), first diff only")
  pacf(TLB_log_rdiff_only, lag.max = 36,
       main = "PACF: log(TLB), first diff only")
  par(mfrow = c(1, 1))
}

png("figures/fig8_tlb_acf_pacf_first_diff.png", width = 1800, height = 900, res = 200)
plot_tlb_acf_pacf_first()
dev.off()
plot_tlb_acf_pacf_first()

# Print values at lags 1-20 for inspection
tfr_acf2  <- acf (TFR_log_rdiff_only, lag.max = 36, plot = FALSE)
tfr_pacf2 <- pacf(TFR_log_rdiff_only, lag.max = 36, plot = FALSE)
tlb_acf2  <- acf (TLB_log_rdiff_only, lag.max = 36, plot = FALSE)
tlb_pacf2 <- pacf(TLB_log_rdiff_only, lag.max = 36, plot = FALSE)

acf_pacf_table_first_diff <- data.frame(
  Lag      = 1:20,
  TFR_ACF  = round(tfr_acf2$acf [2:21], 3),
  TFR_PACF = round(tfr_pacf2$acf[1:20], 3),
  TLB_ACF  = round(tlb_acf2$acf [2:21], 3),
  TLB_PACF = round(tlb_pacf2$acf[1:20], 3)
)

n_eff2 <- length(TFR_log_rdiff_only)
print(acf_pacf_table_first_diff)

# The first-differenced log series still shows significant autocorrelation
# around lag 12. This supports the presence of a 12-year cyclic structure. Since
# the Chinese zodiac repeats every 12 years, this supports testing SARIMA models
# with s=12. This does not prove the Dragon-year effect by itself, but it
# justifies including a 12-year model component.


# PRINT ACF/PACF VALUES FOR LAG IDENTIFICATION
# ----------------------------------------------------------------------------

# Print actual ACF and PACF values at each lag so you can identify which 
# lags exceed the significance bound +/-2/sqrt(n).
n_eff <- length(TFR_log_sdiff_rdiff)
sig_bound <- 2 / sqrt(n_eff)

cat("\n=== ACF/PACF significance bound ===\n")
cat("Effective sample size after differencing: ", n_eff, "\n", sep = "")
cat("Significance bound: +/-", round(sig_bound, 3), "\n", sep = "")
cat("(values outside this band are significant at the 5% level)\n\n")

# Pull ACF/PACF objects without plotting
tfr_acf  <- acf (TFR_log_sdiff_rdiff, lag.max = 36, plot = FALSE)
tfr_pacf <- pacf(TFR_log_sdiff_rdiff, lag.max = 36, plot = FALSE)
tlb_acf  <- acf (TLB_log_sdiff_rdiff, lag.max = 36, plot = FALSE)
tlb_pacf <- pacf(TLB_log_sdiff_rdiff, lag.max = 36, plot = FALSE)

# Build a table of lags 1-20. acf$acf[1] is lag 0 = 1, so index from [2].
# pacf$acf[1] is already lag 1.
acf_pacf_table <- data.frame(
  Lag      = 1:20,
  TFR_ACF  = round(tfr_acf$acf [2:21], 3),
  TFR_PACF = round(tfr_pacf$acf[1:20], 3),
  TLB_ACF  = round(tlb_acf$acf [2:21], 3),
  TLB_PACF = round(tlb_pacf$acf[1:20], 3)
)


cat("(values exceeding +/-", round(sig_bound, 3),
    " are significant)\n\n", sep = "")
print(acf_pacf_table)

# After seasonal and regular differencing, most ACF/PACF values are within the
# significance bounds, suggesting the series is closer to stationary.
# Significant lag-1 values suggest testing low-order AR or MA terms. Remaining
# lag-12 or lag-13 behaviour suggests testing seasonal/cyclic SARIMA terms.
# Therefore, candidate models should include low-order p and q, and a 12-year
# seasonal component.


# SAVE TABLES FOR FINAL REPORT
# ----------------------------------------------------------------------------

saveRDS(kpss_table,                "kpss_table.rds")
saveRDS(kpss_d2,                   "kpss_d2.rds")
saveRDS(acf_pacf_table,            "acf_pacf_table.rds")
saveRDS(acf_pacf_table_first_diff, "acf_pacf_table_first_diff.rds")

