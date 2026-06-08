# Time Series Project: Singapore Births Forecast
This repository contains the final time series analysis for forecasting Singapore's total live birth (TBL) and total fertility rate (TFR) from 1960 to 2024.

This project investigates whether the Chinese zodiac Dragon Year effect remains a useful forecasting signal for Singapore's birth data, especially given that the 2024 Dragon Year did not produce an obvious increase in births.

Author: Kiran Shankaran

Course: Time Series Analysis 
Adelaide University

## Project Overview
Singapore's TFR and TLB have declined substantially since 1960. At the same time, previous research has found that Chinese births in Singapore have historically increased in Dragon years and decreased in Tiger years. Since the Chinese zodiac cycle repeats every 12 years, this project tests whether a 12-year seasonal structure improves forecasts.

The analysis compares non-seasonal ARIMA models against SARIMA models with a 12-year seasonal period. Models are fitted on the 1960--2012 training period and evaluated on the 2013--2024 test period. Model viability is checked using the Ljung-Box test, and forecast accuracy is compared using RMSE.

## Research Question

Does the historical Chinese zodiac Dragon Year effect remain a useful forecasting signal for Singapore's Total Live Births and Total Fertility Rate, given that the 2024 Dragon Year produced no observable birth boost?

## Repository Contents

| File | Description |
|---|---|
| `README.md` | Project overview and file guide |
| `BirthsAndFertilityRates.csv` | Raw SingStat data used for the final report |
| `Analysis.R` | Main R script used for model fitting, diagnostics, forecasts and plots |
| `FinalReport.Rmd` | R Markdown source file for the final report |
| `EDA.Rmd` | Earlier exploratory analysis file |
| `EDA updated from feedback.R` | Updated exploratory code and plot/table saving |


## Data Source

The data was obtained from the Singapore Department of Statistics:  
https://www.singstat.gov.sg



