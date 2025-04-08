# Modeling Climate-Driven Malaria Transmission in Zambia

This R project simulates mosquito population dynamics and malaria transmission in Zambia under various climate scenarios using temperature and rainfall data.

## How It Works

The simulation uses:
- Entomological lifecycle models from White et al. (2011)
- Climate-sensitive parameters (e.g., oviposition, mortality)
- Multiple scenarios: historical climatology, RCP 1.9, 4.5, 8.5

##  Project Structure

- `R/AIDM_Project.R` — Main simulation code
- `data/student_met_data.rds` — Rainfall and temperature data

## Required R Packages

```r
library(deSolve)
library(zoo)
library(tidyverse)
library(dplyr)
library(gridExtra)
library(ggplot2)
library(tidyr)
