---
title: "Francesco Sanna"
output: pdf_document
---

# Summary of my Bachelor Thesis: Environmental Kuznets Curve & Spatial Econometrics

## 1. Project Objective and Implementation

My thesis investigates the Environmental Kuznets Curve (EKC) across 213 European NUTS 2 regions (2014–2022), replicating and extending David Maddison's (2006) spatial econometric framework. The empirical execution employs a First-Difference ($\Delta \ln$) panel estimator and econometric modelling, integrated then with spatial weight matrices ($W$) to model transboundary spillovers, moving to areal data modelling and spatial econometrics. We process and harmonize multi-source data from EDGAR v8.0 and Eurostat to evaluate whether regional economic growth decouples from emissions ($CO_2, CH_4, N_2O, \text{F-gases}$).

## 2. Data Wrangling

- **EDGAR v8.0 (EC-JRC):** NUTS2 emission inventories ($ktons$) across all economic sectors for Fossil $CO_2$, $CH_4$, $N_2O$, and $\text{F-gases}$ ($\text{GWP}_{100}\text{ AR5}$).

- **Eurostat NAMA & DEMO:** Regional $\text{GDP}$ and total resident population.

- **Eurostat GISCO:** $\text{EPSG:3035}$.

Made dataset tidy, found the per-capita value, performed first difference and joint the dataset.

### Final Balanced Dataset Sample ($N=213$ Regions, 852 Observations)

| $\text{NUTS ID}$ | $\text{Subst}$ | $E_{2014}$ | $E_{2022}$ | $\text{GDP}_{2014}$ | $\text{GDP}_{2022}$ | $\Delta \ln(E)$ | $\Delta \ln(\text{GDP})$ | $\Delta (\ln(\text{GDP})^2)$ |
|:------:|:------:|:------:|:------:|:------:|:------:|:------:|:------:|:------:|
| **AT11** | $CO_2$ | $6.593 \times 10^{-3}$ | $6.034 \times 10^{-3}$ | $26,608$ | $34,771$ | $-0.0886$ | $+0.2676$ | $+5.5240$ |
| **AT11** | $CH_4$ | $6.257 \times 10^{-4}$ | $6.087 \times 10^{-4}$ | $26,608$ | $34,771$ | $-0.0275$ | $+0.2676$ | $+5.5240$ |
| **AT11** | $\text{F-gases}$ | $2.562 \times 10^{-4}$ | $2.982 \times 10^{-4}$ | $26,608$ | $34,771$ | $+0.1520$ | $+0.2676$ | $+5.5240$ |
| **AT11** | $N_2O$ | $7.504 \times 10^{-4}$ | $6.677 \times 10^{-4}$ | $26,608$ | $34,771$ | $-0.1168$ | $+0.2676$ | $+5.5240$ |
| **AT12** | $CO_2$ | $1.030 \times 10^{-2}$ | $8.764 \times 10^{-3}$ | $32,016$ | $41,036$ | $-0.1616$ | $+0.2482$ | $+5.2119$ |
| **AT12** | $CH_4$ | $1.371 \times 10^{-3}$ | $1.238 \times 10^{-3}$ | $32,016$ | $41,036$ | $-0.1020$ | $+0.2482$ | $+5.2119$ |
| **AT12** | $\text{F-gases}$ | $3.023 \times 10^{-4}$ | $3.439 \times 10^{-4}$ | $32,016$ | $41,036$ | $+0.1290$ | $+0.2482$ | $+5.2119$ |
| **AT12** | $N_2O$ | $6.734 \times 10^{-4}$ | $5.940 \times 10^{-4}$ | $32,016$ | $41,036$ | $-0.1254$ | $+0.2482$ | $+5.2119$ |

## 3. Econometric Modeling & Diagnostic Testing

Estimation of the model via First Differences using `tidyverse`, `lmtest`, and `tseries`: $$\Delta \ln(E_i) = \beta_1 \Delta \ln(\text{GDP}_i) + \beta_2 \Delta (\ln(\text{GDP}_i)^2) + \Delta \varepsilon_i$$

I implemented the OLS estimation across pollutants and evaluated residual validity against Gauss-Markov assumptions via Jarque-Bera normality tests and White heteroskedasticity tests:

``` r
# Libraries: tidyverse, lmtest, tseries
# OLS estimation in first differences & specification tests per substance
model   <- lm(diff_ln_E ~ diff_ln_GDP + diff_ln_GDP_sq, data = data_sub)
jb_p    <- jarque.bera.test(residuals(model))$p.value
white_p <- bptest(model, ~ diff_ln_GDP * diff_ln_GDP_sq + 
              I(diff_ln_GDP^2) + I(diff_ln_GDP_sq^2), data = data_sub)$p.value
```

While the inverted-U is confirmed ($\beta_1 > 0, \beta_2 < 0$), all specifications reject error normality ($p < 0.05$), indicating structural omission of dependence.


## 4. Spatial Weight & Geometric Modeling

Using `sf` and `spdep` package, I mapped the spatial structure by calculating geometric centroids from official NUTS 2 polygons and constructed three row-standardized spatial weight matrices $W$ ($w_{ij}^* = w_{ij}/\sum_j w_{ij}$): 1. **Inverse Distance with Critical Threshold** ($378\text{ km}$ = 1 -NN max distance): 2. **Inverse Distance Unlimited:** Dense connectivity modeling long-range decay: $w_{ij} = 1/d_{ij}$. 3. $k$**-Nearest Neighbors** ($k=5$): Uniform spatial connectivity ($w_{ij}^* = 0.20$).

``` r
# Libraries: sf, spdep
# 1st Weigth matrix 
coords <- st_coordinates(st_centroid(st_geometry(sf_all %>% filter(Substance == "CO2")),
                        of_largest_polygon = TRUE))
max_1nn <- max(unlist(nbdists(knn2nb(knearneigh(coords, k = 1)), coords))) #k-NN threshold
nb_threshold <- dnearneigh(coords, d1 = 0, d2 = max_1nn)
idw_threshold <- lapply(nbdists(nb_threshold, coords), function(x) 1/x)
W_inv_dist_threshold <- nb2listw(nb_threshold,
                        glist = idw_threshold, style = "W", zero.policy = TRUE)
```

## 5. Spatial test and estimation

Test for spatial autocorrelation and select the appropriate model (SAR vs. SEM). First I performed Moran's $I$ on OLS residuals and then Lagrange Multiplier tests (`lm.LMtests`, including Robust LM). Models were then estimated via Maximum Likelihood and FGLS using `spatialreg`:

``` r
# Libraries: spdep, spatialreg
# 1. Spatial diagnostics on regression residuals
morans_lm <- lm.morantest(ekc_model, W_knn5_rowstd)
morans_lm <- lm.morantest(ekc_model, W_inv_dist_unlimited)
morans_lm <- lm.morantest(ekc_model, W_inv_dist_threshold)
LM_sem    <- lm.LMtests(ekc_model, W_knn5_rowstd, test = c("LMerr", "RLMerr"))
LM_sar    <- lm.LMtests(ekc_model, W_inv_dist_unlimited, test = c("LMlag", "RLMlag"))

# 2. Maximum Likelihood / FGLS spatial estimation
sem_Fgases <- errorsarlm(diff_ln_E ~ diff_ln_GDP + diff_ln_GDP_sq, 
                                  data = df_Fgases, listw = W_knn5_rowstd)
sar_N2O    <- lagsarlm(diff_ln_E ~ diff_ln_GDP + diff_ln_GDP_sq, 
                                  data = df_N2O, listw = W_inv_dist_unlimited)
```

### Comparative Estimation Results

| Model / Substance | $\Delta \ln(\text{GDP})$ [p-val] | $\Delta (\ln\text{GDP})^2$ [p-val] | Spatial Parameter [p-val] |
|:-----------------|:----------------:|:----------------:|:----------------:|
| **OLS:** $CO_2$ | $+3.170$ [$<0.001$] | $-0.153$ [$<0.001$] | — |
| **OLS:** $CH_4$ | $+2.105$ [$<0.001$] | $-0.092$ [$0.003$] | — |
| **SAR:** $N_2O$ *(Unlimited* $W$) | $+2.132$ [$<0.001$] | $-0.100$ [$<0.001$] | $\boldsymbol{\rho = 0.931}$ [$<0.001$] |
| **SEM:** $\text{F-gases}$ *(*$k\text{-NN } 5$) | $+2.271$ [$0.052$] | $-0.058$ [$0.339$] | $\boldsymbol{\lambda = 0.865}$ [$<0.001$] |
