---
title: "eDNAModel: gllvm"
author: "Marwah Soliman, Bert van der Veen"
date: "2026-09-22"
output: 
  rmarkdown::html_vignette:
    toc: true
    toc_depth: 3
    number_sections: true
vignette: >
  %\VignetteIndexEntry{eDNAModel: gllvm}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---




``` r
library(phyloseq)
library(Matrix)
library(eDNAModel)
library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(tibble)
library(gllvm)
library(glmmTMB)
library(reshape2)
library(grid)

```

# 1. Introduction

This vignette demonstrates the updated occupancy–detection model workflow using the `eDNAModel` package and `phyloseq` data objects. The workflow includes fitting binomial and Poisson mixed-effects models across iterations, visualizing occupancy and detection probabilities, using gllvm for occupancy model.

---

# 2. Load and Subset Data

``` r
physeq_new=readRDS("~/Downloads/BGE_Port_sampling.rds")

data("physeq_new",package = "eDNAModel")

physeq_one=physeq_new$`Marine Invasive species Trondheim`
```

# 3. Fit the Occupancy–Detection Model


```
#> No abundance offset detected.
#> OTUs before filtering: 43
#> OTUs after filtering: 42
#> Iteration 1
#> Warning in se.gllvm(out): 3 parameter(s) have negative variance estimates (b x3). The model likely
#> has not converged - consider re-fitting.
#> Warning in finalizeTMB(TMBStruc, obj, fit, h, data.tmb.old): Model convergence problem; false
#> convergence (8). See vignette('troubleshooting'), help('diagnose')
#> Iteration 1 finished in 1.27 seconds
#> Iteration 2
#> Warning in se.gllvm(out): 3 parameter(s) have negative variance estimates (b x3). The model likely
#> has not converged - consider re-fitting.
#> Iteration 2 finished in 1.29 seconds
#> Iteration 3
#> Warning in finalizeTMB(TMBStruc, obj, fit, h, data.tmb.old): Model convergence problem; false
#> convergence (8). See vignette('troubleshooting'), help('diagnose')
#> Iteration 3 finished in 1.66 seconds
#> Iteration 4
#> Warning in finalizeTMB(TMBStruc, obj, fit, h, data.tmb.old): Model convergence problem; false
#> convergence (8). See vignette('troubleshooting'), help('diagnose')
#> Iteration 4 finished in 1.56 seconds
#> Iteration 5
#> Warning in se.gllvm(out): 6 parameter(s) have negative variance estimates (b x6). The model likely
#> has not converged - consider re-fitting.
#> Warning in se.gllvm(out): Model convergence problem; false convergence (8). See
#> vignette('troubleshooting'), help('diagnose')
#> Iteration 5 finished in 1.31 seconds
#> Iteration 6
#> Warning in se.gllvm(out): 53 parameter(s) have negative variance estimates (b x16, lambda x36,
#> sigmaLV x1). The model likely has not converged - consider re-fitting.
#> Warning in se.gllvm(out): Model convergence problem; false convergence (8). See
#> vignette('troubleshooting'), help('diagnose')
#> Iteration 6 finished in 1.84 seconds
#> Iteration 7
#> Warning in se.gllvm(out): 3 parameter(s) have negative variance estimates (b x3). The model likely
#> has not converged - consider re-fitting.
#> Iteration 7 finished in 1.27 seconds
#> Iteration 8
#> Warning in se.gllvm(out): 21 parameter(s) have negative variance estimates (b x3, lambda x17,
#> sigmaLV x1). The model likely has not converged - consider re-fitting.
#> Iteration 8 finished in 1.25 seconds
#> Iteration 9
#> Warning in se.gllvm(out): 3 parameter(s) have negative variance estimates (b x3). The model likely
#> has not converged - consider re-fitting.
#> Iteration 9 finished in 1.3 seconds
#> Iteration 10
#> Warning in se.gllvm(out): 12 parameter(s) have negative variance estimates (b x7, lambda x4,
#> sigmaLV x1). The model likely has not converged - consider re-fitting.
#> Iteration 10 finished in 1.35 seconds
#> Iteration 11
#> Warning in se.gllvm(out): 6 parameter(s) have negative variance estimates (b x6). The model likely
#> has not converged - consider re-fitting.
#> Iteration 11 finished in 2.2 seconds
#> Iteration 12
#> Iteration 12 finished in 1.56 seconds
#> Iteration 13
#> Warning in se.gllvm(out): 12 parameter(s) have negative variance estimates (b x12). The model
#> likely has not converged - consider re-fitting.
#> Iteration 13 finished in 2.12 seconds
#> Iteration 14
#> Warning in se.gllvm(out): 11 parameter(s) have negative variance estimates (b x7, lambda x3,
#> sigmaLV x1). The model likely has not converged - consider re-fitting.
#> Iteration 14 finished in 1.31 seconds
#> Iteration 15
#> Iteration 15 finished in 1.73 seconds
#> Iteration 16
#> Warning in se.gllvm(out): 28 parameter(s) have negative variance estimates (b x4, lambda x23,
#> sigmaLV x1). The model likely has not converged - consider re-fitting.
#> Iteration 16 finished in 1.33 seconds
#> Iteration 17
#> Warning in se.gllvm(out): 57 parameter(s) have negative variance estimates (b x18, lambda x38,
#> sigmaLV x1). The model likely has not converged - consider re-fitting.
#> Iteration 17 finished in 2.17 seconds
#> Iteration 18
#> Warning in se.gllvm(out): 9 parameter(s) have negative variance estimates (b x5, lambda x3, sigmaLV
#> x1). The model likely has not converged - consider re-fitting.
#> Iteration 18 finished in 1.24 seconds
#> Iteration 19
#> Warning in se.gllvm(out): 7 parameter(s) have negative variance estimates (b x2, lambda x4, sigmaLV
#> x1). The model likely has not converged - consider re-fitting.
#> Iteration 19 finished in 1.56 seconds
#> Iteration 20
#> Warning in se.gllvm(out): 3 parameter(s) have negative variance estimates (b x3). The model likely
#> has not converged - consider re-fitting.
#> Iteration 20 finished in 1.29 seconds
#> Iteration 21
#> Warning in se.gllvm(out): 4 parameter(s) have negative variance estimates (b x4). The model likely
#> has not converged - consider re-fitting.
#> Iteration 21 finished in 1.3 seconds
#> Iteration 22
#> Iteration 22 finished in 1.63 seconds
#> Iteration 23
#> Warning in se.gllvm(out): 5 parameter(s) have negative variance estimates (b x5). The model likely
#> has not converged - consider re-fitting.
#> Iteration 23 finished in 1.89 seconds
#> Iteration 24
#> Warning in se.gllvm(out): 39 parameter(s) have negative variance estimates (b x14, lambda x24,
#> sigmaLV x1). The model likely has not converged - consider re-fitting.
#> Iteration 24 finished in 2.2 seconds
#> Iteration 25
#> Iteration 25 finished in 1.64 seconds
#> Iteration 26
#> Warning in se.gllvm(out): 18 parameter(s) have negative variance estimates (b x6, lambda x11,
#> sigmaLV x1). The model likely has not converged - consider re-fitting.
#> Iteration 26 finished in 1.27 seconds
#> Iteration 27
#> Warning in se.gllvm(out): 8 parameter(s) have negative variance estimates (b x8). The model likely
#> has not converged - consider re-fitting.
#> Iteration 27 finished in 1.3 seconds
#> Iteration 28
#> Iteration 28 finished in 1.87 seconds
#> Iteration 29
#> Warning in se.gllvm(out): 9 parameter(s) have negative variance estimates (b x5, lambda x3, sigmaLV
#> x1). The model likely has not converged - consider re-fitting.
#> Iteration 29 finished in 1.26 seconds
#> Iteration 30
#> Warning in se.gllvm(out): 9 parameter(s) have negative variance estimates (b x5, lambda x3, sigmaLV
#> x1). The model likely has not converged - consider re-fitting.
#> Iteration 30 finished in 1.24 seconds
#> Iteration 31
#> Warning in se.gllvm(out): 12 parameter(s) have negative variance estimates (b x12). The model
#> likely has not converged - consider re-fitting.
#> Iteration 31 finished in 2.08 seconds
#> Iteration 32
#> Warning in se.gllvm(out): 11 parameter(s) have negative variance estimates (b x11). The model
#> likely has not converged - consider re-fitting.
#> Iteration 32 finished in 1.73 seconds
#> Iteration 33
#> Warning in se.gllvm(out): 7 parameter(s) have negative variance estimates (b x3, lambda x3, sigmaLV
#> x1). The model likely has not converged - consider re-fitting.
#> Iteration 33 finished in 1.29 seconds
#> Iteration 34
#> Warning in se.gllvm(out): 11 parameter(s) have negative variance estimates (b x11). The model
#> likely has not converged - consider re-fitting.
#> Iteration 34 finished in 2.17 seconds
#> Iteration 35
#> Iteration 35 finished in 1.7 seconds
#> Iteration 36
#> Warning in se.gllvm(out): 40 parameter(s) have negative variance estimates (b x17, lambda x22,
#> sigmaLV x1). The model likely has not converged - consider re-fitting.
#> Iteration 36 finished in 2.46 seconds
#> Iteration 37
#> Warning in se.gllvm(out): 9 parameter(s) have negative variance estimates (b x9). The model likely
#> has not converged - consider re-fitting.
#> Iteration 37 finished in 1.85 seconds
#> Iteration 38
#> Warning in se.gllvm(out): 6 parameter(s) have negative variance estimates (b x6). The model likely
#> has not converged - consider re-fitting.
#> Iteration 38 finished in 1.97 seconds
#> Iteration 39
#> Warning in se.gllvm(out): 5 parameter(s) have negative variance estimates (b x5). The model likely
#> has not converged - consider re-fitting.
#> Iteration 39 finished in 1.33 seconds
#> Iteration 40
#> Warning in se.gllvm(out): 6 parameter(s) have negative variance estimates (b x6). The model likely
#> has not converged - consider re-fitting.
#> Iteration 40 finished in 1.6 seconds
#> Iteration 41
#> Iteration 41 finished in 1.65 seconds
#> Iteration 42
#> Warning in se.gllvm(out): 13 parameter(s) have negative variance estimates (b x13). The model
#> likely has not converged - consider re-fitting.
#> Iteration 42 finished in 1.8 seconds
#> Iteration 43
#> Iteration 43 finished in 1.76 seconds
#> Iteration 44
#> Iteration 44 finished in 1.62 seconds
#> Iteration 45
#> Warning in se.gllvm(out): 9 parameter(s) have negative variance estimates (b x9). The model likely
#> has not converged - consider re-fitting.
#> Iteration 45 finished in 1.31 seconds
#> Iteration 46
#> Warning in se.gllvm(out): 5 parameter(s) have negative variance estimates (lambda x4, sigmaLV x1).
#> The model likely has not converged - consider re-fitting.
#> Iteration 46 finished in 1.32 seconds
#> Iteration 47
#> Warning in se.gllvm(out): 8 parameter(s) have negative variance estimates (b x8). The model likely
#> has not converged - consider re-fitting.
#> Iteration 47 finished in 1.36 seconds
#> Iteration 48
#> Iteration 48 finished in 1.73 seconds
#> Iteration 49
#> Warning in se.gllvm(out): 23 parameter(s) have negative variance estimates (b x8, lambda x14,
#> sigmaLV x1). The model likely has not converged - consider re-fitting.
#> Iteration 49 finished in 1.29 seconds
#> Iteration 50
#> Warning in se.gllvm(out): 7 parameter(s) have negative variance estimates (b x7). The model likely
#> has not converged - consider re-fitting.
#> Iteration 50 finished in 1.28 seconds
```
# 4. Visualization

## 4.1 Occupancy Probability Caterpillar Plot


``` r

library(dplyr)
library(ggplot2)

otu_col <- "OTU"

# ------------------------------------------------------------
# 1. Combine retained GLLVM occupancy predictions
# ------------------------------------------------------------

psi_draws <- dplyr::bind_rows(out$psi_list)

# eta is on the logit scale
psi_draws <- psi_draws %>%
    dplyr::mutate(
        psi = stats::plogis(.data$eta)
    )

# ------------------------------------------------------------
# 2. Summarise occupancy at the OTU level
#
# This averages across sites and retained iterations.
# ------------------------------------------------------------

psi_otu_summary <- psi_draws %>%
    dplyr::group_by(.data[[otu_col]]) %>%
    dplyr::summarise(
        psi_mean = mean(.data$psi, na.rm = TRUE),
        psi_median = stats::median(.data$psi, na.rm = TRUE),
        psi_lwr = stats::quantile(
            .data$psi,
            0.025,
            na.rm = TRUE,
            names = FALSE
        ),
        psi_upr = stats::quantile(
            .data$psi,
            0.975,
            na.rm = TRUE,
            names = FALSE
        ),
        .groups = "drop"
    ) %>%
    dplyr::arrange(.data$psi_mean)

# ------------------------------------------------------------
# 3. Order OTUs by occupancy
# ------------------------------------------------------------

psi_otu_summary[[otu_col]] <- factor(
    psi_otu_summary[[otu_col]],
    levels = psi_otu_summary[[otu_col]]
)

# ------------------------------------------------------------
# 4. Caterpillar plot
# ------------------------------------------------------------

ggplot(
    psi_otu_summary,
    aes(
        x = psi_mean,
        y = .data[[otu_col]]
    )
) +
    geom_point(
        color = "darkgreen",
        size = 2
    ) +
    geom_errorbarh(
        aes(
            xmin = psi_lwr,
            xmax = psi_upr
        ),
        height = 0.2,
        color = "darkgreen"
    ) +
    labs(
        title = "Caterpillar Plot: Occupancy Probability by OTU",
        x = expression("Occupancy probability (" * psi * ")"),
        y = "OTU"
    ) +
    theme_minimal() +
    theme(
        axis.text.y = element_text(size = 8),
        plot.title = element_text(hjust = 0.5)
    )
#> Warning: `geom_errorbarh()` was deprecated in ggplot2 4.0.0.
#> ℹ Please use the `orientation` argument of `geom_errorbar()` instead.
#> This warning is displayed once per session.
#> Call `lifecycle::last_lifecycle_warnings()` to see where this warning was generated.
#> `height` was translated to `width`.
```

<img src="eDNAModel_gllvm_files/figure-html/plot-occupancy-1.png" alt="" width="960" />
## 4.2 Detection Probability Caterpillar Plot


``` r
library(dplyr)
library(ggplot2)

otu_col <- "OTU"

# ------------------------------------------------------------
# 1. Combine retained detection predictors
# ------------------------------------------------------------

p_detect_draws <- dplyr::bind_rows(out$p_detect_list)

# eta is on the cloglog scale
p_detect_draws <- p_detect_draws %>%
    dplyr::mutate(
        p_detect = 1 - exp(-exp(.data$eta))
    )

# ------------------------------------------------------------
# 2. Summarise detection probability by OTU
# ------------------------------------------------------------

p_detect_otu_summary <- p_detect_draws %>%
    dplyr::group_by(.data[[otu_col]]) %>%
    dplyr::summarise(
        p_detect_mean = mean(.data$p_detect, na.rm = TRUE),
        p_detect_median = stats::median(.data$p_detect, na.rm = TRUE),
        p_detect_lwr = stats::quantile(
            .data$p_detect,
            0.025,
            na.rm = TRUE,
            names = FALSE
        ),
        p_detect_upr = stats::quantile(
            .data$p_detect,
            0.975,
            na.rm = TRUE,
            names = FALSE
        ),
        .groups = "drop"
    ) %>%
    dplyr::arrange(.data$p_detect_mean)

# ------------------------------------------------------------
# 3. Order OTUs
# ------------------------------------------------------------

p_detect_otu_summary[[otu_col]] <- factor(
    p_detect_otu_summary[[otu_col]],
    levels = p_detect_otu_summary[[otu_col]]
)

# ------------------------------------------------------------
# 4. Caterpillar plot
# ------------------------------------------------------------

ggplot(
    p_detect_otu_summary,
    aes(
        x = p_detect_mean,
        y = .data[[otu_col]]
    )
) +
    geom_point(
        color = "steelblue",
        size = 2
    ) +
    geom_errorbarh(
        aes(
            xmin = p_detect_lwr,
            xmax = p_detect_upr
        ),
        height = 0.2,
        color = "steelblue"
    ) +
    labs(
        title = "Caterpillar Plot: Detection Probability by OTU",
        x = "Detection Probability",
        y = "OTU"
    ) +
    theme_minimal() +
    theme(
        axis.text.y = element_text(size = 8),
        plot.title = element_text(hjust = 0.5)
    )
#> `height` was translated to `width`.
```

<img src="eDNAModel_gllvm_files/figure-html/plot-detection-1.png" alt="" width="960" />

## 4.2  biplot


