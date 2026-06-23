##################################################
# 03_survival_and_validation_functions.R
#
# Survival analyses and Bulk validation analyses
##################################################
source("00_helper.R")

#' Fit Cox Proportional Hazards Model
#'
#' This function fits a Cox proportional hazards model using `survival::coxph`
#' and returns both the fitted model and a full model summary for reporting.
#'
#' @param data A data.frame containing survival and covariate information.
#' @param formula A survival formula of class `Surv(time, status) ~ predictors`.
#' @param return_model Logical. If TRUE, returns the fitted Cox model object.
#' @param return_summary Logical. If TRUE, returns model summary.
#'
#' @return A list containing:
#'   - model: fitted Cox model (optional)
#'   - summary: Cox model summary (coefficients, HR, p-values, etc.)
#'
#' @details
#' Uses `survival::coxph()` to fit a proportional hazards model.
#' Designed for reproducible survival analysis pipelines.
#'
#' @import survival
#' @export
run_cox_model <- function(data,
                          formula,
                          return_model = TRUE,
                          return_summary = TRUE) {
  
  # -----------------------------
  # 1. Fit Cox model
  # -----------------------------
  fit <- survival::coxph(formula, data = data)
  
  # -----------------------------
  # 2. Create output container
  # -----------------------------
  out <- list()
  
  # -----------------------------
  # 3. Store model (optional)
  # -----------------------------
  if (return_model) {
    out$model <- fit
  }
  
  # -----------------------------
  # 4. Store summary (publication output)
  # -----------------------------
  if (return_summary) {
    out$summary <- summary(fit)
  }
  
  # -----------------------------
  # 5. Attach model formula for traceability
  # -----------------------------
  out$formula <- formula
  
  return(out)
}


############################################################

#' Kaplan–Meier Survival Curve Estimation
#'
#' This function fits a Kaplan–Meier survival curve using `survival::survfit`.
#' It is designed for reproducible survival analysis and downstream plotting.
#'
#' @param data A data.frame containing survival information.
#' @param formula A survival formula of class `Surv(time, status) ~ group`.
#' @param conf.int Logical. If TRUE, includes confidence intervals.
#' @param return_plot Logical. If TRUE, returns a basic KM plot.
#'
#' @return A list containing:
#'   - fit: Kaplan–Meier survival fit object
#'   - plot: ggplot object (optional, if return_plot = TRUE)
#'
#' @details
#' Uses `survival::survfit()` to estimate non-parametric survival curves.
#' Suitable for stratified survival comparisons (e.g., groups, risk strata).
#'
#' @import survival
#' @export
run_kaplan_meier <- function(data,
                             formula,
                             conf.int = TRUE,
                             return_plot = FALSE) {
  
  # -----------------------------
  # 1. Fit Kaplan–Meier model
  # -----------------------------
  fit <- survival::survfit(
    formula,
    data = data,
    conf.int = conf.int
  )
  
  # -----------------------------
  # 2. Prepare output
  # -----------------------------
  out <- list()
  out$fit <- fit
  out$formula <- formula
  
  # -----------------------------
  # 3. Optional plot
  # -----------------------------
  if (return_plot) {
    out$plot <- survminer::ggsurvplot(
      fit,
      data = data,
      risk.table = TRUE,
      conf.int = conf.int
    )
  }
  
  return(out)
}



#' Adjust P-values for Multiple Testing
#'
#' This function performs multiple testing correction on a vector of p-values
#' using standard methods implemented in `stats::p.adjust`.
#'
#' @param pvalues A numeric vector of p-values.
#' @param method Character. Correction method for multiple testing adjustment.
#'   Options include: "bonferroni", "holm", "hochberg", "hommel",
#'   "BH" (Benjamini-Hochberg), "BY", "fdr", "none".
#'
#' @return A named numeric vector of adjusted p-values.
#'
#' @details
#' This function is a lightweight wrapper around `p.adjust()` to ensure
#' reproducibility and consistent reporting in statistical pipelines.
#'
#' @examples
#' adjust_pvalues(c(0.01, 0.04, 0.03), method = "BH")
#'
#' @export
adjust_pvalues <- function(pvalues,
                           method = "bonferroni") {
  
  # -----------------------------
  # Input validation
  # -----------------------------
  if (!is.numeric(pvalues)) {
    stop("pvalues must be a numeric vector.")
  }
  
  valid_methods <- c("bonferroni", "holm", "hochberg", "hommel",
                     "BH", "BY", "fdr", "none")
  
  if (!method %in% valid_methods) {
    stop(paste0("Invalid method. Choose from: ",
                paste(valid_methods, collapse = ", ")))
  }
  
  # -----------------------------
  # Adjustment
  # -----------------------------
  adjusted <- stats::p.adjust(pvalues, method = method)
  
  return(adjusted)
}


############################################################

#' Combine P-values Using Fisher's Method
#'
#' This function combines multiple p-values into a single global p-value
#' using Fisher's method (sum of log-transformed p-values).
#'
#' @param pvalues A numeric vector of p-values (must be between 0 and 1).
#'
#' @return A list containing:
#'   - statistic: Fisher's combined chi-square statistic
#'   - df: degrees of freedom
#'   - p.value: combined p-value
#'
#' @details
#' Fisher's method combines independent p-values using:
#' -2 * sum(log(p_i)) ~ Chi-square with 2k degrees of freedom
#'
#' This method assumes independence between tests.
#'
#' @examples
#' combine_pvalues(c(0.01, 0.04, 0.03))
#'
#' @export
combine_pvalues <- function(pvalues) {
  
  # -----------------------------
  # Input cleaning
  # -----------------------------
  pvalues <- as.numeric(pvalues)
  pvalues <- pvalues[!is.na(pvalues)]
  
  # -----------------------------
  # Validation
  # -----------------------------
  if (length(pvalues) == 0) {
    stop("No valid p-values provided.")
  }
  
  if (any(pvalues <= 0 | pvalues > 1)) {
    stop("All p-values must be in the interval (0, 1].")
  }
  
  # -----------------------------
  # Fisher's statistic
  # -----------------------------
  statistic <- -2 * sum(log(pvalues))
  
  df <- 2 * length(pvalues)
  
  p_combined <- stats::pchisq(statistic, df, lower.tail = FALSE)
  
  # -----------------------------
  # Return structured output
  # -----------------------------
  return(list(
    statistic = statistic,
    df = df,
    p.value = p_combined
  ))
}


############################################################

#' Decision Curve Analysis for Survival Outcomes using Cox Model
#'
#' This function performs Decision Curve Analysis (DCA) for survival data
#' using a Cox proportional hazards model and the `dcurves` package.
#' It evaluates clinical net benefit of a prognostic model at multiple time points.
#'
#' @param data A data.frame containing survival and predictor variables.
#' @param time_var Character. Name of survival time variable (e.g., days).
#' @param status_var Character. Name of event indicator (1 = event, 0 = censored).
#' @param predictor Character. Name of predictor variable in `data`.
#' @param times Numeric vector. Time points for DCA evaluation (same unit as time_var).
#' @param plot Logical. If TRUE, prints DCA plots.
#'
#' @return A list containing:
#'   - cox_model: fitted Cox proportional hazards model
#'   - data: dataset with added predicted risks at each time point
#'   - dca: named list of DCA objects for each time point
#'
#' @details
#' The function:
#' 1. Fits a Cox proportional hazards model
#' 2. Estimates linear predictor risk scores
#' 3. Converts risk to time-specific absolute risk
#' 4. Runs Decision Curve Analysis using `dcurves::dca()`
#'
#' @import survival dcurves
#' @export
run_dca_survival <- function(data,
                             time_var,
                             status_var,
                             predictor,
                             times = c(365, 1095, 1825),
                             plot = TRUE) {
  
  # Ensure required libraries
  # require(survival)
  # require(dcurves)
  
  # -----------------------------
  # 1. Build survival formula
  # -----------------------------
  formula_obj <- as.formula(
    paste0("Surv(", time_var, ", ", status_var, ") ~ ", predictor)
  )
  
  # -----------------------------
  # 2. Fit Cox model
  # -----------------------------
  cox_model <- survival::coxph(formula_obj, data = data)
  
  # Linear predictor (risk score)
  lp <- predict(cox_model, type = "lp")
  
  # Store results
  dca_results <- list()
  data_out <- data
  
  # -----------------------------
  # 3. Loop over time points
  # -----------------------------
  for (t in times) {
    
    # Survival curve at time t
    base_surv <- summary(survival::survfit(cox_model), times = t)$surv
    
    # Avoid NA issues
    if (is.na(base_surv)) next
    
    # Convert to absolute risk
    risk_col <- paste0("risk_", t)
    data_out[[risk_col]] <- 1 - (base_surv ^ exp(lp))
    
    # DCA formula
    dca_formula <- as.formula(
      paste0("Surv(", time_var, ", ", status_var, ") ~ ", risk_col)
    )
    
    # Run DCA
    dca_model <- dcurves::dca(
      formula = dca_formula,
      data = data_out,
      time = t
    )
    
    dca_results[[paste0("DCA_", t)]] <- dca_model
    
    # Optional plot
    if (plot) {
      print(
        plot(dca_model) +
          ggtitle(paste0("Decision Curve Analysis (", t, " days)"))
      )
    }
  }
  
  # -----------------------------
  # 4. Return structured output
  # -----------------------------
  return(list(
    cox_model = cox_model,
    data = data_out,
    dca = dca_results
  ))
}


############################################################

#' Extract Metadata from Filename
#'
#' This function extracts structured information (e.g., patient ID and cluster ID)
#' from standardized filenames using regular expression parsing.
#'
#' Expected filename format:
#'   patientID.CXvCY (e.g., Sample123.C1vC2.txt)
#'
#' @param filename A character string or vector of file paths.
#'
#' @return A data.frame with:
#'   - patient_id: extracted patient/sample identifier
#'   - cluster_id: extracted cluster identifier
#'
#' @details
#' Uses regular expression matching to parse structured filenames.
#' Automatically strips file paths using `basename()`.
#'
#' @examples
#' extract_filename_info("Sample123.C1vC2.txt")
#'
#' @export
extract_filename_info <- function(filename) {
  
  # -----------------------------
  # Load dependency safely
  # -----------------------------
  require(stringr)
  
  # -----------------------------
  # Clean filename path
  # -----------------------------
  filename <- basename(filename)
  
  # -----------------------------
  # Regex extraction
  # -----------------------------
  info <- stringr::str_match(
    filename,
    "^(\\w+)\\.(C\\d+)vC\\d+"
  )
  
  # -----------------------------
  # Validation
  # -----------------------------
  if (all(is.na(info[, 2]))) {
    stop("Filename format not recognized. Expected format: patientID.CXvCY")
  }
  
  # -----------------------------
  # Build output
  # -----------------------------
  out <- data.frame(
    patient_id = info[, 2],
    cluster_id = info[, 3],
    stringsAsFactors = FALSE
  )
  
  return(out)
}
