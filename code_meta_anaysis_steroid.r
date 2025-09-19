# Required libraries
library(patchwork)
library(MetaStan)
library(openxlsx)
library(tidyverse)
library(grid)
library(meta)
library(metafor)
library(ggplot2)
library(gridExtra)
library(bridgesampling)
library(effectsize)
library(dplyr)
library(flextable)
library(officer)
library(scales)
library(stringr)
library(knitr)
library(data.table)    


# Load dataset
load("meta_analysis_steroids_ACP.RData")



######################################################################################################################################################################################################################
######################################################################################################################################################################################################################
######################################################################################################################################################################################################################
######################################################################################################################################################################################################################
######################################################################################################################################################################################################################
######################################################################################################################################################################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
####                    FUNCTIONS                 ####
###############################################################################
###############################################################################
###############################################################################
###############################################################################
######################################################################################################################################################################################################################
######################################################################################################################################################################################################################
######################################################################################################################################################################################################################
######################################################################################################################################################################################################################
######################################################################################################################################################################################################################
######################################################################################################################################################################################################################

###############################################################################
###############################################################################
# FREQUENTIST META-ANALYSIS FUNCTIONS
###############################################################################
###############################################################################

# Main function for meta-analysis of clinical trials
# Function handles both binary outcomes
analysis_meta <- function(
  outcomes = "mortality", 
  data,
  metabin = FALSE,
  use_hartung_knapp = FALSE,
  cont = FALSE,
  event.e = NULL,
  event.c = NULL,
  n.e = NULL,
  n.c = NULL,
  mean.e = NULL,
  sd.e = NULL,
  mean.c = NULL,
  sd.c = NULL,
  studlab = NULL,
  subgroup = NULL,
  xlab_custom = "Risk Ratio (95% CI)",
  include_fixed = FALSE,
  subgroup_remove = NULL,
  exclusion_condition = NULL,
  subgroup_levels = NULL,
  iqwig = FALSE,
  use_reml = FALSE,  
  sm = "RR",
  favors_steroids_right = TRUE,
  extra_columns = NULL,
  rob_column = NULL  #
) {
  
  # Set general meta-analysis parameters 
  # IQWIG METHOD: Uses Paule-Mandel with Hartung-Knapp adjustment
  if (iqwig) {
    settings.meta(
      method.random.ci = "HK",
      adhoc.hakn.ci = "IQWiG6",
      method.tau = "PM",
      tau.common = FALSE,
      MH.exact = FALSE,
      RR.Cochrane = TRUE,
      Q.Cochrane = TRUE,
      prediction = TRUE,
      test.overall = TRUE,
      test.subgroup = !is.null(subgroup),
      test.effect.subgroup = !is.null(subgroup),
      digits.I2 = 0,
      digits.tau2 = 3,
      digits.tau = 4,
      CIbracket = "[",
      CIseparator = ", ",
      header.line = TRUE
    )
  } else {
    # STANDARD SETTINGS: DerSimonian-Laird as default
    settings.meta(
      method.random.ci = if (use_hartung_knapp) "HK" else "classic",
      method.tau = if (use_reml) "REML" else "DL",  # USE REML IF SPECIFIED
      tau.common = FALSE,
      MH.exact = FALSE,
      RR.Cochrane = TRUE,
      Q.Cochrane = TRUE,
      prediction = FALSE,
      test.overall = TRUE,
      test.subgroup = !is.null(subgroup),
      test.effect.subgroup = !is.null(subgroup),
      digits.I2 = 0,
      digits.tau2 = 3,
      digits.tau = 4,
      CIbracket = "[",
      CIseparator = ", ",
      header.line = TRUE
    )
  }

  # Filter outcome of interest
  filter <- data[data$outcomes == outcomes, ]
  
  # Remove specific subgroup if requested
  if (!is.null(subgroup) && !is.null(subgroup_remove)) {
    filter <- filter[filter[[subgroup]] != subgroup_remove, ]
  }
  
  # Apply exclusion criteria if provided
  if (!is.null(exclusion_condition) && exclusion_condition != "") {
    filter <- filter[!(eval(parse(text = exclusion_condition), envir = filter)), ]
  }
  
  # Set subgroup factor levels if specified
  if (!is.null(subgroup) && !is.null(subgroup_levels)) {
    filter[[subgroup]] <- factor(filter[[subgroup]], levels = subgroup_levels)
  }

  # Perform meta-analysis based on data type
  if (cont) {
    # FOR CONTINUOUS OUTCOMES: means and standard deviations
    meta <- metacont(
      n.e = filter[[n.e]],
      mean.e = filter[[mean.e]],
      sd.e = filter[[sd.e]],
      n.c = filter[[n.c]],
      mean.c = filter[[mean.c]],
      sd.c = filter[[sd.c]],
      studlab = filter[[studlab]],
      data = filter,
      sm = "SMD",  # Standardized mean difference
      method.smd = "Cohen",
      comb.random = TRUE,
      comb.fixed = include_fixed,
      method.tau = if (use_reml) "REML" else "DL",  # REML FOR CONTINUOUS DATA
      subgroup = if (!is.null(subgroup)) filter[[subgroup]] else NULL
    )
  } else {
    # FOR BINARY OUTCOMES: events and totals
    meta <- metabin(
      event.e = filter[[event.e]],
      n.e = filter[[n.e]],
      event.c = filter[[event.c]],
      n.c = filter[[n.c]],
      studlab = filter[[studlab]],
      data = filter,
      sm = sm,  
      method = "MH",  # Mantel-Haenszel method
      method.tau = if (iqwig) "PM" else if (use_reml) "REML" else "DL",  # REML FOR BINARY DATA
      common = include_fixed,
      random = TRUE,
      method.random.ci = if (iqwig) "HK" else if (use_hartung_knapp) "HK" else "classic",
      subgroup = if (!is.null(subgroup)) filter[[subgroup]] else NULL
    )
  }

  # Set group labels for the analysis
  if (!cont) {
    meta$label.e <- "Corticosteroids"
    meta$label.c <- "Standard of care"
  }
  meta$label.e <- "Steroids"
  meta$label.c <- "Control"

  # Clean up subgroup names by removing prefix
  if (!is.null(meta$byvar)) {
    levels(meta$byvar) <- gsub("^.*=", "", levels(meta$byvar))
  }

  # Set direction labels for forest plot
  left_label <- if (favors_steroids_right) "Favors corticosteroids" else "Favors control"
  right_label <- if (favors_steroids_right) "Favors control" else "Favors corticosteroids"

  outcome_title <- paste0("**", outcomes, "**")

  # Basic column setup for standard forest plot
  leftcols <- c("studlab", "event.e", "n.e", "event.c", "n.c", "effect.ci", "w.random")
  leftlabs <- c("Study", "Steroids (events/N)", "Control (events/N)", "RR (95% CI)", "Weight (%)")

  # Check if we need enhanced format with extra columns or risk of bias
  if (!is.null(extra_columns) || !is.null(rob_column)) {
    # Add extra columns to meta object
    if (!is.null(extra_columns)) {
      for (col in extra_columns) {
        if (col %in% names(filter)) {
          meta[[col]] <- filter[[col]]
        }
      }
    }
    
    # ADD RISK OF BIAS COLUMN IF SPECIFIED
    if (!is.null(rob_column) && rob_column %in% names(filter)) {
      meta$risk_of_bias <- filter[[rob_column]]
    }
    
    # CREATE COMBINED COLUMNS IN META OBJECT FOR CLEANER PRESENTATION
    meta$steroids_combined <- paste(meta$event.e, meta$n.e, sep="/")
    meta$control_combined <- paste(meta$event.c, meta$n.c, sep="/")
    
    # NEW SIMPLIFIED COLUMNS AND LABELS FOR LEFT SIDE
    leftcols_new <- c("studlab", "steroids_combined", "control_combined", "effect.ci", "w.random")
    leftlabs_new <- c("Study", "Steroids", "Control", "RR (95% CI)", "Weight (%)")
    
    # BUILD RIGHT SIDE COLUMNS FOR ADDITIONAL INFORMATION
    rightcols_final <- extra_columns
    rightlabs_final <- extra_columns
    
    # ADD RISK OF BIAS TO RIGHT SIDE IF SPECIFIED
    if (!is.null(rob_column) && rob_column %in% names(filter)) {
      rightcols_final <- c(rightcols_final, "risk_of_bias")
      rightlabs_final <- c(rightlabs_final, "Risk of Bias")
    }
    
    # Create enhanced forest plot with additional columns on the right
    forest(
      meta,
      layout = "meta",
      comb.random = TRUE,
      comb.fixed = include_fixed,
      fontsize = 9,
      just.studlab = "left",
      lwd = 1,
      lwd.square = 1,
      lwd.diamond = 1,
      col.square = "black",           
      col.square.lines = "black",     
      col.diamond = "#2E8B7A",        
      col.diamond.lines = "#2E8B7A", 
      col.predict = "black",
      col.predict.lines = "black",
      col.lines = "black",
      col.label = "black",
      fill = "transparent",
      type.study = "square",
      type.common = "diamond",
      type.random = "diamond",
      text.addline1 = outcome_title,
      ff.addline1 = "bold",
      colgap.forest.left = "12mm",
      colgap.forest.right = "8mm",
      leftlabs = leftlabs_new,        
      leftcols = leftcols_new,       
      rightcols = rightcols_final,    
      rightlabs = rightlabs_final,   
      addrow.subgroups = TRUE,        
      addrows.below.overall = 2,      
      test.effect.subgroup = FALSE,   
      test.subgroup = TRUE,           
      print.subgroup.name = FALSE,    
      sep.subgroup = "",              
      print.I2 = TRUE,
      print.tau2 = TRUE,
      xlab = xlab_custom,
      label.left = left_label,
      label.right = right_label,
      print.byvar = TRUE
    )
  } else {
    # Standard forest plot with JAMA style layout (original format)
    forest(
      meta,
      layout = "JAMA",
      comb.random = TRUE,
      comb.fixed = include_fixed,
      fontsize = 9,
      just.studlab = "left",
      lwd = 1,
      lwd.square = 1,
      lwd.diamond = 1,
      col.square = "black",
      col.square.lines = "black",
      col.diamond = "black",
      col.diamond.lines = "black",
      col.predict = "black",
      col.predict.lines = "black",
      col.lines = "black",
      col.label = "black",
      fill = "transparent",
      type.study = "square",
      type.common = "diamond",
      type.random = "diamond",
      text.addline1 = outcome_title,
      ff.addline1 = "bold",
      colgap.forest.left = "12mm",
      colgap.forest.right = "8mm",
      leftlabs = leftlabs,
      leftcols = leftcols,
      leftlabs.format = list(
        event.e = function(x, i) paste(meta$event.e[i], "/", meta$n.e[i]),
        event.c = function(x, i) paste(meta$event.c[i], "/", meta$n.c[i])
      ),
      print.I2 = TRUE,
      print.tau2 = TRUE,
      xlab = xlab_custom,
      label.left = left_label,
      label.right = right_label,
      print.byvar = TRUE
    )
  }

  # Capture the plot for return
  plot_grob <- recordPlot()
  return(list(meta = meta, plot = plot_grob))



  attr(meta, "plot") <- plot_grob
  return(meta)
}


analysis_meta3 <- function(
  outcomes = "mortality", 
  data,
  metabin = FALSE,
  use_hartung_knapp = FALSE,
  cont = FALSE,
  event.e = NULL,
  event.c = NULL,
  n.e = NULL,
  n.c = NULL,
  mean.e = NULL,
  sd.e = NULL,
  mean.c = NULL,
  sd.c = NULL,
  studlab = NULL,
  subgroup = NULL,
  xlab_custom = "Risk Ratio (95% CI)",
  include_fixed = FALSE,
  subgroup_remove = NULL,
  exclusion_condition = NULL,
  subgroup_levels = NULL,
  iqwig = FALSE,
  sm = "RR",
  favors_steroids_right = TRUE,
  extra_columns = NULL,
  rob_column = NULL  
) {
  
  # Same general parameter setup as analysis_meta
  if (iqwig) {
    settings.meta(
      method.random.ci = "HK",
      adhoc.hakn.ci = "IQWiG6",
      method.tau = "PM",
      tau.common = FALSE,
      MH.exact = FALSE,
      RR.Cochrane = TRUE,
      Q.Cochrane = TRUE,
      prediction = TRUE,
      test.overall = TRUE,
      test.subgroup = !is.null(subgroup),
      test.effect.subgroup = !is.null(subgroup),
      digits.I2 = 0,
      digits.tau2 = 3,
      digits.tau = 4,
      CIbracket = "[",
      CIseparator = ", ",
      header.line = TRUE
    )
  } else {
    settings.meta(
      method.random.ci = if (use_hartung_knapp) "HK" else "classic",
      method.tau = "DL",
      tau.common = FALSE,
      MH.exact = FALSE,
      RR.Cochrane = TRUE,
      Q.Cochrane = TRUE,
      prediction = FALSE,
      test.overall = TRUE,
      test.subgroup = !is.null(subgroup),
      test.effect.subgroup = !is.null(subgroup),
      digits.I2 = 0,
      digits.tau2 = 3,
      digits.tau = 4,
      CIbracket = "[",
      CIseparator = ", ",
      header.line = TRUE
    )
  }

  # Same filtering approach as the original function
  filter <- data[data$outcomes == outcomes, ]
  if (!is.null(subgroup) && !is.null(subgroup_remove)) {
    filter <- filter[filter[[subgroup]] != subgroup_remove, ]
  }
  if (!is.null(exclusion_condition) && exclusion_condition != "") {
    filter <- filter[!(eval(parse(text = exclusion_condition), envir = filter)), ]
  }
  if (!is.null(subgroup) && !is.null(subgroup_levels)) {
    filter[[subgroup]] <- factor(filter[[subgroup]], levels = subgroup_levels)
  }


  if (cont) {
    meta <- metacont(
      n.e = filter[[n.e]],
      mean.e = filter[[mean.e]],
      sd.e = filter[[sd.e]],
      n.c = filter[[n.c]],
      mean.c = filter[[mean.c]],
      sd.c = filter[[sd.c]],
      studlab = filter[[studlab]],
      data = filter,
      sm = "SMD",
      method.smd = "Cohen",
      random = TRUE,
      common = include_fixed,
      subgroup = if (!is.null(subgroup)) filter[[subgroup]] else NULL
    )
  } else {
    meta <- metabin(
      event.e = filter[[event.e]],
      n.e = filter[[n.e]],
      event.c = filter[[event.c]],
      n.c = filter[[n.c]],
      studlab = filter[[studlab]],
      data = filter,
      sm = sm,
      method = "MH",
      method.tau = if (iqwig) "PM" else "DL",
      common = include_fixed,
      random = TRUE,
      method.random.ci = if (iqwig) "HK" else if (use_hartung_knapp) "HK" else "classic",
      subgroup = if (!is.null(subgroup)) filter[[subgroup]] else NULL
    )
  }

  # Same labeling as original function
  if (!cont) {
    meta$label.e <- "Corticosteroids"
    meta$label.c <- "Standard of care"
  }
  meta$label.e <- "Steroids"
  meta$label.c <- "Control"

  if (!is.null(meta$byvar)) {
    levels(meta$byvar) <- gsub("^.*=", "", levels(meta$byvar))
  }

  # Same direction setup
  left_label <- if (favors_steroids_right) "Favors corticosteroids" else "Favors control"
  right_label <- if (favors_steroids_right) "Favors control" else "Favors corticosteroids"
  outcome_title <- paste0("**", outcomes, "**")

  # Basic columns setup
  leftcols <- c("studlab", "event.e", "n.e", "event.c", "n.c", "effect.ci", "w.random")
  leftlabs <- c("Study", "Steroids (events/N)", "Control (events/N)", "RR (95% CI)", "Weight (%)")

  # Handle extra columns or risk of bias column
  if (!is.null(extra_columns) || !is.null(rob_column)) {
    # Add extra columns to meta object
    if (!is.null(extra_columns)) {
      for (col in extra_columns) {
        if (col %in% names(filter)) {
          meta[[col]] <- filter[[col]]
        }
      }
    }
    
    # ADD RISK OF BIAS COLUMN IF SPECIFIED
    if (!is.null(rob_column) && rob_column %in% names(filter)) {
      meta$risk_of_bias <- filter[[rob_column]]
    }
    
    # CREATE COMBINED COLUMNS 
    if (cont) {
      # For continuous outcomes - show mean (SD)
      meta$steroids_combined <- paste(formatC(meta$mean.e, format = "f", digits = 1), 
                                     formatC(meta$sd.e, format = "f", digits = 1), sep=" (") 
      meta$steroids_combined <- paste0(meta$steroids_combined, ")")
      meta$control_combined <- paste(formatC(meta$mean.c, format = "f", digits = 1), 
                                    formatC(meta$sd.c, format = "f", digits = 1), sep=" (")
      meta$control_combined <- paste0(meta$control_combined, ")")
    } else {
      # For binary outcomes 
      meta$steroids_combined <- paste(meta$event.e, meta$n.e, sep="/")
      meta$control_combined <- paste(meta$event.c, meta$n.c, sep="/")
    }
    
    # NEW COLUMNS AND LABELS with adaptation for continuous data
    leftcols_new <- c("studlab", "steroids_combined", "control_combined", "effect.ci", "w.random")
    if (cont) {
      leftlabs_new <- c("Study", "Steroids", "Control", "SMD (95% CI)", "Weight (%)")
      xlab_custom <- "SMD (95% CI)"
    } else {
      leftlabs_new <- c("Study", "Steroids", "Control", "RR (95% CI)", "Weight (%)")
    }
    
    # BUILD RIGHT COLUMNS
    rightcols_final <- extra_columns
    rightlabs_final <- extra_columns
    
    # ADD RISK OF BIAS IF SPECIFIED
    if (!is.null(rob_column) && rob_column %in% names(filter)) {
      rightcols_final <- c(rightcols_final, "risk_of_bias")
      rightlabs_final <- c(rightlabs_final, "Risk of Bias")
    }
    
    # Forest plot with right columns
    forest(
      meta,
      layout = "meta",
      random = TRUE,
      common = include_fixed,
      fontsize = 9,
      just.studlab = "left",
      lwd = 1,
      lwd.square = 1,
      lwd.diamond = 1,
      col.square = "black",           
      col.square.lines = "black",     
      col.diamond = "#2E8B7A",       
      col.diamond.lines = "#2E8B7A",  
      col.predict = "black",
      col.predict.lines = "black",
      col.lines = "black",
      col.label = "black",
      fill = "transparent",
      type.study = "square",
      type.common = "diamond",
      type.random = "diamond",
      text.addline1 = outcome_title,
      ff.addline1 = "bold",
      colgap.forest.left = "12mm",
      colgap.forest.right = "8mm",
      leftlabs = leftlabs_new,        
      leftcols = leftcols_new,        
      rightcols = rightcols_final,    
      rightlabs = rightlabs_final,   
      addrow.subgroups = TRUE,        
      addrows.below.overall = 2,      
      test.effect.subgroup = FALSE,   
      test.subgroup = TRUE,           
      print.subgroup.name = FALSE,    
      sep.subgroup = "",              
      print.I2 = TRUE,
      print.tau2 = TRUE,
      xlab = xlab_custom,
      label.left = left_label,
      label.right = right_label,
      print.byvar = TRUE
    )
  } else {
    # Standard forest plot with JAMA format
    forest(
      meta,
      layout = "JAMA",
      random = TRUE,
      common = include_fixed,
      fontsize = 9,
      just.studlab = "left",
      lwd = 1,
      lwd.square = 1,
      lwd.diamond = 1,
      col.square = "black",
      col.square.lines = "black",
      col.diamond = "black",
      col.diamond.lines = "black",
      col.predict = "black",
      col.predict.lines = "black",
      col.lines = "black",
      col.label = "black",
      fill = "transparent",
      type.study = "square",
      type.common = "diamond",
      type.random = "diamond",
      text.addline1 = outcome_title,
      ff.addline1 = "bold",
      colgap.forest.left = "12mm",
      colgap.forest.right = "8mm",
      leftlabs = leftlabs,
      leftcols = leftcols,
      leftlabs.format = list(
        event.e = function(x, i) paste(meta$event.e[i], "/", meta$n.e[i]),
        event.c = function(x, i) paste(meta$event.c[i], "/", meta$n.c[i])
      ),
      print.I2 = TRUE,
      print.tau2 = TRUE,
      xlab = xlab_custom,
      label.left = left_label,
      label.right = right_label,
      print.byvar = TRUE
    )
  }

  # Return results
  plot_grob <- recordPlot()
  return(list(meta = meta, plot = plot_grob))

  # Alternative return
  plot_grob <- recordPlot()
  attr(meta, "plot") <- plot_grob
  return(meta)
}


###############################################################################
###############################################################################
# BAYESIAN
###############################################################################
###############################################################################

# Bayesian function
bayesian <- function(outcomes = "neuro",
                     data = dt,
                     event_placebo_col = "event_placebo",
                     total_placebo_col = "total_placebo",
                     event_steroid_col = "event_steroid",
                     total_steroid_col = "total_steroid",
                     study_label_col = NULL,
                     thetaprior = c(-0.5, 0.5),
                     returnstan = TRUE,
                     return_summary = FALSE,  
                     tauprior = 0.5,
                     adapt_delta = 0.95) {

  if (!requireNamespace("MetaStan", quietly = TRUE) || !requireNamespace("gridExtra", quietly = TRUE)) {
    stop("Packages MetaStan and gridExtra must be installed.")
  }

  filter <- data[data$outcomes == outcomes, ]

  required_cols <- c(event_placebo_col, total_placebo_col, event_steroid_col, total_steroid_col)
  if (!all(required_cols %in% colnames(filter))) {
    stop("One or more specified columns do not exist in the dataset.")
  }

  filter$r1 <- filter[[event_placebo_col]]
  filter$n1 <- filter[[total_placebo_col]]
  filter$r2 <- filter[[event_steroid_col]]
  filter$n2 <- filter[[total_steroid_col]]

  filter <- filter[complete.cases(filter[, c("r1", "n1", "r2", "n2")]), ]
  if (nrow(filter) == 0) {
    stop("No complete cases available after removing missing values.")
  }

  prep_data <- create_MetaStan_dat(dat = filter,
                                   armVars = c(responders = "r", sampleSize = "n"))

  set.seed(1234)
  stan_object <- meta_stan(data = prep_data,
                           likelihood = "binomial",
                           re = TRUE,
                           tau_prior = tauprior,
                           ncp = TRUE,
                           mu_prior = c(0, 10),
                           theta_prior = thetaprior,
                           tau_prior_dist = "half-normal",
                           adapt_delta = adapt_delta)

  if (return_summary) {
    post_theta <- rstan::extract(stan_object$fit, pars = "theta")$theta
    post_tau <- rstan::extract(stan_object$fit, pars = "tau")$tau

    summary_df <- data.frame(
      Mean_theta = mean(post_theta),
      CrI_2.5 = quantile(post_theta, 0.025),
      CrI_97.5 = quantile(post_theta, 0.975),
      Tau_median = median(post_tau),
      Prob_benefit = mean(post_theta < 0)
    )

    return(summary_df)
  }

  if (!is.null(study_label_col)) {
    if (!(study_label_col %in% colnames(filter))) {
      stop(paste0("Column '", study_label_col, "' not found for study labels."))
    }
    forest <- forest_plot(stan_object, xlab = "log-OR", labels = filter[[study_label_col]])
  } else {
    forest <- forest_plot(stan_object, xlab = "log-OR")
  }

  if (returnstan) {
    return(stan_object)
  } else {
    return(forest)
  }
}

# Plot function for prior and posterior distributions
plot_prior_posterior <- function(prior_mean, prior_sd, posterior_samples, 
                                 prior_label,  
                                 xlim = c(-1, 2.5), 
                                 ylim = NULL,  
                                 prior_color = "#B3B3B3", 
                                 posterior_color = "#BC3C29", 
                                 ci_color = "#E15759",  
                                 zero_line_color = "#2A475E") {  
    # Check parameters
    if (!all(is.finite(xlim))) stop("Invalid xlim: ", paste(xlim, collapse = ", "))
    if (!is.finite(prior_mean) || !is.finite(prior_sd) || prior_sd <= 0) {
        stop("Invalid prior parameters: mean = ", prior_mean, ", sd = ", prior_sd)
    }
    
    # Density of prior
    theta_prior <- seq(xlim[1], xlim[2], length.out = 1000)
    prior_density <- dnorm(theta_prior, mean = prior_mean, sd = prior_sd)
    prior_density_df <- data.frame(theta = theta_prior, density = prior_density)
    
    # Posterior density
    posterior_density <- density(posterior_samples, from = xlim[1], to = xlim[2])
    posterior_density_df <- data.frame(
        theta = posterior_density$x,
        density = posterior_density$y
    )
    
    # Calculate 95% credible intervals
    ci_posterior <- quantile(posterior_samples, probs = c(0.025, 0.975))
    ci_prior <- c(prior_mean - 1.96 * prior_sd, prior_mean + 1.96 * prior_sd)
    
    # Define Y-axis limits
    if (is.null(ylim)) {
        ylim <- c(0, max(c(prior_density_df$density, posterior_density_df$density)) * 1.05)
    }
    
    # Plot
    plot <- ggplot() +
        geom_line(data = prior_density_df, aes(x = theta, y = density), 
                  color = prior_color, size = 1.2, linetype = "dashed") +
        geom_ribbon(data = prior_density_df[prior_density_df$theta >= ci_prior[1] & prior_density_df$theta <= ci_prior[2], ],
                    aes(x = theta, ymin = 0, ymax = density), fill = prior_color, alpha = 0.3) +
        geom_line(data = posterior_density_df, aes(x = theta, y = density), 
                  color = posterior_color, size = 1.5) +
        geom_ribbon(data = posterior_density_df[posterior_density_df$theta >= ci_posterior[1] & posterior_density_df$theta <= ci_posterior[2], ],
                    aes(x = theta, ymin = 0, ymax = density), fill = posterior_color, alpha = 0.3) +
        geom_vline(xintercept = ci_posterior, linetype = "dotted", color = ci_color, size = 0.8) +
        geom_vline(xintercept = 0, linetype = "dashed", color = zero_line_color, size = 1.0) +
        labs(
            subtitle = prior_label,  
            x = expression(theta),
            y = NULL
        ) +
        scale_x_continuous(limits = xlim, breaks = seq(xlim[1], xlim[2], by = 0.5)) +
        scale_y_continuous(limits = ylim, expand = c(0, 0)) +
        theme_minimal(base_size = 14) +
        theme(
            plot.subtitle = element_text(hjust = 0.5, size = 8, face = "italic"),  
            axis.title.x = element_text(size = 12),
            axis.text.x = element_text(size = 10),
            axis.line.x = element_line(color = "black"),  
            axis.text.y = element_blank(),               
            axis.ticks.y = element_blank(),              
            axis.line.y = element_blank(),               
            panel.grid = element_blank(),                
            panel.background = element_rect(fill = "white", color = NA),  
            plot.background = element_rect(fill = "white", color = NA),   
            legend.position = "none"                     
        )
    
    return(plot)

}


# Function to create multiple prior/posterior plots
make_prior_plots <- function(liste, prior_labels, titre) {
  plots <- list()
  
  for (i in seq_along(liste)) {
    element <- liste[[i]]
    
    # Extract samples and prior parameters
    theta_samples <- rstan::extract(element$fit, pars = "theta", permuted = TRUE)$theta
    prior_mean <- element$stanDat$theta_prior[1]
    prior_sd <- sqrt(element$stanDat$theta_prior[2])
    
    # Check parameters
    if (!all(is.finite(theta_samples))) {
      message("Skipping element ", i, ": Invalid posterior samples")
      next
    }
    if (!is.finite(prior_mean) || !is.finite(prior_sd) || prior_sd <= 0) {
      message("Skipping element ", i, ": Invalid prior parameters (mean = ", prior_mean, 
              ", sd = ", prior_sd, ")")
      next
    }
    
    # Create plot
    plot <- plot_prior_posterior(prior_mean, prior_sd, theta_samples, prior_labels[i])
    
    plots[[i]] <- plot
  }
  
  # Combine with patchwork
  big_plot <- (plots[[1]] + plots[[2]]) /
              (plots[[3]] + plots[[4]]) /
              (plots[[5]] + plot_spacer()) +
    plot_layout(heights = c(1,1,1)) +
    plot_annotation(
      title = titre,
      theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
    )
  
  return(big_plot)
}

# Calcul I2 with Bayesian model
calc_I2_bayes <- function(fit, data, tau_param = "tau") {
  # Extraction des données
  a <- data$event_steroid
  b <- data$total_steroid - data$event_steroid  
  c <- data$event_placebo
  d <- data$total_placebo - data$event_placebo
  
  # Correction de continuité si nécessaire
  zero_cells <- (a == 0 | b == 0 | c == 0 | d == 0)
  if(any(zero_cells, na.rm=TRUE)) {
    message(paste("Correction de continuité appliquée à", sum(zero_cells, na.rm=TRUE), "études"))
    a[zero_cells] <- a[zero_cells] + 0.5
    b[zero_cells] <- b[zero_cells] + 0.5  
    c[zero_cells] <- c[zero_cells] + 0.5
    d[zero_cells] <- d[zero_cells] + 0.5
  }
  
  # Calcul  de la variance et I2
  vi <- 1/a + 1/b + 1/c + 1/d
  vbar <- mean(vi, na.rm=TRUE)
  
  tau <- rstan::extract(fit, tau_param)[[1]]
  I2 <- 100 * tau^2 / (tau^2 + vbar)
  
  # Stat
  q <- quantile(I2, c(0.5, 0.025, 0.975), na.rm = TRUE)
  cat(sprintf("I² = %.1f%% [95%% CrI: %.1f%% - %.1f%%]\n", 
              q[1], q[2], q[3]))
  
  return(invisible(list(
    I2_median = q[1],
    I2_CI = c(q[2], q[3]),
    I2_samples = I2,
    tau_median = median(tau),
    vbar = vbar
  )))
}




######################################################################################################################################################################################################################
######################################################################################################################################################################################################################
######################################################################################################################################################################################################################
######################################################################################################################################################################################################################
######################################################################################################################################################################################################################
######################################################################################################################################################################################################################
##############################################################################################################################################################
##############################################################################################################################################################
##############################################################################################################################################################
####                     INFECTION                                         ####
##############################################################################################################################################################
##############################################################################################################################################################
##############################################################################################################################################################
######################################################################################################################################################################################################################
######################################################################################################################################################################################################################
######################################################################################################################################################################################################################
######################################################################################################################################################################################################################
######################################################################################################################################################################################################################
######################################################################################################################################################################################################################

###############################################################################
###############################################################################
###############################################################################
####                    Hospital-acquired infections                       ####
###############################################################################
###############################################################################
###############################################################################


# Primary analysis (low dose, short duration) frequentist : Paule-Mandel + Hartung-Knapp
infection_ha <- analysis_meta(
  outcomes = "Hospital-acquired infections",
  data = dt[dt$outcomes == "Hospital-acquired infections" & !is.na(dt$event_placebo) & dt$CI==TRUE, ]
  [order(dt[dt$outcomes == "Hospital-acquired infections" & !is.na(dt$event_placebo) & dt$CI==TRUE, ]$steroid), ], # CI = low dose, short duration,timing < 7 days
  metabin = TRUE,
  use_hartung_knapp = T,
  iqwig = T,
  sm = "RR",
  event.e = "event_steroid",
  n.e = "total_steroid",
  event.c = "event_placebo",
  n.c = "total_placebo",
  studlab = "Authors+O44LA1:R41",
  subgroup = "type2",
  subgroup_levels = c("Non severe pneumonia", "Severe pneumonia", "ARDS"),
  extra_columns = c( "Oxygen", "ICU", "Steroid","Daily Dose (mg)", "Duration (days)"),      
  rob_column = "rob2_patterns"       # AJOUTE LA COLONNE ROB2
)


# Secondary analysis (without restriction on dose, duration or timing) frequentist : Paule-Mandel + Hartung-Knapp
infection_ha_all <- analysis_meta(
  outcomes = "Hospital-acquired infections",
  data = dt[dt$outcomes == "Hospital-acquired infections" & !is.na(dt$event_placebo), ]
 [order(dt[dt$outcomes == "Hospital-acquired infections" & !is.na(dt$event_placebo), ]$steroid), ],
  metabin = TRUE,
  use_hartung_knapp = T,
  iqwig = T,
  sm = "RR",
  event.e = "event_steroid",
  n.e = "total_steroid",
  event.c = "event_placebo",
  n.c = "total_placebo",
  studlab = "Authors+O44LA1:R41",
  subgroup = "type2",
  subgroup_levels = c("Non severe pneumonia", "Severe pneumonia", "ARDS"),
  extra_columns = c( "Oxygen", "ICU", "Steroid","Daily Dose (mg)", "Duration (days)"),       
  rob_column = "rob2_patterns"      
)


# Bayesian analysis (low dose, short duration) 
ha_bayes <- bayesian(
  outcomes = "Hospital-acquired infections",
  data = dt[dt$outcomes == "Hospital-acquired infections" & !is.na(dt$event_placebo), ],
  study_label_col = "Authors+O44LA1:R41",
  thetaprior = c(0.96, 0.5), # prior Blum et al.
  returnstan = TRUE, # FALSE for posteriori probability et exp(log(OR)) + CrI
  tauprior = 0.5,
  adapt_delta = 0.95,
  return_summary = F
)
# Base


# Bayesian analysis (without restriction on dose, duration or timing)
ha_bayes_primary <- bayesian(
  outcomes = "Hospital-acquired infections",
  data = dt[dt$outcomes == "Hospital-acquired infections" & !is.na(dt$event_placebo)& dt$CI==TRUE, ],
  study_label_col = "Authors+O44LA1:R41",
  thetaprior = c(0.96, 0.5),# prior Blum et al.
  returnstan = TRUE, # FALSE for posteriori probability et exp(log(OR)) + CrI
  tauprior = 0.5,
  adapt_delta = 0.95,
  return_summary = F
)


###############################################################################
###############################################################################
###############################################################################
####                    Choc septique                                        ####
###############################################################################
###############################################################################
###############################################################################

# Primary analysis (low dose, short duration) frequentist : Paule-Mandel + Hartung-Knapp
infection_ss <- analysis_meta(
  outcomes = "Occurrence of septic shock",
  data = dt[dt$outcomes == "Occurrence of septic shock" & !is.na(dt$event_placebo) & dt$CI==TRUE, ]
  [order(dt[dt$outcomes == "Occurrence of septic shock" & !is.na(dt$event_placebo) & dt$CI==TRUE, ]$steroid), ],
  metabin = TRUE,
  use_hartung_knapp = T,
  iqwig = T,
  sm = "RR",
  event.e = "event_steroid",
  n.e = "total_steroid",
  event.c = "event_placebo",
  n.c = "total_placebo",
  studlab = "Authors+O44LA1:R41",
  subgroup = "type2",
  subgroup_levels = c("Non severe pneumonia", "Severe pneumonia", "ARDS"),
  extra_columns = c( "Oxygen", "ICU", "Steroid","Daily Dose (mg)", "Duration (days)"),        
  rob_column = "rob2_patterns"       # AJOUTE LA COLONNE ROB2
)


# Secondary analysis (without restriction on dose, duration or timing) frequentist : Paule-Mandel + Hartung-Knapp
infection_ss_all <- analysis_meta(
  outcomes = "Occurrence of septic shock",
  data = dt[dt$outcomes == "Occurrence of septic shock" & !is.na(dt$event_placebo), ]
  [order(dt[dt$outcomes == "Occurrence of septic shock" & !is.na(dt$event_placebo), ]$steroid), ],
  metabin = TRUE,
  use_hartung_knapp = T,
  iqwig = T,
  sm = "RR",
  event.e = "event_steroid",
  n.e = "total_steroid",
  event.c = "event_placebo",
  n.c = "total_placebo",
  studlab = "Authors+O44LA1:R41",
  subgroup = "type2",
  subgroup_levels = c("Non severe pneumonia", "Severe pneumonia", "ARDS"),
  extra_columns = c( "Oxygen", "ICU", "Steroid","Daily Dose (mg)", "Duration (days)"),        
  rob_column = "rob2_patterns"       # AJOUTE LA COLONNE ROB2
)

# Bayesian analysis (low dose, short duration)  
ss_bayes <- bayesian(
  outcomes = "Occurrence of septic shock",
  data = dt[dt$outcomes == "Occurrence of septic shock" & !is.na(dt$event_placebo), ],
  study_label_col = "Authors+O44LA1:R41",
  thetaprior = c(0.69, 0.5),
  returnstan = TRUE,
  tauprior = 0.3,
  adapt_delta = 0.95,
  return_summary = F
)

# Bayesian analysis (without restriction on dose, duration or timing)
ss_bayes_primary <- bayesian(
  outcomes = "Occurrence of septic shock",
  data = dt[dt$outcomes == "Occurrence of septic shock" & !is.na(dt$event_placebo)& dt$CI==TRUE, ],
  study_label_col = "Authors+O44LA1:R41",
  thetaprior = c(0.69, 0.5),
  returnstan = TRUE,
  tauprior = 0.3,
  adapt_delta = 0.95,
  return_summary = F
)





###############################################################################
###############################################################################
###############################################################################
####                    Pneumonie                                        ####
###############################################################################
###############################################################################
###############################################################################


# Primary analysis (low dose, short duration) frequentist : Paule-Mandel + Hartung-Knapp
infection_sp <- analysis_meta(
  outcomes = "Occurrence of secondary pneumonia",
  data = dt[dt$outcomes == "Occurrence of secondary pneumonia" & !is.na(dt$event_placebo) & dt$CI==TRUE, ]
  [order(dt[dt$outcomes == "Occurrence of secondary pneumonia" & !is.na(dt$event_placebo) & dt$CI==TRUE, ]$steroid), ],
  metabin = TRUE,
  use_hartung_knapp = T,
  iqwig = T,
  sm = "RR",
  event.e = "event_steroid",
  n.e = "total_steroid",
  event.c = "event_placebo",
  n.c = "total_placebo",
  studlab = "Authors+O44LA1:R41",
  subgroup = "type2",
  subgroup_levels = c("Non severe pneumonia", "Severe pneumonia", "ARDS"),
  extra_columns = c( "Oxygen", "ICU", "Steroid","Daily Dose (mg)", "Duration (days)"),       
  rob_column = "rob2_patterns"       # AJOUTE LA COLONNE ROB2
)


# Secondary analysis (without restriction on dose, duration or timing) frequentist : Paule-Mandel + Hartung-Knapp
infection_sp_all <- analysis_meta(
  outcomes = "Occurrence of secondary pneumonia",
  data = dt[dt$outcomes == "Occurrence of secondary pneumonia" & !is.na(dt$event_placebo), ]
  [order(dt[dt$outcomes == "Occurrence of secondary pneumonia" & !is.na(dt$event_placebo), ]$steroid), ],
  metabin = TRUE,
  use_hartung_knapp = T,
  iqwig = T,
  sm = "OR",
  event.e = "event_steroid",
  n.e = "total_steroid",
  event.c = "event_placebo",
  n.c = "total_placebo",
  studlab = "Authors+O44LA1:R41",
  subgroup = "type2",
  subgroup_levels = c("Non severe pneumonia", "Severe pneumonia", "ARDS"),
  extra_columns = c( "Oxygen", "ICU", "Steroid","Daily Dose (mg)", "Duration (days)"),        
  rob_column = "rob2_patterns"       
)



# Bayesian analysis (low dose, short duration) 
sp_bayes <- bayesian(
  outcomes = "Occurrence of secondary pneumonia",
  data = dt[dt$outcomes == "Occurrence of secondary pneumonia" & !is.na(dt$event_placebo), ],
  study_label_col = "Authors+O44LA1:R41",
  thetaprior = c(0.85, 0.5), # priori Snidjers et al.
  returnstan = TRUE,
  tauprior = 0.5,
  adapt_delta = 0.95,
  return_summary = F
)

# Bayesian analysis (without restriction on dose, duration or timing)
sp_bayes_primary <- bayesian(
  outcomes = "Occurrence of secondary pneumonia",
  data = dt[dt$outcomes == "Occurrence of secondary pneumonia" & !is.na(dt$event_placebo)& dt$CI==TRUE, ],
  study_label_col = "Authors+O44LA1:R41",
  thetaprior = c(0.85, 0.5),# priori Snidjers et al.
  returnstan = TRUE,
  tauprior = 0.5,
  adapt_delta = 0.95,
  return_summary = F
)




###############################################################################
###############################################################################
###############################################################################
####                    Catether                                       ####
###############################################################################
###############################################################################
###############################################################################



# Primary analysis (low dose, short duration) frequentist : Paule-Mandel + Hartung-Knapp
infection_cat <- analysis_meta(
  outcomes = "Bloodstream/catheter related related infection (catheter)",
  data = dt[dt$outcomes == "Bloodstream/catheter related related infection (catheter)" & !is.na(dt$event_placebo) & dt$CI==TRUE, ]
  [order(dt[dt$outcomes == "Bloodstream/catheter related related infection (catheter)" & !is.na(dt$event_placebo) & dt$CI==TRUE, ]$steroid), ],
  metabin = TRUE,
  use_hartung_knapp = T,
  iqwig = T,
  sm = "RR",
  event.e = "event_steroid",
  n.e = "total_steroid",
  event.c = "event_placebo",
  n.c = "total_placebo",
  studlab = "Authors+O44LA1:R41",
  subgroup = "type2",
  subgroup_levels = c("Non severe pneumonia", "Severe pneumonia", "ARDS"),
  extra_columns = c( "Oxygen", "ICU", "Steroid","Daily Dose (mg)", "Duration (days)"),        
  rob_column = "rob2_patterns"       
)
# 10 6
# Secondary analysis (without restriction on dose, duration or timing) frequentist : Paule-Mandel + Hartung-Knapp
infection_cat_all <- analysis_meta(
  outcomes = "Bloodstream/catheter related related infection (catheter)",
  data = dt[dt$outcomes == "Bloodstream/catheter related related infection (catheter)" & !is.na(dt$event_placebo), ]
  [order(dt[dt$outcomes == "Bloodstream/catheter related related infection (catheter)" & !is.na(dt$event_placebo) , ]$steroid), ],
  metabin = TRUE,
  use_hartung_knapp = T,
  iqwig = T,
  sm = "RR",
  event.e = "event_steroid",
  n.e = "total_steroid",
  event.c = "event_placebo",
  n.c = "total_placebo",
  studlab = "Authors+O44LA1:R41",
  subgroup = "type2",
  subgroup_levels = c("Non severe pneumonia", "Severe pneumonia", "ARDS"),
  extra_columns = c( "Oxygen", "ICU", "Steroid","Daily Dose (mg)", "Duration (days)"),      
  rob_column = "rob2_patterns"       
)

# Bayesian analysis (low dose, short duration) 
cat_bayes <- bayesian(
  outcomes = "Bloodstream/catheter related related infection (catheter)",
  data = dt[dt$outcomes == "Bloodstream/catheter related related infection (catheter)" & !is.na(dt$event_placebo), ],
  study_label_col = "Authors+O44LA1:R41",
  thetaprior = c(0.44, 0.5),
  returnstan = TRUE,
  tauprior = 0.5,
  adapt_delta = 0.95,
  return_summary = F
)

# Bayesian analysis (without restriction on dose, duration or timing)
cat_bayes_primary <- bayesian(
  outcomes = "Bloodstream/catheter related related infection (catheter)",
  data = dt[dt$outcomes == "Bloodstream/catheter related related infection (catheter)" & !is.na(dt$event_placebo)& dt$CI==TRUE, ],
  study_label_col = "Authors+O44LA1:R41",
  thetaprior = c(0.44, 0.5),
  returnstan = TRUE,
  tauprior = 0.5,
  adapt_delta = 0.95,
  return_summary = F
)


###############################################################################
###############################################################################
###############################################################################
####                    Blood                                       ####
###############################################################################
###############################################################################
###############################################################################

# Primary analysis (low dose, short duration) frequentist : Paule-Mandel + Hartung-Knapp
infection_bld <- analysis_meta(
  outcomes = "Bloodstream/catheter related related infection (bloodstream)",
  data = dt[dt$outcomes == "Bloodstream/catheter related related infection (bloodstream)" & !is.na(dt$event_placebo) & dt$CI==TRUE, ]
  [order(dt[dt$outcomes == "Bloodstream/catheter related related infection (bloodstream)" & !is.na(dt$event_placebo) & dt$CI==TRUE, ]$steroid), ],
  metabin = TRUE,
  use_hartung_knapp = T,
  iqwig = T,
  sm = "RR",
  event.e = "event_steroid",
  n.e = "total_steroid",
  event.c = "event_placebo",
  n.c = "total_placebo",
  studlab = "Authors+O44LA1:R41",
  subgroup = "type2",
  subgroup_levels = c("Non severe pneumonia", "Severe pneumonia", "ARDS"),
  extra_columns = c( "Oxygen", "ICU", "Steroid","Daily Dose (mg)", "Duration (days)"),        
  rob_column = "rob2_patterns"       # AJOUTE LA COLONNE ROB2
)

# Secondary analysis (without restriction on dose, duration or timing) frequentist : Paule-Mandel + Hartung-Knapp
infection_bld_all <- analysis_meta(
  outcomes = "Bloodstream/catheter related related infection (bloodstream)",
  data = dt[dt$outcomes == "Bloodstream/catheter related related infection (bloodstream)" & !is.na(dt$event_placebo), ]
  [order(dt[dt$outcomes == "Bloodstream/catheter related related infection (bloodstream)" & !is.na(dt$event_placebo), ]$steroid), ],
  metabin = TRUE,
  use_hartung_knapp = T,
  iqwig = T,
  sm = "RR",
  event.e = "event_steroid",
  n.e = "total_steroid",
  event.c = "event_placebo",
  n.c = "total_placebo",
  studlab = "Authors+O44LA1:R41",
  subgroup = "type2",
  subgroup_levels = c("Non severe pneumonia", "Severe pneumonia", "ARDS"),
  extra_columns = c( "Oxygen", "ICU", "Steroid","Daily Dose (mg)", "Duration (days)"),       
  rob_column = "rob2_patterns"       
)

# Bayesian analysis (low dose, short duration) 
bld_bayes <- bayesian(
  outcomes = "Bloodstream/catheter related related infection (bloodstream)",
  data = dt[dt$outcomes == "Bloodstream/catheter related related infection (bloodstream)" & !is.na(dt$event_placebo), ],
  study_label_col = "Authors+O44LA1:R41",
  thetaprior = c(0.69, 0.5),
  returnstan = TRUE,
  tauprior = 0.5,
  adapt_delta = 0.95,
  return_summary = F
)

# Bayesian analysis (without restriction on dose, duration or timing)
bld_bayes_primary <- bayesian(
  outcomes = "Bloodstream/catheter related related infection (bloodstream)",
  data = dt[dt$outcomes == "Bloodstream/catheter related related infection (bloodstream)" & !is.na(dt$event_placebo) &dt$CI==TRUE, ],
  study_label_col = "Authors+O44LA1:R41",
  thetaprior = c(0.69, 0.5),
  returnstan = TRUE,
  tauprior = 0.5,
  adapt_delta = 0.95,
  return_summary = F
)





###############################################################################
###############################################################################
###############################################################################
####                    DISTRIBUTION                                       ####
###############################################################################
###############################################################################
###############################################################################


# ---- Bayesian list ----
liste_secondary <- list(ha_bayes, ss_bayes, sp_bayes, cat_bayes, bld_bayes)
liste_primary   <- list(ha_bayes_primary, ss_bayes_primary, sp_bayes_primary, cat_bayes_primary, bld_bayes_primary)

prior_labels <- c(
  "Hospital-acquired infections",
  "Secondary shock",
  "Secondary pneumonia",
  "Catheter-related infection",
  "Bloodstream infection"
)

# ---- Création des plots ----
bayesian_plot         <- make_prior_plots(liste_secondary, prior_labels,
                                     "Prior and Posterior distributions for Bayesian analyses")

bayesian_plot_primary <- make_prior_plots(liste_primary, prior_labels,
                                     "Low dose, short duration : Prior and Posterior distributions for Bayesian analyses")


######################################################################################################################################################################################################################
######################################################################################################################################################################################################################
######################################################################################################################################################################################################################
######################################################################################################################################################################################################################
######################################################################################################################################################################################################################
######################################################################################################################################################################################################################
###############################################################################
###############################################################################
####                       Mortality Analysis                              ####
###############################################################################
###############################################################################
######################################################################################################################################################################################################################
######################################################################################################################################################################################################################
######################################################################################################################################################################################################################
######################################################################################################################################################################################################################
######################################################################################################################################################################################################################
######################################################################################################################################################################################################################


###############################################################################
# Low dose low duration (primary analysis)
# Paule-Manle + Hartung-Knapp
###############################################################################

# Filter primary data
primary <- dt[
  dt$outcomes == "Short-term mortality" & 
  !is.na(dt$event_placebo) & 
  dt$CI == TRUE & # CI = low dose, short duration,timing < 7 days
  dt$duplicata == 0, 
]

# Remove NA in 'type2' and sort by steroid (type2 = "ARDS", or "severe pneumonia" or "non severe pneumonia")
primary <- primary[!is.na(primary$type2), ]
primary <- primary[order(primary$Steroid), ]

# Main primary analysis
primary_analysis <- analysis_meta(
  outcomes = "Short-term mortality",
  data = primary,
  metabin = TRUE,
  use_hartung_knapp = T,
  iqwig = T,
  sm = "RR",
  event.e = "event_steroid",
  n.e     = "total_steroid",
  event.c = "event_placebo",
  n.c     = "total_placebo",
  studlab = "Authors+O44LA1:R41",
  subgroup = "type2",
  subgroup_remove = "",
  exclusion_condition = NULL,
  subgroup_levels = c("Non severe pneumonia", "Severe pneumonia", "ARDS"),
  extra_columns = c( "Oxygen", "ICU", "Steroid","Daily Dose (mg)", "Duration (days)"),        
  rob_column = "rob2_patterns"       
)


###############################################################################
# All studies (secondary analysis)
###############################################################################

# Filter secondary data
secondary <- dt[
  dt$outcomes == "Short-term mortality" & 
  !is.na(dt$event_placebo) & 
  dt$duplicata < 0.8, # information for remove duplicate 
]

# Remove NA and sort
secondary <- secondary[!is.na(secondary$type2), ]
secondary <- secondary[order(secondary$Steroid), ]

secondary_analysis <- analysis_meta(
  outcomes = "Short-term mortality",
  data = secondary,
  metabin = TRUE,
  use_hartung_knapp = T,
  iqwig = T,
  sm = "RR",
  event.e = "event_steroid",
  n.e     = "total_steroid",
  event.c = "event_placebo",
  n.c     = "total_placebo",
  studlab = "Authors+O44LA1:R41",
  subgroup = "type2",
  subgroup_remove = "",
  exclusion_condition = NULL,
  subgroup_levels = c("Non severe pneumonia", "Severe pneumonia", "ARDS"),
  extra_columns = c( "Oxygen", "ICU", "Steroid","Daily Dose (mg)", "Duration (days)"),        
  rob_column = "rob2_patterns"       
)

###############################################################################
# Risk of bias analysis (Rob2 = 1 only)
###############################################################################

# Paule-Mandel + Hartung-Knapp
# Primary low risk of bias
lowrisk <- analysis_meta(
  outcomes = "Short-term mortality",
  data = primary[primary$Rob2 == 1,], # selection of low risk of bias study
  metabin = TRUE,
  use_hartung_knapp = T,
  iqwig = T,
  sm = "RR",
  event.e = "event_steroid",
  n.e     = "total_steroid",
  event.c = "event_placebo",
  n.c     = "total_placebo",
  studlab = "Authors+O44LA1:R41",
  subgroup = "type2",
  subgroup_remove = "Rob2",
  exclusion_condition = 0,
  subgroup_levels = c("Non severe pneumonia", "Severe pneumonia", "ARDS"),
  extra_columns = c( "Oxygen", "ICU", "Steroid","Daily Dose (mg)", "Duration (days)"),        
  rob_column = "rob2_patterns"       
)

# Secondary low risk of bias
secondarylowrisk <- analysis_meta(
  outcomes = "Short-term mortality",
  data = secondary[secondary$Rob2 == 1,], # selection of low risk of bias study
  metabin = TRUE,
  use_hartung_knapp = T,
  iqwig = T,
  sm = "RR",
  event.e = "event_steroid",
  n.e     = "total_steroid",
  event.c = "event_placebo",
  n.c     = "total_placebo",
  studlab = "Authors+O44LA1:R41",
  subgroup = "type2",
  subgroup_remove = "Rob2",
  exclusion_condition = 0,
  subgroup_levels = c("Non severe pneumonia", "Severe pneumonia", "ARDS"),
  extra_columns = c( "Oxygen", "ICU", "Steroid","Daily Dose (mg)", "Duration (days)"),        
  rob_column = "rob2_patterns"       
)

###############################################################################
# DerSimonian-Laird method (DL)
###############################################################################

# Primary with DL method
primary_analysis_dl <- analysis_meta(
  outcomes = "Short-term mortality",
  data = primary,
  metabin = TRUE,
  use_hartung_knapp = F,
  iqwig = F,
  sm = "RR",
  event.e = "event_steroid",
  n.e     = "total_steroid",
  event.c = "event_placebo",
  n.c     = "total_placebo",
  studlab = "Authors+O44LA1:R41",
  subgroup = "type2",
  subgroup_remove = "",
  exclusion_condition = NULL,
  subgroup_levels = c("Non severe pneumonia", "Severe pneumonia", "ARDS"),
  extra_columns = c( "Oxygen", "ICU", "Steroid","Daily Dose (mg)", "Duration (days)"),        
  rob_column = "rob2_patterns"       
)



# All studies with DL method
sensitivity_analysis_all_dl <- analysis_meta(
  outcomes = "Short-term mortality",
  data = secondary,
  metabin = TRUE,
  use_hartung_knapp = F,
  iqwig = F,
  sm = "RR",
  event.e = "event_steroid",
  n.e     = "total_steroid",
  event.c = "event_placebo",
  n.c     = "total_placebo",
  studlab = "Authors+O44LA1:R41",
  subgroup = "type2",
  subgroup_remove = "",
  exclusion_condition = NULL,
  subgroup_levels = c("Non severe pneumonia", "Severe pneumonia", "ARDS"),
  extra_columns = c( "Oxygen", "ICU", "Steroid","Daily Dose (mg)", "Duration (days)"),        
  rob_column = "rob2_patterns"       
)

###############################################################################
# Bayesian analyses
###############################################################################

# Prior = estimation of REMAP-CAP study
# Secondary ARDS bayesian
secondary_bayesien_ards <- bayesian(outcomes = "Short-term mortality", 
                                   data = secondary[secondary$type2 == "ARDS", ],
                                   study_label_col = "Authors+O44LA1:R41",
                                   thetaprior = c(0.49, 0.5), # priori REMAP-CAP
                                   returnstan = T,
                                   tauprior = 0.5,
                                   adapt_delta = 0.95, 
                                   return_summary = T)

# Secondary severe pneumonia bayesian                                   
secondary_bayesien_severe <- bayesian(outcomes = "Short-term mortality", 
                                     data = secondary[secondary$type2 == "Severe pneumonia", ],
                                     study_label_col = "Authors+O44LA1:R41",
                                     thetaprior = c(0.49, 0.5), # priori REMAP-CAP
                                     returnstan = T,
                                     tauprior = 0.5,
                                     adapt_delta = 0.95, 
                                     return_summary = T)

# Secondary non-severe pneumonia bayesian
secondary_bayesien_nonsevere <- bayesian(outcomes = "Short-term mortality", 
                                        data = secondary[secondary$type2 == "Non severe pneumonia", ],
                                        study_label_col = "Authors+O44LA1:R41",
                                        thetaprior = c(0.405, 0.5), # priori REMAP-CAP
                                        returnstan = T,
                                        tauprior = 0.5,
                                        adapt_delta = 0.95, 
                                        return_summary = T)

# Primary ARDS bayesian
primary_bayesien_ards <- bayesian(outcomes = "Short-term mortality", 
                                 data = primary[primary$type2 == "ARDS", ],
                                 study_label_col = "Authors+O44LA1:R41",
                                 thetaprior = c(0.49, 0.5), # priori REMAP-CAP
                                 returnstan = T,
                                 tauprior = 0.5,
                                 adapt_delta = 0.99, 
                                 return_summary = F)

# Primary severe pneumonia bayesian                     
primary_bayesien_severe <- bayesian(outcomes = "Short-term mortality", 
                                   data = primary[primary$type2 == "Severe pneumonia", ],
                                   study_label_col = "Authors+O44LA1:R41",
                                   thetaprior = c(0.49, 0.5), # priori REMAP-CAP
                                   returnstan = T,
                                   tauprior = 0.5,
                                   adapt_delta = 0.95, 
                                   return_summary = T)


# Other function for plots
make_prior_plots2 <- function(liste, prior_labels, titre) {
  plots <- list()
  
  for (i in seq_along(liste)) {
    element <- liste[[i]]
    
    # Extraction des échantillons et paramètres du prior
    theta_samples <- rstan::extract(element$fit, pars = "theta", permuted = TRUE)$theta
    prior_mean <- element$stanDat$theta_prior[1]
    prior_sd <- sqrt(element$stanDat$theta_prior[2])
    
    # Vérification des paramètres
    if (!all(is.finite(theta_samples))) {
      message("Skipping element ", i, ": Invalid posterior samples")
      next
    }
    if (!is.finite(prior_mean) || !is.finite(prior_sd) || prior_sd <= 0) {
      message("Skipping element ", i, ": Invalid prior parameters (mean = ", prior_mean, 
              ", sd = ", prior_sd, ")")
      next
    }
    
    # Création du graphique avec la légende
    plot <- plot_prior_posterior(prior_mean, prior_sd, theta_samples, prior_labels[i])
    plots[[i]] <- plot
  }
  
  # Assemblage patchwork selon le nombre de plots
  if (length(plots) == 2) {
    big_plot <- plots[[1]] + plots[[2]]
  } else if (length(plots) == 3) {
    big_plot <- (plots[[1]] + plots[[2]]) /
                (plots[[3]] + plot_spacer())
  } else {
    stop("⚠️ La fonction attend 2 ou 3 plots, pas ", length(plots))
  }
  
  big_plot <- big_plot +
    plot_layout(heights = c(1, 1)) +
    plot_annotation(
      title = titre,
      theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
    )
  
  return(big_plot)

}


big_plot         <- make_prior_plots2(secondary_bayes_plots, prior_labels,
                                     "Prior and Posterior distributions for Bayesian analyses")

big_plot_primary <- make_prior_plots2(# ---- Liste des analyses bayésiennes primaires ----
primary_bayes_plots , prior_labels,
                                     "Low dose, short duration : Prior and Posterior distributions for Bayesian analyses")



######################################################################################################################################################################################################################
######################################################################################################################################################################################################################
######################################################################################################################################################################################################################
######################################################################################################################################################################################################################
######################################################################################################################################################################################################################
######################################################################################################################################################################################################################
###############################################################################
###############################################################################
####                       Exploration heterogeneity                              ####
###############################################################################
###############################################################################
######################################################################################################################################################################################################################
######################################################################################################################################################################################################################
######################################################################################################################################################################################################################
######################################################################################################################################################################################################################
######################################################################################################################################################################################################################
######################################################################################################################################################################################################################


#########################################################################################################################################################
######### Leave one out analysis
#########################################################################################################################################################

plot_loo_subgroup <- function(meta_object, subgroup_name, return = T) {
  m_sub <- update(meta_object, subset = (meta_object$byvar == subgroup_name))
  inf <- metainf(m_sub, pooled = "random")
  
  if (return) {
    # Return forest plot directly
    return(
      forest(inf, 
             col.square = "black", 
             col.square.lines = "black", 
             col.diamond = "#2E8B7A", 
             col.diamond.lines = "#2E8B7A")
    )
  } else {
    # Return grob for save_forest
    return(
      gridGraphics::echoGrob(function() {
        forest(inf, 
               col.square = "black", 
               col.square.lines = "black", 
               col.diamond = "#2E8B7A", 
               col.diamond.lines = "#2E8B7A")
      })
    )
  }
}

# Generate leave-one-out plots
loo_primary_severe <- plot_loo_subgroup(primary_analysis$meta, "Severe pneumonia")
loo_primary_ARDS <- plot_loo_subgroup(primary_analysis$meta, "ARDS")
loo_secondary_severe <- plot_loo_subgroup(secondary_analysis$meta, "Severe pneumonia")
loo_secondary_ARDS <- plot_loo_subgroup(secondary_analysis$meta, "ARDS")

#########################################################################################################################################################
######### Cumulative plot
#########################################################################################################################################################

cumu_plot <- function(meta, var, val, return_only = TRUE) {
  keep_indices <- which(meta$data[[var]] == val)
  years_numeric <- as.numeric(meta$data$years[keep_indices])
  sorted_order <- order(years_numeric)
  final_indices <- keep_indices[sorted_order]
  meta_sub <- update(meta, subset = final_indices)
  cum <- metacum(meta_sub, pooled = "random")
  
  if (return_only) {
    class(cum) <- c("auto_plot", class(cum))
    return(cum)
  } else {
    forest(cum, col.square = "black", col.square.lines = "black", 
           col.diamond = "#2E8B7A", col.diamond.lines = "#2E8B7A")
  }
}

# Print method for auto plots
print.auto_plot <- function(x, ...) {
  forest(x, col.square = "black", col.square.lines = "black", 
         col.diamond = "#2E8B7A", col.diamond.lines = "#2E8B7A")
}

# Generate cumulative plots
cum_primary_severe <- cumu_plot(primary_analysis$meta, "type2", "Severe pneumonia")
cum_primary_ARDS <- cumu_plot(primary_analysis$meta, "type2", "ARDS")
cum_secondary_severe <- cumu_plot(secondary_analysis$meta, "type2", "Severe pneumonia")
cum_secondary_ARDS <- cumu_plot(secondary_analysis$meta, "type2", "ARDS")

#########################################################################################################################################################
######### Baujat plot
#########################################################################################################################################################

plot_baujat <- function(meta_object, subgroup_vector, value, return_only = T) {
  m <- update(meta_object, subset = (subgroup_vector == value))
  
  if (return_only == T) {
    # Store meta object, not baujat result
    class(m) <- c("auto_plot_baujat", class(m))
    return(m)
  } else {
    baujat(m,
           col = "#2E8B7A",
           pch = 19,
           cex = 1.2,
           las = 1)
  }
}

print.auto_plot_baujat <- function(x, ...) {
    baujat(x, 
           col = "#2E8B7A", 
           pch = 19, 
           cex = 1.2, 
           las = 1)
}

# Generate Baujat plots
bj_primary_severe <- plot_baujat(primary_analysis$meta, primary$type2, "Severe pneumonia")
bj_primary_ARDS <- plot_baujat(primary_analysis$meta, primary$type2, "ARDS")
bj_secondary_severe <- plot_baujat(secondary_analysis$meta, secondary$type2, "Severe pneumonia")
bj_secondary_ARDS <- plot_baujat(secondary_analysis$meta, secondary$type2, "ARDS")

#########################################################################################################################################################
# Meta-regression analyses
#########################################################################################################################################################

# Helper function for bubble plots
create_bubble_plot <- function(metareg_obj, xlim = NULL, ylim = NULL, return = FALSE) {
  if (return) {
    # Return bubble plot directly
    return(
      meta::bubble(
        metareg_obj,
        xlim = xlim,
        ylim = ylim,
        pch = 21,              
        col = "#2E8B7A",       
        bg  = "#2E8B7A",       
        col.line = "#2E8B7A",  
        studlab = TRUE
      )
    )
  } else {
    # Return grob for save_forest
    return(
      gridGraphics::echoGrob(function() {
        meta::bubble(
          metareg_obj,
          xlim = xlim,
          ylim = ylim,
          pch = 21,
          col = "#2E8B7A",
          bg  = "#2E8B7A",
          col.line = "#2E8B7A",
          studlab = TRUE
        )
      })
    )
  }
}

# Primary analysis meta-regressions
primary$type2_bin <- ifelse(primary$type2 == "Severe pneumonia", 1, 0)
primary$type2_factor <- factor(primary$type2, levels = c("ARDS", "Severe pneumonia"))
primary_analysis$meta$data$type2_factor <- primary$type2_factor
meta_reg_primary_group <- metareg(primary_analysis$meta, ~ type2_factor)

# Create bubble plots for primary analysis
bubble_reg_primary_group <- create_bubble_plot(meta_reg_primary_group, xlim = c(-1, 2), ylim = c(0.3, 3), return = T)

# Meta regression dose - severe pneumonia
m_severe <- update(primary_analysis$meta, subset = (primary$type2 == "Severe pneumonia"))
meta_reg_severe_dose <- metareg(m_severe, ~ dose_moyenne)
bubble_reg_severe_dose <- create_bubble_plot(meta_reg_severe_dose, xlim = c(0, 100), ylim = c(0.3, 3), return = T)

# Meta regression dose - ARDS
m_ARDS <- update(primary_analysis$meta, subset = (primary$type2 == "ARDS"))
meta_reg_ARDS_dose <- metareg(m_ARDS, ~ dose_moyenne)
bubble_reg_ARDS_dose <- create_bubble_plot(meta_reg_ARDS_dose, xlim = c(0, 100), ylim = c(0.3, 3), return = T)

# Meta regression duration - severe pneumonia
meta_reg_severe_duree <- metareg(m_severe, ~ duration_max)
bubble_reg_severe_duree <- create_bubble_plot(meta_reg_severe_duree, xlim = c(-8, 30), ylim = c(0.3, 3),  return = T)

# Meta regression duration - ARDS
meta_reg_ARDS_duree <- metareg(m_ARDS, ~ duration_max)
bubble_reg_ARDS_duree <- create_bubble_plot(meta_reg_ARDS_duree, xlim = c(-8, 30), ylim = c(0.3, 3),  return = T)

# Interaction analysis
complete_cases <- complete.cases(primary_analysis$meta$TE, primary_analysis$data$type2)
primary_interaction <- update(primary_analysis$meta, 
       TE = primary_analysis$meta$TE[complete_cases],
       seTE = primary_analysis$meta$seTE[complete_cases],
       studlab = primary_analysis$meta$studlab[complete_cases],
       data = primary_analysis$meta$data[complete_cases, ],
       subgroup = primary_analysis$meta$data$type2[complete_cases])

#########################################################################################################################################################
# Secondary analysis meta-regression
#########################################################################################################################################################

# Remove Mikami 2007 study
secondary_clean <- secondary[!secondary$`Authors+O44LA1:R41` == "Mikami 2007",] # 0 events for MIKAMI et al. 

secondary_analysis_without_mikami <- analysis_meta( # 0 events for MIKAMI et al. 
  outcomes = "Short-term mortality",
  data = secondary_clean,
  metabin = TRUE,
  use_hartung_knapp = T,
  iqwig = T,
  sm = "RR",
  event.e = "event_steroid",
  n.e     = "total_steroid",
  event.c = "event_placebo",
  n.c     = "total_placebo",
  studlab = "Authors+O44LA1:R41",
  subgroup = "type2",
  subgroup_remove = "",
  exclusion_condition = NULL,
  subgroup_levels = c("Non severe pneumonia", "Severe pneumonia", "ARDS"),
  extra_columns = c( "Oxygen", "ICU", "Steroid","Daily Dose (mg)", "Duration (days)"),        
  rob_column = "rob2_patterns"       
)

# Meta regression secondary - subgroups
secondary_clean$type2_factor <- factor(secondary_clean$type2, levels = c("ARDS", "Severe pneumonia", "Non severe pneumonia"))
secondary_analysis_without_mikami$meta$data$type2_factor <- secondary_clean$type2_factor
meta_reg_secondary_group <- metareg(secondary_analysis_without_mikami$meta, ~ type2_factor)
bubble_reg_secondary_group <- create_bubble_plot(meta_reg_secondary_group, xlim = c(-1, 1), ylim = c(0.3, 3), return =T)

# Meta regression dose secondary - severe pneumonia
m_severe_secondary <- update(secondary_analysis_without_mikami$meta, subset = (secondary_clean$type2 == "Severe pneumonia"))
meta_reg_secondary_severe_dose <- metareg(m_severe_secondary, ~ dose_moyenne)
bubble_reg_secondary_severe_dose <- create_bubble_plot(meta_reg_secondary_severe_dose, xlim = c(0, 100), ylim = c(0.3, 3), return =T)

# Meta regression dose secondary - ARDS
m_ARDS_secondary <- update(secondary_analysis_without_mikami$meta, subset = (secondary_clean$type2 == "ARDS"))
meta_reg_secondary_ARDS_dose <- metareg(m_ARDS_secondary, ~ dose_moyenne)
bubble_reg_secondary_ARDS_dose <- create_bubble_plot(meta_reg_secondary_ARDS_dose, xlim = c(0, 100), ylim = c(0.3, 3))

# Meta regression duration secondary - severe pneumonia  
m_severe_secondary_duree <- update(secondary_analysis_without_mikami$meta, subset = (secondary_clean$type2 == "Severe pneumonia"))
meta_reg_secondary_severe_duree <- metareg(m_severe_secondary_duree, ~ duration_max)
bubble_reg_secondary_severe_duree <- create_bubble_plot(meta_reg_secondary_severe_duree, xlim = c(-8, 30), ylim = c(0.3, 3))

# Meta regression duration secondary - ARDS
m_ARDS_secondary_duree <- update(secondary_analysis_without_mikami$meta, subset = (secondary_clean$type2 == "ARDS"))
meta_reg_secondary_ARDS_duree <- metareg(m_ARDS_secondary_duree, ~ duration_max)
bubble_reg_secondary_ARDS_duree <- create_bubble_plot(meta_reg_secondary_ARDS_duree, xlim = c(-8, 30), ylim = c(0.3, 3))

# Interaction analysis for secondary
complete_cases_secondary <- complete.cases(secondary_analysis_without_mikami$meta$TE, secondary_analysis_without_mikami$data$type2)
secondary_interaction <- update(secondary_analysis_without_mikami$meta, 
       TE = secondary_analysis_without_mikami$meta$TE[complete_cases_secondary],
       seTE = secondary_analysis_without_mikami$meta$seTE[complete_cases_secondary],
       studlab = secondary_analysis_without_mikami$meta$studlab[complete_cases_secondary],
       data = secondary_analysis_without_mikami$meta$data[complete_cases_secondary, ],
       subgroup = secondary_analysis_without_mikami$meta$data$type2[complete_cases_secondary])



######################################################################################################################################################################################################################
######################################################################################################################################################################################################################
######################################################################################################################################################################################################################
######################################################################################################################################################################################################################
######################################################################################################################################################################################################################
######################################################################################################################################################################################################################
###############################################################################
###############################################################################
####                       Exploration heterogeneity                              ####
###############################################################################
###############################################################################
######################################################################################################################################################################################################################
######################################################################################################################################################################################################################
######################################################################################################################################################################################################################
######################################################################################################################################################################################################################
######################################################################################################################################################################################################################



###############################################################################
###############################################################################
# Subgroup mortality analyses
###############################################################################
###############################################################################

###############################################################################
# Shock-related analyses
###############################################################################

# Primary analysis excluding baseline shock patients
primary_shock <- analysis_meta(
  outcomes = "Short-term mortality",
  data = primary[(primary$shock <= 0.6 | is.na(primary$shock)),],
  metabin = TRUE,
  use_hartung_knapp = T,
  iqwig = T,
  sm = "RR",
  event.e = "event_steroid",
  n.e     = "total_steroid",
  event.c = "event_placebo",
  n.c     = "total_placebo",
  studlab = "Authors+O44LA1:R41",
  subgroup = "type2",
  subgroup_remove = "",
  exclusion_condition = NULL,
  subgroup_levels = c("Non severe pneumonia", "Severe pneumonia", "ARDS"),
  extra_columns = c( "Oxygen", "ICU", "Steroid","Daily Dose (mg)", "Duration (days)"),        
  rob_column = "rob2_patterns"       
)

# Secondary analysis excluding baseline shock patients
secondary_shock <- analysis_meta(
  outcomes = "Short-term mortality",
  data = secondary[(secondary$shock <= 0.6 | is.na(secondary$shock)),],
  metabin = TRUE,
  use_hartung_knapp = T,
  iqwig = T,
  sm = "RR",
  event.e = "event_steroid",
  n.e     = "total_steroid",
  event.c = "event_placebo",
  n.c     = "total_placebo",
  studlab = "Authors+O44LA1:R41",
  subgroup = "type2",
  subgroup_remove = "",
  exclusion_condition = NULL,
  subgroup_levels = c("Non severe pneumonia", "Severe pneumonia", "ARDS"),
  extra_columns = c( "Oxygen", "ICU", "Steroid","Daily Dose (mg)", "Duration (days)"),        
  rob_column = "rob2_patterns"       
)

# Pneumonia-only analysis
only_pneumonia_names <- c(
  "Heming 2023. (sub group with ARDS)",
  "Tongyoo 2016 (only pneumonia)"
)

only_pneumonia <- analysis_meta(
  outcomes = "Short-term mortality",
  data = dt[dt$`Authors+O44LA1:R41`%in% only_pneumonia_names,],
  metabin = TRUE,
  use_hartung_knapp = T,
  iqwig = T,
  sm = "RR",
  event.e = "event_steroid",
  n.e     = "total_steroid",
  event.c = "event_placebo",
  n.c     = "total_placebo",
  studlab = "Authors+O44LA1:R41",
  subgroup = "type2",
  subgroup_remove = "",
  exclusion_condition = NULL,
  subgroup_levels = c("Non severe pneumonia", "Severe pneumonia", "ARDS"),
  extra_columns = c( "Oxygen", "ICU", "Steroid","Daily Dose (mg)", "Duration (days)"),        
  rob_column = "rob2_patterns"       
)

###############################################################################
# Steroid type analyses - PRIMARY
###############################################################################

# Hydrocortisone
hydrocortisone_primary <- primary[primary$steroid == "HC",]
hydrocortisone_analysis_primary <- analysis_meta(
  outcomes = "Short-term mortality",
  data = hydrocortisone_primary,
  metabin = TRUE,
  use_hartung_knapp = T,
  iqwig = T,
  sm = "RR",
  event.e = "event_steroid",
  n.e     = "total_steroid",
  event.c = "event_placebo",
  n.c     = "total_placebo",
  studlab = "Authors+O44LA1:R41",
  subgroup = "type2",
  subgroup_remove = "",
  exclusion_condition = NULL,
  subgroup_levels = c("Non severe pneumonia", "Severe pneumonia", "ARDS"),
  extra_columns = c( "Oxygen", "ICU", "Steroid","Daily Dose (mg)", "Duration (days)"),
  rob_column = "rob2_patterns"
)

# Prednisone
prednisone_primary <- primary[primary$steroid == "PN", ]
prednisone_analysis_primary <- analysis_meta(
  outcomes = "Short-term mortality",
  data = prednisone_primary,
  metabin = TRUE,
  use_hartung_knapp = TRUE,
  iqwig = TRUE,
  sm = "RR",
  event.e = "event_steroid",
  n.e     = "total_steroid",
  event.c = "event_placebo",
  n.c     = "total_placebo",
  studlab = "Authors+O44LA1:R41",
  subgroup = "type2",
  subgroup_remove = "",
  exclusion_condition = NULL,
  subgroup_levels = c("Non severe pneumonia", "Severe pneumonia", "ARDS"),
  extra_columns = c( "Oxygen", "ICU", "Steroid","Daily Dose (mg)", "Duration (days)"),
  rob_column = "rob2_patterns"
)

# Methylprednisolone
methylprednisone_primary <- primary[primary$steroid == "MP", ]
methylprednisone_analysis_primary <- analysis_meta(
  outcomes = "Short-term mortality",
  data = methylprednisone_primary,
  metabin = TRUE,
  use_hartung_knapp = TRUE,
  iqwig = TRUE,
  sm = "RR",
  event.e = "event_steroid",
  n.e     = "total_steroid",
  event.c = "event_placebo",
  n.c     = "total_placebo",
  studlab = "Authors+O44LA1:R41",
  subgroup = "type2",
  subgroup_remove = "",
  exclusion_condition = NULL,
  subgroup_levels = c("Non severe pneumonia", "Severe pneumonia", "ARDS"),
  extra_columns = c( "Oxygen", "ICU", "Steroid","Daily Dose (mg)", "Duration (days)"),
  rob_column = "rob2_patterns"
)

# Dexamethasone
dexamethasone_primary <- primary[primary$steroid == "DXM", ]
dexamethasone_analysis_primary <- analysis_meta(
  outcomes = "Short-term mortality",
  data = dexamethasone_primary,
  metabin = TRUE,
  use_hartung_knapp = TRUE,
  iqwig = TRUE,
  sm = "RR",
  event.e = "event_steroid",
  n.e     = "total_steroid",
  event.c = "event_placebo",
  n.c     = "total_placebo",
  studlab = "Authors+O44LA1:R41",
  subgroup = "type2",
  subgroup_remove = "",
  exclusion_condition = NULL,
  subgroup_levels = c("Non severe pneumonia", "Severe pneumonia", "ARDS"),
  extra_columns = c( "Oxygen", "ICU", "Steroid","Daily Dose (mg)", "Duration (days)"),
  rob_column = "rob2_patterns"
)

###############################################################################
# Duration analyses - PRIMARY
###############################################################################

# Duration ≤7 days
inf_seven_primary <- primary[primary$duration_max <= 7, ]
primary_analysis_le7days <- analysis_meta(
  outcomes = "Short-term mortality",
  data = inf_seven_primary,
  metabin = TRUE,
  use_hartung_knapp = T,
  iqwig = T,
  sm = "RR",
  event.e = "event_steroid",
  n.e     = "total_steroid",
  event.c = "event_placebo",
  n.c     = "total_placebo",
  studlab = "Authors+O44LA1:R41",
  subgroup = "type2",
  subgroup_remove = "",
  exclusion_condition = NULL,
  subgroup_levels = c("Non severe pneumonia", "Severe pneumonia", "ARDS"),
  extra_columns = c( "Oxygen", "ICU", "Steroid","Daily Dose (mg)", "Duration (days)"),
  rob_column = "rob2_patterns"
)

# Duration >7 days
sup_seven_primary <- primary[primary$duration_max > 7, ]
primary_analysis_gt7days <- analysis_meta(
  outcomes = "Short-term mortality",
  data = sup_seven_primary,
  metabin = TRUE,
  use_hartung_knapp = T,
  iqwig = T,
  sm = "RR",
  event.e = "event_steroid",
  n.e     = "total_steroid",
  event.c = "event_placebo",
  n.c     = "total_placebo",
  studlab = "Authors+O44LA1:R41",
  subgroup = "type2",
  subgroup_remove = "",
  exclusion_condition = NULL,
  subgroup_levels = c("Non severe pneumonia", "Severe pneumonia", "ARDS"),
  extra_columns = c( "Oxygen", "ICU", "Steroid","Daily Dose (mg)", "Duration (days)"),
  rob_column = "rob2_patterns"
)

###############################################################################
# Dose analyses - PRIMARY
###############################################################################

# Low dose ≤50mg
dose_basse_primary <- primary[primary$dose_moyenne <= 50, ]
primary_analysis_le50mg <- analysis_meta(
  outcomes = "Short-term mortality",
  data = dose_basse_primary,
  metabin = TRUE,
  use_hartung_knapp = T,
  iqwig = T,
  sm = "RR",
  event.e = "event_steroid",
  n.e     = "total_steroid",
  event.c = "event_placebo",
  n.c     = "total_placebo",
  studlab = "Authors+O44LA1:R41",
  subgroup = "type2",
  subgroup_remove = "",
  exclusion_condition = NULL,
  subgroup_levels = c("Non severe pneumonia", "Severe pneumonia", "ARDS"),
  extra_columns = c( "Oxygen", "ICU", "Steroid","Daily Dose (mg)", "Duration (days)"),
  rob_column = "rob2_patterns"
)

# High dose >50mg
dose_haute_primary <- primary[primary$dose_moyenne > 50, ]
primary_analysis_gt50mg <- analysis_meta(
  outcomes = "Short-term mortality",
  data = dose_haute_primary,
  metabin = TRUE,
  use_hartung_knapp = T,
  iqwig = T,
  sm = "RR",
  event.e = "event_steroid",
  n.e     = "total_steroid",
  event.c = "event_placebo",
  n.c     = "total_placebo",
  studlab = "Authors+O44LA1:R41",
  subgroup = "type2",
  subgroup_remove = "",
  exclusion_condition = NULL,
  subgroup_levels = c("Non severe pneumonia", "Severe pneumonia", "ARDS"),
  extra_columns = c( "Oxygen", "ICU", "Steroid","Daily Dose (mg)", "Duration (days)"),
  rob_column = "rob2_patterns"
)

###############################################################################
# Steroid type analyses - SECONDARY
###############################################################################

# Hydrocortisone
hydrocortisone_secondary <- secondary[secondary$steroid == "HC", ]
hydrocortisone_analysis_secondary <- analysis_meta(
  outcomes = "Short-term mortality",
  data = hydrocortisone_secondary,
  metabin = TRUE,
  use_hartung_knapp = TRUE,
  iqwig = TRUE,
  sm = "RR",
  event.e = "event_steroid",
  n.e     = "total_steroid",
  event.c = "event_placebo",
  n.c     = "total_placebo",
  studlab = "Authors+O44LA1:R41",
  subgroup = "type2",
  subgroup_remove = "",
  exclusion_condition = NULL,
  subgroup_levels = c("Non severe pneumonia", "Severe pneumonia", "ARDS"),
  extra_columns = c( "Oxygen", "ICU", "Steroid","Daily Dose (mg)", "Duration (days)"),
  rob_column = "rob2_patterns"
)

# Prednisone
prednisone_secondary <- secondary[secondary$steroid == "PN", ]
prednisone_analysis_secondary <- analysis_meta(
  outcomes = "Short-term mortality",
  data = prednisone_secondary,
  metabin = TRUE,
  use_hartung_knapp = TRUE,
  iqwig = TRUE,
  sm = "RR",
  event.e = "event_steroid",
  n.e     = "total_steroid",
  event.c = "event_placebo",
  n.c     = "total_placebo",
  studlab = "Authors+O44LA1:R41",
  subgroup = "type2",
  subgroup_remove = "",
  exclusion_condition = NULL,
  subgroup_levels = c("Non severe pneumonia", "Severe pneumonia", "ARDS"),
  extra_columns = c( "Oxygen", "ICU", "Steroid","Daily Dose (mg)", "Duration (days)"),
  rob_column = "rob2_patterns"
)

# Methylprednisolone
methylprednisone_secondary <- secondary[secondary$steroid == "MP", ]
methylprednisone_analysis_secondary <- analysis_meta(
  outcomes = "Short-term mortality",
  data = methylprednisone_secondary,
  metabin = TRUE,
  use_hartung_knapp = TRUE,
  iqwig = TRUE,
  sm = "RR",
  event.e = "event_steroid",
  n.e     = "total_steroid",
  event.c = "event_placebo",
  n.c     = "total_placebo",
  studlab = "Authors+O44LA1:R41",
  subgroup = "type2",
  subgroup_remove = "",
  exclusion_condition = NULL,
  subgroup_levels = c("Non severe pneumonia", "Severe pneumonia", "ARDS"),
  extra_columns = c( "Oxygen", "ICU", "Steroid","Daily Dose (mg)", "Duration (days)"),
  rob_column = "rob2_patterns"
)

# Dexamethasone
dexamethasone_secondary <- secondary[secondary$steroid == "DXM", ]
dexamethasone_analysis_secondary <- analysis_meta(
  outcomes = "Short-term mortality",
  data = dexamethasone_secondary,
  metabin = TRUE,
  use_hartung_knapp = TRUE,
  iqwig = TRUE,
  sm = "RR",
  event.e = "event_steroid",
  n.e     = "total_steroid",
  event.c = "event_placebo",
  n.c     = "total_placebo",
  studlab = "Authors+O44LA1:R41",
  subgroup = "type2",
  subgroup_remove = "",
  exclusion_condition = NULL,
  subgroup_levels = c("Non severe pneumonia", "Severe pneumonia", "ARDS"),
  extra_columns = c( "Oxygen", "ICU", "Steroid","Daily Dose (mg)", "Duration (days)"),
  rob_column = "rob2_patterns"
)

###############################################################################
# ICU analyses
###############################################################################

# Primary ICU patients only
ICU_primary <- primary[!is.na(primary$ICU) & primary$ICU == "●", ]
primary_ICU_analysis <- analysis_meta(
  outcomes = "Short-term mortality",
  data = ICU_primary,
  metabin = TRUE,
  use_hartung_knapp = TRUE,
  iqwig = TRUE,
  sm = "RR",
  event.e = "event_steroid",
  n.e     = "total_steroid",
  event.c = "event_placebo",
  n.c     = "total_placebo",
  studlab = "Authors+O44LA1:R41",
  subgroup = "type2",
  subgroup_remove = "",
  exclusion_condition = NULL,
  subgroup_levels = c("Severe pneumonia", "ARDS"),
  extra_columns = c( "Oxygen", "ICU", "Steroid","Daily Dose (mg)", "Duration (days)"),
  rob_column = "rob2_patterns"
)

# Secondary ICU patients only
ICU_secondary <- secondary[!is.na(secondary$ICU) & secondary$ICU == "●", ]
secondary_ICU_analysis <- analysis_meta(
  outcomes = "Short-term mortality",
  data = ICU_secondary,
  metabin = TRUE,
  use_hartung_knapp = TRUE,
  iqwig = TRUE,
  sm = "RR",
  event.e = "event_steroid",
  n.e     = "total_steroid",
  event.c = "event_placebo",
  n.c     = "total_placebo",
  studlab = "Authors+O44LA1:R41",
  subgroup = "type2",
  subgroup_remove = "",
  exclusion_condition = NULL,
  subgroup_levels = c("Severe pneumonia", "ARDS"),
  extra_columns = c( "Oxygen", "ICU", "Steroid","Daily Dose (mg)", "Duration (days)"),
  rob_column = "rob2_patterns"
)

###############################################################################
# All pneumonia analyses (excluding ARDS)
###############################################################################

# All pneumonia patients
all_patient <- dt[dt$type2 %in% c("Severe pneumonia","Non severe pneumonia") &!is.na(dt$event_placebo)  & dt$duplicata == 0,]
all_patient_primary <- all_patient[all_patient$CI == 1 | all_patient$type2 == "Non severe pneumonia", ]

all_pneumonia <- analysis_meta(
  outcomes = "Short-term mortality",
  data = all_patient ,
  metabin = TRUE,
  use_hartung_knapp = TRUE,
  iqwig = TRUE,
  sm = "RR",
  event.e = "event_steroid",
  n.e     = "total_steroid",
  event.c = "event_placebo",
  n.c     = "total_placebo",
  studlab = "Authors+O44LA1:R41",
  subgroup = "type2",
  subgroup_remove = "",
  exclusion_condition = NULL,
  subgroup_levels = c("Severe pneumonia", "Non severe pneumonia"),
  extra_columns = c( "Oxygen", "ICU", "Steroid","Daily Dose (mg)", "Duration (days)"),        
  rob_column = "rob2_patterns"       
)

primary_all_pneumonia <- analysis_meta(
  outcomes = "Short-term mortality",
  data = all_patient_primary ,
  metabin = TRUE,
  use_hartung_knapp = TRUE,
  iqwig = TRUE,
  sm = "RR",
  event.e = "event_steroid",
  n.e     = "total_steroid",
  event.c = "event_placebo",
  n.c     = "total_placebo",
  studlab = "Authors+O44LA1:R41",
  subgroup = "type2",
  subgroup_remove = "",
  exclusion_condition = NULL,
  subgroup_levels = c("Severe pneumonia", "Non severe pneumonia"),
  extra_columns = c( "Oxygen", "ICU", "Steroid","Daily Dose (mg)", "Duration (days)"),        
  rob_column = "rob2_patterns"       
)

secondary_all_pneumonia <- analysis_meta(
  outcomes = "Short-term mortality",
  data = secondary[!secondary$type2 == "ARDS", ],
  metabin = TRUE,
  use_hartung_knapp = TRUE,
  iqwig = TRUE,
  sm = "RR",
  event.e = "event_steroid",
  n.e     = "total_steroid",
  event.c = "event_placebo",
  n.c     = "total_placebo",
  studlab = "Authors+O44LA1:R41",
  subgroup = "type2",
  subgroup_remove = "",
  exclusion_condition = NULL,
  subgroup_levels = c("Severe pneumonia", "Non severe pneumonia"),
  extra_columns = c( "Oxygen", "ICU", "Steroid","Daily Dose (mg)", "Duration (days)"),        
  rob_column = "rob2_patterns"       
)

###############################################################################
# Duration analyses - SECONDARY
###############################################################################

# Duration ≤7 days
inf_seven_secondary <- secondary[secondary$duration_max <= 7, ]
secondary_analysis_le7days <- analysis_meta(
  outcomes = "Short-term mortality",
  data = inf_seven_secondary,
  metabin = TRUE,
  use_hartung_knapp = T,
  iqwig = T,
  sm = "RR",
  event.e = "event_steroid",
  n.e     = "total_steroid",
  event.c = "event_placebo",
  n.c     = "total_placebo",
  studlab = "Authors+O44LA1:R41",
  subgroup = "type2",
  subgroup_remove = "",
  exclusion_condition = NULL,
  subgroup_levels = c("Non severe pneumonia", "Severe pneumonia", "ARDS"),
  extra_columns = c( "Oxygen", "ICU", "Steroid","Daily Dose (mg)", "Duration (days)"),
  rob_column = "rob2_patterns"
)

# Duration >7 days
sup_seven_secondary <- secondary[secondary$duration_max > 7, ]
secondary_analysis_gt7days <- analysis_meta(
  outcomes = "Short-term mortality",
  data = sup_seven_secondary,
  metabin = TRUE,
  use_hartung_knapp = T,
  iqwig = T,
  sm = "RR",
  event.e = "event_steroid",
  n.e     = "total_steroid",
  event.c = "event_placebo",
  n.c     = "total_placebo",
  studlab = "Authors+O44LA1:R41",
  subgroup = "type2",
  subgroup_remove = "",
  exclusion_condition = NULL,
  subgroup_levels = c("Non severe pneumonia", "Severe pneumonia", "ARDS"),
  extra_columns = c( "Oxygen", "ICU", "Steroid","Daily Dose (mg)", "Duration (days)"),
  rob_column = "rob2_patterns"
)

###############################################################################
# Dose analyses - SECONDARY
###############################################################################

# Low dose ≤50mg
dose_basse_secondary <- secondary[secondary$dose_moyenne <= 50, ]
secondary_analysis_le50mg <- analysis_meta(
  outcomes = "Short-term mortality",
  data = dose_basse_secondary,
  metabin = TRUE,
  use_hartung_knapp = T,
  iqwig = T,
  sm = "RR",
  event.e = "event_steroid",
  n.e     = "total_steroid",
  event.c = "event_placebo",
  n.c     = "total_placebo",
  studlab = "Authors+O44LA1:R41",
  subgroup = "type2",
  subgroup_remove = "",
  exclusion_condition = NULL,
  subgroup_levels = c("Non severe pneumonia", "Severe pneumonia", "ARDS"),
  extra_columns = c( "Oxygen", "ICU", "Steroid","Daily Dose (mg)", "Duration (days)"),
  rob_column = "rob2_patterns"
)

# High dose >50mg
dose_haute_secondary <- secondary[secondary$dose_moyenne > 50, ]
secondary_analysis_gt50mg <- analysis_meta(
  outcomes = "Short-term mortality",
  data = dose_haute_secondary,
  metabin = TRUE,
  use_hartung_knapp = T,
  iqwig = T,
  sm = "RR",
  event.e = "event_steroid",
  n.e     = "total_steroid",
  event.c = "event_placebo",
  n.c     = "total_placebo",
  studlab = "Authors+O44LA1:R41",
  subgroup = "type2",
  subgroup_remove = "",
  exclusion_condition = NULL,
  subgroup_levels = c("Non severe pneumonia", "Severe pneumonia", "ARDS"),
  extra_columns = c( "Oxygen", "ICU", "Steroid","Daily Dose (mg)", "Duration (days)"),
  rob_column = "rob2_patterns"
)


######################################################################################################################################################################################################################
######################################################################################################################################################################################################################
######################################################################################################################################################################################################################
######################################################################################################################################################################################################################
######################################################################################################################################################################################################################
######################################################################################################################################################################################################################
###############################################################################
###############################################################################
####                    Others outcomes  AND  continuous outcomed                   ####
###############################################################################
###############################################################################
######################################################################################################################################################################################################################
######################################################################################################################################################################################################################
######################################################################################################################################################################################################################
######################################################################################################################################################################################################################
######################################################################################################################################################################################################################




    

# Transform data frame to data.table
setDT(dt)              

###############################################################################
# Function to parse median/IQR data
###############################################################################

parse_stats <- function(x){
  # x: character vector "val1 (val2)" or "val1 (val2–val3)"
  res <- lapply(x, function(s){
    if (is.na(s) || trimws(s) == "") return(list(mean = NA_real_, sd = NA_real_))
    
    # extract all numbers (decimal, comma or dot)
    nums <- str_extract_all(s, "[0-9]+[\\.,]?[0-9]*")[[1]]
    nums <- as.numeric(str_replace(nums, ",", "."))   # comma → dot
    if (length(nums) < 2)    return(list(mean = NA_real_, sd = NA_real_))
    
    if (length(nums) == 2){                 # case mean (SD)
      list(mean = nums[1],
           sd   = nums[2])
    } else {                                # case median (Q1–Q3)
      list(mean = nums[1],                  # mean ≈ median
           sd   = round((nums[3] - nums[2]) / 1.35, 3))
    }
  })
  rbindlist(res)                            # return data.table with 2 columns
}

 # Apply to placebo group
 tmp <- parse_stats(dtmedian_placebo(IQR)`)
 dt[, `:=`(mean_placebo_est = tmp$mean,
          sd_placebo_est   = tmp$sd )]

 # Apply to steroid group
 tmp <- parse_stats(dtmedian_steroid(IQR)`)
 dt[, `:=`(mean_steroid_est = tmp$mean,
          sd_steroid_est   = tmp$sd )]

 # Prepare continuous data
 continu <- dt[, .(
  outcomes,
  total_placebo,
  total_steroid,
  `median_placebo(IQR)`,
  mean_placebo_est,
  sd_placebo_est,
  `median_steroid(IQR)`,
  mean_steroid_est,
  sd_steroid_est,
  `Authors+O44LA1:R41`,
  type2, 
  CI, 
  Oxygen, 
  ICU, 
  Steroid,
  `Daily Dose (mg)`, 
  `Duration (days)`
 )]

 continu <- continu[
  !is.na(continumedian_placebo(IQR)`) & continumedian_placebo(IQR)` != "median_placebo(IQR)",
]

par(ask = FALSE)

###############################################################################
# Continuous outcomes - PRIMARY : Paule-Mandel + Hartung-Knapp
###############################################################################

# Time to clinical cure
ttc_primary <- analysis_meta3(
  outcomes = "Time to clinical cure",
  data = continu[continu$outcomes == "Time to clinical cure" & continu$CI == 1, ],
  cont = TRUE,
  iqwig = TRUE,
  use_hartung_knapp = TRUE,
  n.e = "total_steroid", n.c = "total_placebo",
  mean.e = "mean_steroid_est", sd.e = "sd_steroid_est",
  mean.c = "mean_placebo_est", sd.c = "sd_placebo_est",
  studlab = "Authors+O44LA1:R41",
  subgroup = "type2",
  subgroup_levels = c("Severe pneumonia", "Non severe pneumonia", "ARDS"),
  extra_columns = c("Oxygen", "ICU", "Steroid","Daily Dose (mg)", "Duration (days)"),
  rob_column = "rob2_patterns"
)

# Duration of mechanical ventilation
dmv_primary <- analysis_meta3(
  outcomes = "Duration of mechanical ventilation",
  data = continu[continu$outcomes == "Duration of mechanical ventilation" & continu$CI == 1, ],
  cont = TRUE,
  iqwig = TRUE,
  use_hartung_knapp = TRUE,
  n.e = "total_steroid", n.c = "total_placebo",
  mean.e = "mean_steroid_est", sd.e = "sd_steroid_est",
  mean.c = "mean_placebo_est", sd.c = "sd_placebo_est",
  studlab = "Authors+O44LA1:R41",
  subgroup = "type2",
  subgroup_levels = c("Severe pneumonia", "Non severe pneumonia", "ARDS"),
  extra_columns = c("Oxygen", "ICU", "Steroid","Daily Dose (mg)", "Duration (days)"),
  rob_column = "rob2_patterns"
)

# Length of stay in ICU
los_icu_primary <- analysis_meta3(
  outcomes = "Length of stay in the intensive care unit",
  data = continu[continu$outcomes == "Length of stay in the intensive care unit" & continu$CI == 1, ],
  cont = TRUE,
  iqwig = TRUE,
  use_hartung_knapp = TRUE,
  n.e = "total_steroid", n.c = "total_placebo",
  mean.e = "mean_steroid_est", sd.e = "sd_steroid_est",
  mean.c = "mean_placebo_est", sd.c = "sd_placebo_est",
  studlab = "Authors+O44LA1:R41",
  subgroup = "type2",
  subgroup_levels = c("Severe pneumonia", "Non severe pneumonia", "ARDS"),
  extra_columns = c("Oxygen", "ICU", "Steroid","Daily Dose (mg)", "Duration (days)"),
  rob_column = "rob2_patterns"
)

# Length of stay in hospital
los_hosp_primary <- analysis_meta3(
  outcomes = "Length of stay in the hospital",
  data = continu[continu$outcomes == "Length of stay in the hospital" & continu$CI == 1, ],
  cont = TRUE,
  iqwig = TRUE,
  use_hartung_knapp = TRUE,
  n.e = "total_steroid", n.c = "total_placebo",
  mean.e = "mean_steroid_est", sd.e = "sd_steroid_est",
  mean.c = "mean_placebo_est", sd.c = "sd_placebo_est",
  studlab = "Authors+O44LA1:R41",
  subgroup = "type2",
  subgroup_levels = c("Severe pneumonia", "Non severe pneumonia", "ARDS"),
  extra_columns = c("Oxygen", "ICU", "Steroid","Daily Dose (mg)", "Duration (days)"),
  rob_column = "rob2_patterns"
)

###############################################################################
# Continuous outcomes - SECONDARY
###############################################################################

# Time to clinical cure
ttc_secondary <- analysis_meta3(
  outcomes = "Time to clinical cure",
  data = continu[continu$outcomes == "Time to clinical cure", ],
  cont = TRUE,
  iqwig = TRUE,
  use_hartung_knapp = TRUE,
  n.e = "total_steroid", n.c = "total_placebo",
  mean.e = "mean_steroid_est", sd.e = "sd_steroid_est",
  mean.c = "mean_placebo_est", sd.c = "sd_placebo_est",
  studlab = "Authors+O44LA1:R41",
  subgroup = "type2",
  subgroup_levels = c("Severe pneumonia", "Non severe pneumonia", "ARDS"),
  extra_columns = c("Oxygen", "ICU", "Steroid","Daily Dose (mg)", "Duration (days)"),
  rob_column = "rob2_patterns"
)

# Duration of mechanical ventilation
dmv_secondary <- analysis_meta3(
  outcomes = "Duration of mechanical ventilation",
  data = continu[continu$outcomes == "Duration of mechanical ventilation", ],
  cont = TRUE,
  iqwig = TRUE,
  use_hartung_knapp = TRUE,
  n.e = "total_steroid", n.c = "total_placebo",
  mean.e = "mean_steroid_est", sd.e = "sd_steroid_est",
  mean.c = "mean_placebo_est", sd.c = "sd_placebo_est",
  studlab = "Authors+O44LA1:R41",
  subgroup = "type2",
  subgroup_levels = c("Severe pneumonia", "Non severe pneumonia", "ARDS"),
  extra_columns = c("Oxygen", "ICU", "Steroid","Daily Dose (mg)", "Duration (days)"),
  rob_column = "rob2_patterns"
)

# Length of stay in ICU
los_icu_secondary <- analysis_meta3(
  outcomes = "Length of stay in the intensive care unit",
  data = continu[continu$outcomes == "Length of stay in the intensive care unit", ],
  cont = TRUE,
  iqwig = TRUE,
  use_hartung_knapp = TRUE,
  n.e = "total_steroid", n.c = "total_placebo",
  mean.e = "mean_steroid_est", sd.e = "sd_steroid_est",
  mean.c = "mean_placebo_est", sd.c = "sd_placebo_est",
  studlab = "Authors+O44LA1:R41",
  subgroup = "type2",
  subgroup_levels = c("Severe pneumonia", "Non severe pneumonia", "ARDS"),
  extra_columns = c("Oxygen", "ICU", "Steroid","Daily Dose (mg)", "Duration (days)"),
  rob_column = "rob2_patterns"
)

# Length of stay in hospital
los_hosp_secondary <- analysis_meta3(
  outcomes = "Length of stay in the hospital",
  data = continu[continu$outcomes == "Length of stay in the hospital", ],
  cont = TRUE,
  iqwig = TRUE,
  use_hartung_knapp = TRUE,
  n.e = "total_steroid", n.c = "total_placebo",
  mean.e = "mean_steroid_est", sd.e = "sd_steroid_est",
  mean.c = "mean_placebo_est", sd.c = "sd_placebo_est",
  studlab = "Authors+O44LA1:R41",
  subgroup = "type2",
  subgroup_levels = c("Severe pneumonia", "Non severe pneumonia", "ARDS"),
  extra_columns = c("Oxygen", "ICU", "Steroid","Daily Dose (mg)", "Duration (days)"),
  rob_column = "rob2_patterns"
)

###############################################################################
# Binary outcomes - PRIMARY
###############################################################################

# Need for invasive ventilation
invasive_vent_primary <- analysis_meta(
  outcomes = "Need for invasive ventilation",
  data = dt[dt$outcomes == "Need for invasive ventilation" & !is.na(dt$event_placebo) & dt$CI == 1, ],
  metabin = TRUE, 
  use_hartung_knapp = TRUE, 
  iqwig = TRUE, 
  sm = "RR",
  event.e = "event_steroid", n.e = "total_steroid",
  event.c = "event_placebo", n.c = "total_placebo",
  studlab = "Authors+O44LA1:R41",
  subgroup = "type2", subgroup_remove = "", exclusion_condition = NULL,
  subgroup_levels = c("Severe pneumonia", "Non severe pneumonia", "ARDS"),
  extra_columns = c("Oxygen", "ICU", "Steroid","Daily Dose (mg)", "Duration (days)"),
  rob_column = "rob2_patterns"
)

# Development of respiratory failure
resp_failure_primary <- analysis_meta(
  outcomes = "Development of respiratory failure",
  data = dt[dt$outcomes == "Development of respiratory failure" & !is.na(dt$event_placebo) & dt$CI == 1, ],
  metabin = TRUE, 
  use_hartung_knapp = TRUE, 
  iqwig = TRUE, 
  sm = "RR",
  event.e = "event_steroid", n.e = "total_steroid",
  event.c = "event_placebo", n.c = "total_placebo",
  studlab = "Authors+O44LA1:R41", subgroup = "type2", subgroup_remove = "",
  exclusion_condition = NULL, 
  subgroup_levels = c("Severe pneumonia", "Non severe pneumonia", "ARDS"),
  extra_columns = c("Oxygen", "ICU", "Steroid","Daily Dose (mg)", "Duration (days)"),
  rob_column = "rob2_patterns"
)

# Transfer to ICU
transfer_icu_primary <- analysis_meta(
  outcomes = "Transfer to ICU",
  data = dt[dt$outcomes == "Transfer to ICU" & !is.na(dt$event_placebo) & dt$CI == 1, ],
  metabin = TRUE, 
  use_hartung_knapp = TRUE, 
  iqwig = TRUE, 
  sm = "RR",
  event.e = "event_steroid", n.e = "total_steroid",
  event.c = "event_placebo", n.c = "total_placebo",
  studlab = "Authors+O44LA1:R41", subgroup = "type2", subgroup_remove = "",
  exclusion_condition = NULL, 
  subgroup_levels = c("Severe pneumonia", "Non severe pneumonia", "ARDS"),
  extra_columns = c("Oxygen", "ICU", "Steroid","Daily Dose (mg)", "Duration (days)"),
  rob_column = "rob2_patterns"
)

# ICU-acquired neuromyopathy
neuromyopathy_primary <- analysis_meta(
  outcomes = "Intensive care unit-acquired neuromyopathy",
  data = dt[dt$outcomes == "Intensive care unit-acquired neuromyopathy" & !is.na(dt$event_placebo) & dt$CI == 1, ],
  metabin = TRUE, 
  use_hartung_knapp = TRUE, 
  iqwig = TRUE, 
  sm = "RR",
  event.e = "event_steroid", n.e = "total_steroid",
  event.c = "event_placebo", n.c = "total_placebo",
  studlab = "Authors+O44LA1:R41", subgroup = "type2", subgroup_remove = "",
  exclusion_condition = NULL, 
  subgroup_levels = c("Severe pneumonia", "Non severe pneumonia", "ARDS"),
  extra_columns = c("Oxygen", "ICU", "Steroid","Daily Dose (mg)", "Duration (days)"),
  rob_column = "rob2_patterns"
)

# Gastrointestinal bleeding
gi_bleeding_primary <- analysis_meta(
  outcomes = "Gastrointestinal bleeding",
  data = dt[dt$outcomes == "Gastrointestinal bleeding" & !is.na(dt$event_placebo) & dt$CI == 1, ],
  metabin = TRUE, 
  use_hartung_knapp = TRUE, 
  iqwig = TRUE, 
  sm = "RR",
  event.e = "event_steroid", n.e = "total_steroid",
  event.c = "event_placebo", n.c = "total_placebo",
  studlab = "Authors+O44LA1:R41", subgroup = "type2", subgroup_remove = "",
  exclusion_condition = NULL, 
  subgroup_levels = c("Severe pneumonia", "Non severe pneumonia", "ARDS"),
  extra_columns = c("Oxygen", "ICU", "Steroid","Daily Dose (mg)", "Duration (days)"),
  rob_column = "rob2_patterns"
)

# Hyperglycemia
hyperglycemia_primary <- analysis_meta(
  outcomes = "Hyperglycemia",
  data = dt[dt$outcomes == "Hyperglycemia" & !is.na(dt$event_placebo) & dt$CI == 1, ],
  metabin = TRUE, 
  use_hartung_knapp = TRUE, 
  iqwig = TRUE, 
  sm = "RR",
  event.e = "event_steroid", n.e = "total_steroid",
  event.c = "event_placebo", n.c = "total_placebo",
  studlab = "Authors+O44LA1:R41", subgroup = "type2", subgroup_remove = "",
  exclusion_condition = NULL, 
  subgroup_levels = c("Severe pneumonia", "Non severe pneumonia", "ARDS"),
  extra_columns = c("Oxygen", "ICU", "Steroid","Daily Dose (mg)", "Duration (days)"),
  rob_column = "rob2_patterns"
)

# Cardiac complications
cardiac_comp_primary <- analysis_meta(
  outcomes = "Cardiac complications",
  data = dt[dt$outcomes == "Cardiac complications" & !is.na(dt$event_placebo) & dt$CI == 1, ],
  metabin = TRUE, 
  use_hartung_knapp = TRUE, 
  iqwig = TRUE, 
  sm = "RR",
  event.e = "event_steroid", n.e = "total_steroid",
  event.c = "event_placebo", n.c = "total_placebo",
  studlab = "Authors+O44LA1:R41", subgroup = "type2", subgroup_remove = "",
  exclusion_condition = NULL, 
  subgroup_levels = c("Severe pneumonia", "Non severe pneumonia", "ARDS"),
  extra_columns = c("Oxygen", "ICU", "Steroid","Daily Dose (mg)", "Duration (days)"),
  rob_column = "rob2_patterns"
)

# Long-term mortality
longterm_mort_primary <- analysis_meta(
  outcomes = "Long-term mortality",
  data = dt[dt$outcomes == "Long-term mortality" & !is.na(dt$event_placebo) & dt$CI == 1, ],
  metabin = TRUE, 
  use_hartung_knapp = TRUE, 
  iqwig = TRUE, 
  sm = "RR",
  event.e = "event_steroid", n.e = "total_steroid",
  event.c = "event_placebo", n.c = "total_placebo",
  studlab = "Authors+O44LA1:R41", subgroup = "type2", subgroup_remove = "",
  exclusion_condition = NULL, 
  subgroup_levels = c("Severe pneumonia", "Non severe pneumonia", "ARDS"),
  extra_columns = c("Oxygen", "ICU", "Steroid","Daily Dose (mg)", "Duration (days)"),
  rob_column = "rob2_patterns"
)

###############################################################################
# Binary outcomes - SECONDARY
###############################################################################

# Need for invasive ventilation
invasive_vent_secondary <- analysis_meta(
  outcomes = "Need for invasive ventilation",
  data = dt[dt$outcomes == "Need for invasive ventilation" & !is.na(dt$event_placebo), ],
  metabin = TRUE, 
  use_hartung_knapp = TRUE, 
  iqwig = TRUE, 
  sm = "RR",
  event.e = "event_steroid", n.e = "total_steroid",
  event.c = "event_placebo", n.c = "total_placebo",
  studlab = "Authors+O44LA1:R41",
  subgroup = "type2", subgroup_remove = "", exclusion_condition = NULL,
  subgroup_levels = c("Severe pneumonia", "Non severe pneumonia", "ARDS"),
  extra_columns = c("Oxygen", "ICU", "Steroid","Daily Dose (mg)", "Duration (days)"),
  rob_column = "rob2_patterns"
)

# Development of respiratory failure
resp_failure_secondary <- analysis_meta(
  outcomes = "Development of respiratory failure",
  data = dt[dt$outcomes == "Development of respiratory failure" & !is.na(dt$event_placebo), ],
  metabin = TRUE, 
  use_hartung_knapp = TRUE, 
  iqwig = TRUE, 
  sm = "RR",
  event.e = "event_steroid", n.e = "total_steroid",
  event.c = "event_placebo", n.c = "total_placebo",
  studlab = "Authors+O44LA1:R41", subgroup = "type2", subgroup_remove = "",
  exclusion_condition = NULL, 
  subgroup_levels = c("Severe pneumonia", "Non severe pneumonia", "ARDS"),
  extra_columns = c("Oxygen", "ICU", "Steroid","Daily Dose (mg)", "Duration (days)"),
  rob_column = "rob2_patterns"
)

# Transfer to ICU
transfer_icu_secondary <- analysis_meta(
  outcomes = "Transfer to ICU",
  data = dt[dt$outcomes == "Transfer to ICU" & !is.na(dt$event_placebo), ],
  metabin = TRUE, 
  use_hartung_knapp = TRUE, 
  iqwig = TRUE, 
  sm = "RR",
  event.e = "event_steroid", n.e = "total_steroid",
  event.c = "event_placebo", n.c = "total_placebo",
  studlab = "Authors+O44LA1:R41", subgroup = "type2", subgroup_remove = "",
  exclusion_condition = NULL, 
  subgroup_levels = c("Severe pneumonia", "Non severe pneumonia", "ARDS"),
  extra_columns = c("Oxygen", "ICU", "Steroid","Daily Dose (mg)", "Duration (days)"),
  rob_column = "rob2_patterns"
)

# ICU-acquired neuromyopathy
neuromyopathy_secondary <- analysis_meta(
  outcomes = "Intensive care unit-acquired neuromyopathy",
  data = dt[dt$outcomes == "Intensive care unit-acquired neuromyopathy" & !is.na(dt$event_placebo), ],
  metabin = TRUE, 
  use_hartung_knapp = TRUE, 
  iqwig = TRUE, 
  sm = "RR",
  event.e = "event_steroid", n.e = "total_steroid",
  event.c = "event_placebo", n.c = "total_placebo",
  studlab = "Authors+O44LA1:R41", subgroup = "type2", subgroup_remove = "",
  exclusion_condition = NULL, 
  subgroup_levels = c("Severe pneumonia", "Non severe pneumonia", "ARDS"),
  extra_columns = c("Oxygen", "ICU", "Steroid","Daily Dose (mg)", "Duration (days)"),
  rob_column = "rob2_patterns"
)

# Gastrointestinal bleeding
gi_bleeding_secondary <- analysis_meta(
  outcomes = "Gastrointestinal bleeding",
  data = dt[dt$outcomes == "Gastrointestinal bleeding" & !is.na(dt$event_placebo), ],
  metabin = TRUE, 
  use_hartung_knapp = TRUE, 
  iqwig = TRUE, 
  sm = "RR",
  event.e = "event_steroid", n.e = "total_steroid",
  event.c = "event_placebo", n.c = "total_placebo",
  studlab = "Authors+O44LA1:R41", subgroup = "type2", subgroup_remove = "",
  exclusion_condition = NULL, 
  subgroup_levels = c("Severe pneumonia", "Non severe pneumonia", "ARDS"),
  extra_columns = c("Oxygen", "ICU", "Steroid","Daily Dose (mg)", "Duration (days)"),
  rob_column = "rob2_patterns"
)

# Hyperglycemia
hyperglycemia_secondary <- analysis_meta(
  outcomes = "Hyperglycemia",
  data = dt[dt$outcomes == "Hyperglycemia" & !is.na(dt$event_placebo), ],
  metabin = TRUE, 
  use_hartung_knapp = TRUE, 
  iqwig = TRUE, 
  sm = "RR",
  event.e = "event_steroid", n.e = "total_steroid",
  event.c = "event_placebo", n.c = "total_placebo",
  studlab = "Authors+O44LA1:R41", subgroup = "type2", subgroup_remove = "",
  exclusion_condition = NULL, 
  subgroup_levels = c("Severe pneumonia", "Non severe pneumonia", "ARDS"),
  extra_columns = c("Oxygen", "ICU", "Steroid","Daily Dose (mg)", "Duration (days)"),
  rob_column = "rob2_patterns"
)

# Cardiac complications
cardiac_comp_secondary <- analysis_meta(
  outcomes = "Cardiac complications",
  data = dt[dt$outcomes == "Cardiac complications" & !is.na(dt$event_placebo), ],
  metabin = TRUE, 
  use_hartung_knapp = TRUE, 
  iqwig = TRUE, 
  sm = "RR",
  event.e = "event_steroid", n.e = "total_steroid",
  event.c = "event_placebo", n.c = "total_placebo",
  studlab = "Authors+O44LA1:R41", subgroup = "type2", subgroup_remove = "",
  exclusion_condition = NULL, 
  subgroup_levels = c("Severe pneumonia", "Non severe pneumonia", "ARDS"),
  extra_columns = c("Oxygen", "ICU", "Steroid","Daily Dose (mg)", "Duration (days)"),
  rob_column = "rob2_patterns"
)

# Long-term mortality
longterm_mort_secondary <- analysis_meta(
  outcomes = "Long-term mortality",
  data = dt[dt$outcomes == "Long-term mortality" & !is.na(dt$event_placebo), ],
  metabin = TRUE, 
  use_hartung_knapp = TRUE, 
  iqwig = TRUE, 
  sm = "RR",
  event.e = "event_steroid", n.e = "total_steroid",
  event.c = "event_placebo", n.c = "total_placebo",
  studlab = "Authors+O44LA1:R41", subgroup = "type2", subgroup_remove = "",
  exclusion_condition = NULL, 
  subgroup_levels = c("Severe pneumonia", "Non severe pneumonia", "ARDS"),
  extra_columns = c("Oxygen", "ICU", "Steroid","Daily Dose (mg)", "Duration (days)"),
  rob_column = "rob2_patterns"
)




######################################################################################################################################################################################################################
######################################################################################################################################################################################################################
######################################################################################################################################################################################################################
######################################################################################################################################################################################################################
######################################################################################################################################################################################################################
######################################################################################################################################################################################################################
###############################################################################
###############################################################################
####                    FUNNEL PLOT                  ####
###############################################################################
###############################################################################
######################################################################################################################################################################################################################
######################################################################################################################################################################################################################
######################################################################################################################################################################################################################
######################################################################################################################################################################################################################
######################################################################################################################################################################################################################


# Function
funnel_fonction_colored_final <- function(meta_object, outcome_name, title = "Funnel plot", 
                                          xlim = NULL, ylim = NULL, show_labels = TRUE) {
  # Étape 1 : Extraire les noms des études dans l'objet méta
  study_labels <- meta_object$studlab

  # Étape 2 : Associer le sous-groupe depuis dt
  subgroup_vector <- sapply(study_labels, function(name) {
    match_row <- dt[dt$`Authors+O44LA1:R41` == name & dt$outcomes == outcome_name, ]
    if (nrow(match_row) > 0) {
      return(as.character(match_row$type2[1]))
    } else {
      return(NA)
    }
  })

  # Étape 3 : Couleurs JAMA
  jama_colors <- c("ARDS" = "#0072B2", "Severe pneumonia" = "#D55E00", "Non severe pneumonia" = "#009E73")
  point_colors <- jama_colors[subgroup_vector]

  # Étape 4 : Affichage du funnel
  funnel(meta_object,
         bg = point_colors,
         pch = 21,
         studlab = show_labels,
         contour = c(0.9, 0.95, 0.99),
         col.contour = c("lightsteelblue", "powderblue", "aliceblue"),
         level = 0.95,
         xlim = xlim,
         ylim = ylim)

  # Étape 5 : Légende
  legend("topright",
         legend = names(jama_colors),
         fill = jama_colors,
         border = "black",
         title = "Subgroup",
         cex = 0.8)

  # Étape 6 : Titre
  title(main = title)
}



# Primary analsyis : low dose, short duration



# Funnel plots – Infections
funnel_fonction_colored_final(
  meta_object = infection_ha$meta,
  outcome_name = "Hospital-acquired infections",
  title = "Funnel Plot – Hospital-Acquired Infections",
  show_labels = FALSE
)

funnel_fonction_colored_final(
  meta_object = infection_ss$meta,
  outcome_name = "Occurrence of septic shock",
  title = "Funnel Plot – Shock",
  show_labels = FALSE
)

funnel_fonction_colored_final(
  meta_object = infection_sp$meta,
  outcome_name = "Secondary pneumonia",
  title = "Funnel Plot – Secondary Pneumonia",
  show_labels = FALSE
)

funnel_fonction_colored_final(
  meta_object = infection_cat$meta,
  outcome_name = "Catheter-related infection",
  title = "Funnel Plot – Catheter-Related Infection",
  show_labels = FALSE
)

funnel_fonction_colored_final(
  meta_object = infection_bld$meta,
  outcome_name = "Bacteremia",
  title = "Funnel Plot – Bacteremia",
  show_labels = FALSE
)

# Funnel plot – Time to clinical cure (continuous)
funnel_fonction_colored_final(
  meta_object  = ttc_primary$meta,
  outcome_name = "Time to clinical cure",
  title        = "Funnel Plot – Time to Clinical Cure",
  show_labels  = FALSE
)

# Funnel plot – Duration of mechanical ventilation (continuous)
funnel_fonction_colored_final(
  meta_object  = dmv_primary$meta,
  outcome_name = "Duration of mechanical ventilation",
  title        = "Funnel Plot – Duration of Mechanical Ventilation",
  show_labels  = FALSE
)

# Funnel plot – ICU length of stay (continuous)
funnel_fonction_colored_final(
  meta_object  = los_icu_primary$meta,
  outcome_name = "ICU length of stay",
  title        = "Funnel Plot – Length of Stay in ICU",
  show_labels  = FALSE
)

# Funnel plot – Hospital length of stay (continuous)
funnel_fonction_colored_final(
  meta_object  = los_hosp_primary$meta,
  outcome_name = "Hospital length of stay",
  title        = "Funnel Plot – Length of Stay in Hospital",
  show_labels  = FALSE
)

# Funnel plot – Need for invasive ventilation (binary)
funnel_fonction_colored_final(
  meta_object  = invasive_vent_primary$meta,
  outcome_name = "Need for invasive ventilation",
  title        = "Funnel Plot – Need for Invasive Ventilation",
  show_labels  = FALSE
)

# Funnel plot – Development of respiratory failure (binary)
funnel_fonction_colored_final(
  meta_object  = resp_failure_primary$meta,
  outcome_name = "Development of respiratory failure",
  title        = "Funnel Plot – Development of Respiratory Failure",
  show_labels  = FALSE
)

# Funnel plot – Transfer to ICU (binary)
funnel_fonction_colored_final(
  meta_object  = transfer_icu_primary$meta,
  outcome_name = "Transfer to ICU",
  title        = "Funnel Plot – Transfer to ICU",
  show_labels  = FALSE
)

# Funnel plot – ICU-acquired neuromyopathy (binary)
funnel_fonction_colored_final(
  meta_object  = neuromyopathy_primary$meta,
  outcome_name = "Intensive care unit-acquired neuromyopathy",
  title        = "Intensive care unit-acquired neuromyopathy",
  show_labels  = FALSE
)

# Funnel plot – Gastrointestinal bleeding (binary)
funnel_fonction_colored_final(
  meta_object  = gi_bleeding_primary$meta,
  outcome_name = "Gastrointestinal bleeding",
  title        = "Funnel Plot – Gastrointestinal Bleeding",
  show_labels  = FALSE
)

# Funnel plot – Hyperglycemia (binary)
funnel_fonction_colored_final(
  meta_object  = hyperglycemia_primary$meta,
  outcome_name = "Hyperglycemia",
  title        = "Funnel Plot – Hyperglycemia",
  show_labels  = FALSE
)

# Funnel plot – Cardiac complications (binary)
funnel_fonction_colored_final(
  meta_object  = cardiac_comp_primary$meta,
  outcome_name = "Cardiac complications",
  title        = "Funnel Plot – Cardiac Complications",
  show_labels  = FALSE
)

# Funnel plot – Long-term mortality (binary)
funnel_fonction_colored_final(
  meta_object  = longterm_mort_primary$meta,
  outcome_name = "Long-term mortality",
  title        = "Funnel Plot – Long-Term Mortality",
  show_labels  = FALSE
)



# Secondary analysis without restriction on dose, duration or timinginginginginginginginginginginginginginginging

# Funnel plots – Infections
funnel_fonction_colored_final(
  meta_object = infection_ha_all$meta,
  outcome_name = "Hospital-acquired infections",
  title = "Funnel Plot – Hospital-Acquired Infections",
  show_labels = FALSE
)


funnel_fonction_colored_final(
  meta_object = infection_ss_all$meta,
  outcome_name = "Occurrence of septic shock",
  title = "Funnel Plot – Occurrence of septic shock",
  show_labels = FALSE
)



funnel_fonction_colored_final(
  meta_object = infection_sp_all$meta,
  outcome_name = "Secondary pneumonia",
  title = "Funnel Plot – Secondary Pneumonia",
  show_labels = FALSE
)

funnel_fonction_colored_final(
  meta_object = infection_cat_all$meta,
  outcome_name = "Catheter-related infection",
  title = "Funnel Plot – Catheter-Related Infection",
  show_labels = FALSE
)

funnel_fonction_colored_final(
  meta_object = infection_bld_all$meta,
  outcome_name = "Bacteremia",
  title = "Funnel Plot – Bacteremia",
  show_labels = FALSE
)

# Funnel plot – Time to clinical cure (continuous)
funnel_fonction_colored_final(
  meta_object  = ttc_secondary$meta,,
  outcome_name = "Time to clinical cure",
  title        = "Funnel Plot – Time to Clinical Cure",
  show_labels  = FALSE
)

# Funnel plot – Duration of mechanical ventilation (continuous)
funnel_fonction_colored_final(
  meta_object  = dmv_secondary$meta,
  outcome_name = "Duration of mechanical ventilation",
  title        = "Funnel Plot – Duration of Mechanical Ventilation",
  show_labels  = FALSE
)



# Funnel plot – ICU length of stay (continuous)
funnel_fonction_colored_final(
  meta_object  = los_icu_secondary$meta,
  outcome_name = "ICU length of stay",
  title        = "Funnel Plot – Length of Stay in ICU",
  show_labels  = FALSE
)

# Funnel plot – Hospital length of stay (continuous)
funnel_fonction_colored_final(
  meta_object  = los_hosp_secondary$meta,
  outcome_name = "Hospital length of stay",
  title        = "Funnel Plot – Length of Stay in Hospital",
  show_labels  = FALSE
)

# Funnel plot – Need for invasive ventilation (binary)
funnel_fonction_colored_final(
  meta_object  = invasive_vent_secondary$meta,
  outcome_name = "Need for invasive ventilation",
  title        = "Funnel Plot – Need for Invasive Ventilation",
  show_labels  = FALSE
)

# Funnel plot – Development of respiratory failure (binary)
funnel_fonction_colored_final(
  meta_object  = resp_failure_secondary$meta,
  outcome_name = "Development of respiratory failure",
  title        = "Funnel Plot – Development of Respiratory Failure",
  show_labels  = FALSE
)

# Funnel plot – Transfer to ICU (binary)
funnel_fonction_colored_final(
  meta_object  = transfer_icu_secondary$meta,
  outcome_name = "Transfer to ICU",
  title        = "Funnel Plot – Transfer to ICU",
  show_labels  = FALSE
)

# Funnel plot – ICU-acquired neuromyopathy (binary)
funnel_fonction_colored_final(
  meta_object  = neuromyopathy_secondary$meta,
  outcome_name = "Intensive care unit-acquired neuromyopathy",
  title        = "Funnel Plot – ICU-Acquired Neuromyopathy",
  show_labels  = FALSE
)

# Funnel plot – Gastrointestinal bleeding (binary)
funnel_fonction_colored_final(
  meta_object  = gi_bleeding_secondary$meta,
  outcome_name = "Gastrointestinal bleeding",
  title        = "Funnel Plot – Gastrointestinal Bleeding",
  show_labels  = FALSE
)

# Funnel plot – Hyperglycemia (binary)
funnel_fonction_colored_final(
  meta_object  = hyperglycemia_secondary$meta,
  outcome_name = "Hyperglycemia",
  title        = "Funnel Plot – Hyperglycemia",
  show_labels  = FALSE
)

# Funnel plot – Cardiac complications (binary)
funnel_fonction_colored_final(
  meta_object  = cardiac_comp_secondary$meta,
  outcome_name = "Cardiac complications",
  title        = "Funnel Plot – Cardiac Complications",
  show_labels  = FALSE
)

# Funnel plot – Long-term mortality (binary)
funnel_fonction_colored_final(
  meta_object  = longterm_mort_secondary$meta,
  outcome_name = "Long-term mortality",
  title        = "Funnel Plot – Long-Term Mortality",
  show_labels  = FALSE
)








