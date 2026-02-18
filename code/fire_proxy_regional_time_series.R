# =============================================================================
# Regional time series of North American fire proxy data
# =============================================================================

# Clear the environment
rm(list = ls(all.names = TRUE))
cat("\014")

library(tidyverse)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggrepel)
library(locfit)

# =============================== Set up paths ================================ ####

base_dir <- file.path('/Users/Nick/Google_Drive/RESEARCH_PROJECTS/Paleo_aerosols/paleo-aerosols')
data_dir <- file.path(base_dir, "data")
output_dir <- file.path(base_dir, "output/working_fig_pieces")

dir.create(output_dir, recursive = TRUE)

# ================================= Switches ================================== ####
printOn <- TRUE
smooth_scar_data <- TRUE
iceMode <- 'raw'  # Options: 'inversion' or 'raw'
charSmoothType <- 'narrow' #Options: 'wide' or 'narrow'
iceSmoothType <- 'all' #Options: 'all' or 'indiv'

text_size_multiplier <- 0.8
fig_w <- 3
fig_h <- 2.5

# Smoothing parameters (example values - adjust as needed)
smooth_params <- list(
  hw1 = 5,             # Half-width in years for narrow smoothing
  hw2 = 15,            # Half-width in years for wide smoothing
  degree = 1,           # Degree of polynomial (1 = local linear)
  family = "qrgauss"   # Family for locfit ("gaussian" is standard)
)

# ========================= Define custom functions =========================== ####

# ---------------------------- smooth_scar_series -------------------------------- #

# Function to smooth a single scar dataset with two smoothing levels
smooth_scar_series <- function(scar_data, hw1, hw2, degree = 1, family = "gaussian") {
  
  # Extract x (year) and y (percent_scarred)
  x <- scar_data$year_ce
  y <- scar_data$percent_scarred
  
  # Remove any NA values
  valid_idx <- !is.na(x) & !is.na(y)
  x_valid <- x[valid_idx]
  y_valid <- y[valid_idx]
  
  if(length(x_valid) < 3) {
    warning("Not enough data points for smoothing")
    return(scar_data)
  }
  
  # Initialize smoothed columns with NA
  scar_data$percent_scarred_smooth1 <- NA
  scar_data$percent_scarred_smooth2 <- NA
  
  # Fit narrow LOWESS curve using only valid data
  loc_fit_narrow <- locfit(y_valid ~ lp(x_valid, deg = degree, h = hw1), family = family)
  
  # Fit wide LOWESS curve using only valid data
  loc_fit_wide <- locfit(y_valid ~ lp(x_valid, deg = degree, h = hw2), family = family)
  
  # Predict only for valid indices
  smoothed_narrow <- predict(loc_fit_narrow, newdata = data.frame(x_valid = x[valid_idx]))
  smoothed_wide <- predict(loc_fit_wide, newdata = data.frame(x_valid = x[valid_idx]))
  
  # Add smoothed values only where we have valid data
  scar_data$percent_scarred_smooth1[valid_idx] <- smoothed_narrow
  scar_data$percent_scarred_smooth2[valid_idx] <- smoothed_wide
  
  return(scar_data)
}

# ---------------------------- zscore_normalize --------------------------------- #

# Z-score normalization function
zscore_normalize <- function(data_list) {
  lapply(data_list, function(df) {
    if(!is.null(df) && "rBC_deposition" %in% names(df)) {
      # Calculate mean and sd for the entire record
      mean_val <- mean(df$rBC_deposition, na.rm = TRUE)
      sd_val <- sd(df$rBC_deposition, na.rm = TRUE)
      
      # Create z-scored version
      df$rBC_deposition_zscore <- (df$rBC_deposition - mean_val) / sd_val
      
      # Replace the original column with z-scored values
      df$rBC_deposition <- df$rBC_deposition_zscore
      df$rBC_deposition_zscore <- NULL  # Remove temporary column
    }
    return(df)
  })
}

# ------------------------- smooth_ice_core_series ------------------------------ #

# Smooth ice core time series function
smooth_ice_core_series <- function(data_list, hw1, hw2, degree = 1, family = "gaussian") {
  lapply(data_list, function(df) {
    if(!is.null(df) && "rBC_deposition" %in% names(df)) {
      
      # Extract x (year) and y (rBC_deposition)
      x <- df$year_ce
      y <- df$rBC_deposition
      
      # Remove any NA values
      valid_idx <- !is.na(x) & !is.na(y)
      x_valid <- x[valid_idx]
      y_valid <- y[valid_idx]
      
      if(length(x_valid) < 3) {
        warning("Not enough data points for smoothing")
        return(df)
      }
      
      # Initialize smoothed columns with NA
      df$rBC_deposition_smooth1 <- NA
      df$rBC_deposition_smooth2 <- NA
      
      # Fit narrow LOWESS curve using only valid data
      loc_fit_narrow <- locfit(y_valid ~ lp(x_valid, deg = degree, h = hw1), family = family)
      
      # Fit wide LOWESS curve using only valid data
      loc_fit_wide <- locfit(y_valid ~ lp(x_valid, deg = degree, h = hw2), family = family)
      
      # Predict only for valid indices
      smoothed_narrow <- predict(loc_fit_narrow, newdata = data.frame(x_valid = x[valid_idx]))
      smoothed_wide <- predict(loc_fit_wide, newdata = data.frame(x_valid = x[valid_idx]))
      
      # Add smoothed values only where we have valid data
      df$rBC_deposition_smooth1[valid_idx] <- smoothed_narrow
      df$rBC_deposition_smooth2[valid_idx] <- smoothed_wide
    }
    return(df)
  })
}

# ------------------------- smooth_ice_core_collection -------------------------- #

# Function to fit a single smooth curve to a collection of z-scored ice core time series
smooth_ice_core_collection <- function(data_list, hw1, hw2, degree = 1, family = "gaussian") {
  
  # Step 1: Pool all data points from all sites
  all_years <- c()
  all_values <- c()
  
  for(site_name in names(data_list)) {
    df <- data_list[[site_name]]
    
    if(!is.null(df) && "rBC_deposition" %in% names(df)) {
      # Extract x (year) and y (rBC_deposition - already z-scored)
      x <- df$year_ce
      y <- df$rBC_deposition
      
      # Remove any NA values
      valid_idx <- !is.na(x) & !is.na(y)
      x_valid <- x[valid_idx]
      y_valid <- y[valid_idx]
      
      # Add to pooled data
      all_years <- c(all_years, x_valid)
      all_values <- c(all_values, y_valid)
    }
  }
  
  if(length(all_years) < 3) {
    warning("Not enough pooled data points for smoothing")
    return(data_list)
  }
  
  cat("Pooled", length(all_years), "data points from", length(data_list), "sites\n")
  
  # Step 2: Fit smoothed curves to the pooled data
  loc_fit_narrow <- locfit(all_values ~ lp(all_years, deg = degree, h = hw1), family = family)
  loc_fit_wide <- locfit(all_values ~ lp(all_years, deg = degree, h = hw2), family = family)
  
  # Step 3: Create a common time grid for predictions
  year_range <- range(all_years)
  common_years <- seq(floor(year_range[1]), ceiling(year_range[2]), by = 1)
  
  # Predict smooth curves on common grid
  smoothed_narrow <- predict(loc_fit_narrow, newdata = data.frame(all_years = common_years))
  smoothed_wide <- predict(loc_fit_wide, newdata = data.frame(all_years = common_years))
  
  # Create a dataframe with the smoothed curves
  smooth_df <- data.frame(
    year_ce = common_years,
    rBC_deposition_smooth1 = smoothed_narrow,
    rBC_deposition_smooth2 = smoothed_wide
  )
  
  # Step 4: Add smoothed values to each site's dataframe (matching by year)
  for(site_name in names(data_list)) {
    df <- data_list[[site_name]]
    
    if(!is.null(df) && "rBC_deposition" %in% names(df)) {
      # Initialize smooth columns
      df$rBC_deposition_smooth1 <- NA
      df$rBC_deposition_smooth2 <- NA
      
      # Match years and add smooth values
      for(i in 1:nrow(df)) {
        year <- df$year_ce[i]
        match_idx <- which(smooth_df$year_ce == year)
        
        if(length(match_idx) > 0) {
          df$rBC_deposition_smooth1[i] <- smooth_df$rBC_deposition_smooth1[match_idx]
          df$rBC_deposition_smooth2[i] <- smooth_df$rBC_deposition_smooth2[match_idx]
        }
      }
      
      data_list[[site_name]] <- df
    }
  }
  
  return(data_list)
}

# ----------------------------- plot_timeseries --------------------------------- #

# plot_timeseries <- function(data_list, pltParams,
#                             age_var = "year",           # Name of age/time column
#                             value_var = "value",        # Name of value column to plot
#                             ci_upper_var = NULL,        # Name of upper CI column (optional)
#                             ci_lower_var = NULL,        # Name of lower CI column (optional)
#                             smooth1_var = NULL,         # Name of first smooth curve column (optional)
#                             smooth2_var = NULL,         # Name of second smooth curve column (optional)
#                             plot_smooth1 = FALSE,       # Whether to plot first smooth curve
#                             plot_smooth2 = FALSE,       # Whether to plot second smooth curve
#                             smooth1_label = "smooth",   # Label for first smooth curve in legend
#                             smooth2_label = "smooth2",  # Label for second smooth curve in legend
#                             ylab = "Fire emissions",    # Y-axis label
#                             save_path = NULL,
#                             fname_prefix = "timeseries",
#                             save_transparent = TRUE,
#                             save_variants = TRUE) {
# 
#   # Set up y-limits and x-limits
#   xlims <- c(pltParams[["xmin"]], pltParams[["xmax"]])
#   ylims <- c(pltParams[["ymin"]], pltParams[["ymax"]])
# 
#   # Create print directory
#   if(is.null(save_path)) {
#     print_dir <- output_dir
#   } else {
#     print_dir <- save_path
#   }
#   dir.create(print_dir, recursive = TRUE, showWarnings = FALSE)
# 
#   # Function to create a single plot with options
#   create_single_plot <- function(show_legend = TRUE,
#                                  show_xaxis = TRUE,
#                                  yaxis_right = FALSE) {
# 
#     # Set margins based on options
#     if (yaxis_right && !show_xaxis) {
#       par(mar = c(1, 1, 1, 4), cex = 0.75 * text_size_multiplier)
#     } else if (!show_xaxis) {
#       par(mar = c(1, 4, 1, 2), cex = 0.75 * text_size_multiplier)
#     } else if (yaxis_right) {
#       par(mar = c(5, 2, 4, 4), cex = 0.75 * text_size_multiplier)
#     } else {
#       par(mar = c(5, 4, 4, 2), cex = 0.75 * text_size_multiplier)
#     }
# 
#     # Determine axis parameters
#     yaxt_param <- if(yaxis_right) "n" else "s"
#     xaxt_param <- if(!show_xaxis) "n" else "s"
#     bty_param <- if(!show_xaxis) "n" else "o"
# 
#     plot(NA, type = "n", xlim = xlims, ylim = ylims,
#          xlab = if(show_xaxis) "Year (CE)" else "",
#          ylab = if(!yaxis_right) ylab else "",
#          main = "",
#          yaxt = yaxt_param,
#          xaxt = xaxt_param,
#          bty = bty_param)
# 
#     # Add y-axis on right if needed
#     if(yaxis_right) {
#       axis(4, las = 1)
#       mtext(ylab, side = 4, line = 3, cex = 0.75 * text_size_multiplier)
#     }
# 
#     # For plots without x-axis, manually draw only the y-axis line
#     if(!show_xaxis) {
#       usr <- par("usr")
# 
#       if(yaxis_right) {
#         lines(c(usr[2], usr[2]), c(usr[3], usr[4]), col = "black", lwd = 1)
#       } else {
#         lines(c(usr[1], usr[1]), c(usr[3], usr[4]), col = "black", lwd = 1)
#       }
#     }
# 
#     # Add major and minor ticks to x-axis if shown
#     if(show_xaxis) {
#       major_ticks <- pretty(xlims)
#       minor_ticks <- unlist(lapply(1:(length(major_ticks) - 1), function(k) {
#         seq(major_ticks[k], major_ticks[k + 1], length.out = 6)[-c(1, 6)]
#       }))
# 
#       axis(side = 1, at = major_ticks, labels = TRUE)
#       axis(side = 1, at = minor_ticks, labels = FALSE, tcl = -0.3)
#       axis(side = 3, at = major_ticks, labels = FALSE, tcl = 0.5)
#       axis(side = 3, at = minor_ticks, labels = FALSE, tcl = 0.25)
#     }
# 
#     # Plot data
#     for(dataset_name in names(data_list)) {
#       data <- data_list[[dataset_name]]
# 
#       # Get colors for raw and smooth curves
#       color_raw <- pltParams[["colors"]]["raw", dataset_name]
#       color_smooth1 <- pltParams[["colors"]]["smooth1", dataset_name]
#       color_smooth2 <- pltParams[["colors"]]["smooth2", dataset_name]
# 
#       # Get alpha values (default to 1.0 if not set)
#       raw_alpha <- if("raw_alpha" %in% names(pltParams)) pltParams[["raw_alpha"]] else 1.0
#       smooth_alpha <- if("smooth_alpha" %in% names(pltParams)) pltParams[["smooth_alpha"]] else 1.0
# 
#       # Extract variables using flexible column names
#       age <- data[[age_var]]
#       value <- data[[value_var]]
# 
#       # Plot confidence interval if available
#       if(!is.null(ci_upper_var) && !is.null(ci_lower_var)) {
#         if(ci_upper_var %in% names(data) && ci_lower_var %in% names(data)) {
#           upper <- data[[ci_upper_var]]
#           lower <- data[[ci_lower_var]]
# 
#           ci_color <- adjustcolor(color_raw, alpha.f = pltParams[["alpha"]])
#           polygon(c(age, rev(age)),
#                   c(upper, rev(lower)),
#                   col = ci_color, border = NA)
#         }
#       }
# 
#       # Plot the raw line with alpha
#       raw_color_with_alpha <- adjustcolor(color_raw, alpha.f = raw_alpha)
#       lines(age, value, col = raw_color_with_alpha, lwd = pltParams[["linewidths"]][1], lty = 1)
# 
#       # Plot first smooth curve if requested
#       if(plot_smooth1 && !is.null(smooth1_var)) {
#         if(smooth1_var %in% names(data)) {
#           smooth1_value <- data[[smooth1_var]]
#           smooth1_color_with_alpha <- adjustcolor(color_smooth1, alpha.f = smooth_alpha)
#           lines(age, smooth1_value, col = smooth1_color_with_alpha, lwd = pltParams[["linewidths"]][2], lty = 1)
#         }
#       }
# 
#       # Plot second smooth curve if requested
#       if(plot_smooth2 && !is.null(smooth2_var)) {
#         if(smooth2_var %in% names(data)) {
#           smooth2_value <- data[[smooth2_var]]
#           smooth2_color_with_alpha <- adjustcolor(color_smooth2, alpha.f = smooth_alpha)
#           lines(age, smooth2_value, col = smooth2_color_with_alpha, lwd = pltParams[["linewidths"]][3], lty = 1)
#         }
#       }
#     }
# 
#     # Add legend if requested
#     if (show_legend) {
#       plot_xlim <- par("usr")[1:2]
#       plot_ylim <- par("usr")[3:4]
#       x_range <- diff(plot_xlim)
#       y_range <- diff(plot_ylim)
# 
#       # Build legend based on what's being plotted
#       legend_labels <- names(data_list)
#       legend_colors <- pltParams[["colors"]]["raw", names(data_list)]
#       legend_lwd <- rep(pltParams[["linewidths"]][1], length(data_list))
# 
#       if(plot_smooth1 && !is.null(smooth1_var)) {
#         # Add first smooth curve entries to legend
#         smooth1_labels <- paste0(names(data_list), " (", smooth1_label, ")")
#         legend_labels <- c(legend_labels, smooth1_labels)
#         legend_colors <- c(legend_colors, pltParams[["colors"]]["smooth1", names(data_list)])
#         legend_lwd <- c(legend_lwd, rep(pltParams[["linewidths"]][2], length(data_list)))
#       }
# 
#       if(plot_smooth2 && !is.null(smooth2_var)) {
#         # Add second smooth curve entries to legend
#         smooth2_labels <- paste0(names(data_list), " (", smooth2_label, ")")
#         legend_labels <- c(legend_labels, smooth2_labels)
#         legend_colors <- c(legend_colors, pltParams[["colors"]]["smooth2", names(data_list)])
#         legend_lwd <- c(legend_lwd, rep(pltParams[["linewidths"]][3], length(data_list)))
#       }
# 
#       legend(x = plot_xlim[1] + 0.02*x_range,
#              y = plot_ylim[4] - 0.02*y_range,
#              legend = legend_labels,
#              col = legend_colors,
#              lwd = legend_lwd,
#              bty = "n",
#              cex = text_size_multiplier)
#     }
#   }
# 
#   # Create plot for display
#   create_single_plot(show_legend = TRUE, show_xaxis = TRUE, yaxis_right = FALSE)
# 
#   # Save all variants if requested
#   if (printOn && save_variants) {
# 
#     variants <- list(
#       list(name = "", show_legend = TRUE, show_xaxis = TRUE, yaxis_right = FALSE),
#       list(name = "_noLeg", show_legend = FALSE, show_xaxis = TRUE, yaxis_right = FALSE),
#       list(name = "_noX", show_legend = FALSE, show_xaxis = FALSE, yaxis_right = FALSE),
#       list(name = "_noX_Yright", show_legend = FALSE, show_xaxis = FALSE, yaxis_right = TRUE),
#       list(name = "_noLeg_Yright", show_legend = FALSE, show_xaxis = TRUE, yaxis_right = TRUE)
#     )
# 
#     for (variant in variants) {
#       plot_filename <- file.path(print_dir, paste0(fname_prefix, variant$name, ".png"))
# 
#       if (save_transparent) {
#         png(plot_filename, width = fig_w, height = fig_h, units = "in", res = 600, bg = "transparent")
#       } else {
#         png(plot_filename, width = fig_w, height = fig_h, units = "in", res = 600, bg = "white")
#       }
# 
#       create_single_plot(show_legend = variant$show_legend,
#                          show_xaxis = variant$show_xaxis,
#                          yaxis_right = variant$yaxis_right)
#       dev.off()
#     }
# 
#     message("Saved ", length(variants), " variants")
#   }
# }

plot_timeseries <- function(data_list, pltParams,
                            age_var = "year",           # Name of age/time column
                            value_var = "value",        # Name of value column to plot
                            ci_upper_var = NULL,        # Name of upper CI column (optional)
                            ci_lower_var = NULL,        # Name of lower CI column (optional)
                            smooth1_var = NULL,         # Name of first smooth curve column (optional)
                            smooth2_var = NULL,         # Name of second smooth curve column (optional)
                            plot_smooth1 = FALSE,       # Whether to plot first smooth curve
                            plot_smooth2 = FALSE,       # Whether to plot second smooth curve
                            smooth1_label = "smooth",   # Label for first smooth curve in legend
                            smooth2_label = "smooth2",  # Label for second smooth curve in legend
                            ylab = "Fire emissions",    # Y-axis label
                            save_path = NULL,
                            fname_prefix = "timeseries",
                            save_transparent = TRUE,
                            save_variants = TRUE) {
  
  # Set up y-limits and x-limits
  xlims <- c(pltParams[["xmin"]], pltParams[["xmax"]])
  ylims <- c(pltParams[["ymin"]], pltParams[["ymax"]])
  
  # Create print directory
  if(is.null(save_path)) {
    print_dir <- output_dir
  } else {
    print_dir <- save_path
  }
  dir.create(print_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Function to create a single plot with options
  create_single_plot <- function(show_legend = TRUE,
                                 show_xaxis = TRUE,
                                 yaxis_right = FALSE) {
    
    # Set margins based on options
    if (yaxis_right && !show_xaxis) {
      par(mar = c(1, 1, 1, 4), cex = 0.75 * text_size_multiplier)
    } else if (!show_xaxis) {
      par(mar = c(1, 4, 1, 2), cex = 0.75 * text_size_multiplier)
    } else if (yaxis_right) {
      par(mar = c(5, 2, 4, 4), cex = 0.75 * text_size_multiplier)
    } else {
      par(mar = c(5, 4, 4, 2), cex = 0.75 * text_size_multiplier)
    }
    
    # Determine axis parameters
    yaxt_param <- if(yaxis_right) "n" else "s"
    xaxt_param <- if(!show_xaxis) "n" else "s"
    bty_param <- if(!show_xaxis) "n" else "o"
    
    plot(NA, type = "n", xlim = xlims, ylim = ylims,
         xlab = if(show_xaxis) "Year (CE)" else "",
         ylab = if(!yaxis_right) ylab else "",
         main = "",
         yaxt = yaxt_param,
         xaxt = xaxt_param,
         bty = bty_param)
    
    # Add y-axis on right if needed
    if(yaxis_right) {
      axis(4, las = 1)
      mtext(ylab, side = 4, line = 3, cex = 0.75 * text_size_multiplier)
    }
    
    # For plots without x-axis, manually draw only the y-axis line
    if(!show_xaxis) {
      usr <- par("usr")
      
      if(yaxis_right) {
        lines(c(usr[2], usr[2]), c(usr[3], usr[4]), col = "black", lwd = 1)
      } else {
        lines(c(usr[1], usr[1]), c(usr[3], usr[4]), col = "black", lwd = 1)
      }
    }
    
    # Add major and minor ticks to x-axis if shown
    if(show_xaxis) {
      major_ticks <- pretty(xlims)
      minor_ticks <- unlist(lapply(1:(length(major_ticks) - 1), function(k) {
        seq(major_ticks[k], major_ticks[k + 1], length.out = 6)[-c(1, 6)]
      }))
      
      axis(side = 1, at = major_ticks, labels = TRUE)
      axis(side = 1, at = minor_ticks, labels = FALSE, tcl = -0.3)
      axis(side = 3, at = major_ticks, labels = FALSE, tcl = 0.5)
      axis(side = 3, at = minor_ticks, labels = FALSE, tcl = 0.25)
    }
    
    # Plot data
    for(dataset_name in names(data_list)) {
      data <- data_list[[dataset_name]]
      
      # Get colors for raw and smooth curves
      color_raw <- pltParams[["colors"]]["raw", dataset_name]
      color_smooth1 <- pltParams[["colors"]]["smooth1", dataset_name]
      color_smooth2 <- pltParams[["colors"]]["smooth2", dataset_name]
      
      # Get alpha values (default to 1.0 if not set)
      raw_alpha <- if("raw_alpha" %in% names(pltParams)) pltParams[["raw_alpha"]] else 1.0
      smooth_alpha <- if("smooth_alpha" %in% names(pltParams)) pltParams[["smooth_alpha"]] else 1.0
      
      # Extract variables using flexible column names
      age <- data[[age_var]]
      value <- data[[value_var]]
      
      # Plot confidence interval if available
      if(!is.null(ci_upper_var) && !is.null(ci_lower_var)) {
        if(ci_upper_var %in% names(data) && ci_lower_var %in% names(data)) {
          upper <- data[[ci_upper_var]]
          lower <- data[[ci_lower_var]]
          
          # Filter out NAs for CI polygon
          valid_idx <- !is.na(age) & !is.na(upper) & !is.na(lower)
          if(sum(valid_idx) > 0) {
            age_valid <- age[valid_idx]
            upper_valid <- upper[valid_idx]
            lower_valid <- lower[valid_idx]
            
            ci_color <- adjustcolor(color_raw, alpha.f = pltParams[["alpha"]])
            polygon(c(age_valid, rev(age_valid)),
                    c(upper_valid, rev(lower_valid)),
                    col = ci_color, border = NA)
          }
        }
      }
      
      # Plot the raw line with alpha - filter out NAs
      valid_idx <- !is.na(age) & !is.na(value)
      if(sum(valid_idx) > 0) {
        age_valid <- age[valid_idx]
        value_valid <- value[valid_idx]
        
        raw_color_with_alpha <- adjustcolor(color_raw, alpha.f = raw_alpha)
        lines(age_valid, value_valid, col = raw_color_with_alpha, 
              lwd = pltParams[["linewidths"]][1], lty = 1)
      }
      
      # Plot first smooth curve if requested - filter out NAs
      if(plot_smooth1 && !is.null(smooth1_var)) {
        if(smooth1_var %in% names(data)) {
          smooth1_value <- data[[smooth1_var]]
          
          valid_idx <- !is.na(age) & !is.na(smooth1_value)
          if(sum(valid_idx) > 0) {
            age_valid <- age[valid_idx]
            smooth1_valid <- smooth1_value[valid_idx]
            
            smooth1_color_with_alpha <- adjustcolor(color_smooth1, alpha.f = smooth_alpha)
            lines(age_valid, smooth1_valid, col = smooth1_color_with_alpha, 
                  lwd = pltParams[["linewidths"]][2], lty = 1)
          }
        }
      }
      
      # Plot second smooth curve if requested - filter out NAs
      if(plot_smooth2 && !is.null(smooth2_var)) {
        if(smooth2_var %in% names(data)) {
          smooth2_value <- data[[smooth2_var]]
          
          valid_idx <- !is.na(age) & !is.na(smooth2_value)
          if(sum(valid_idx) > 0) {
            age_valid <- age[valid_idx]
            smooth2_valid <- smooth2_value[valid_idx]
            
            smooth2_color_with_alpha <- adjustcolor(color_smooth2, alpha.f = smooth_alpha)
            lines(age_valid, smooth2_valid, col = smooth2_color_with_alpha, 
                  lwd = pltParams[["linewidths"]][3], lty = 1)
          }
        }
      }
    }
    
    # Add legend if requested
    if (show_legend) {
      plot_xlim <- par("usr")[1:2]
      plot_ylim <- par("usr")[3:4]
      x_range <- diff(plot_xlim)
      y_range <- diff(plot_ylim)
      
      # Build legend based on what's being plotted
      legend_labels <- names(data_list)
      legend_colors <- pltParams[["colors"]]["raw", names(data_list)]
      legend_lwd <- rep(pltParams[["linewidths"]][1], length(data_list))
      
      if(plot_smooth1 && !is.null(smooth1_var)) {
        # Add first smooth curve entries to legend
        smooth1_labels <- paste0(names(data_list), " (", smooth1_label, ")")
        legend_labels <- c(legend_labels, smooth1_labels)
        legend_colors <- c(legend_colors, pltParams[["colors"]]["smooth1", names(data_list)])
        legend_lwd <- c(legend_lwd, rep(pltParams[["linewidths"]][2], length(data_list)))
      }
      
      if(plot_smooth2 && !is.null(smooth2_var)) {
        # Add second smooth curve entries to legend
        smooth2_labels <- paste0(names(data_list), " (", smooth2_label, ")")
        legend_labels <- c(legend_labels, smooth2_labels)
        legend_colors <- c(legend_colors, pltParams[["colors"]]["smooth2", names(data_list)])
        legend_lwd <- c(legend_lwd, rep(pltParams[["linewidths"]][3], length(data_list)))
      }
      
      legend(x = plot_xlim[1] + 0.02*x_range,
             y = plot_ylim[4] - 0.02*y_range,
             legend = legend_labels,
             col = legend_colors,
             lwd = legend_lwd,
             bty = "n",
             cex = text_size_multiplier)
    }
  }
  
  # Create plot for display
  create_single_plot(show_legend = TRUE, show_xaxis = TRUE, yaxis_right = FALSE)
  
  # Save all variants if requested
  if (printOn && save_variants) {
    
    variants <- list(
      list(name = "", show_legend = TRUE, show_xaxis = TRUE, yaxis_right = FALSE),
      list(name = "_noLeg", show_legend = FALSE, show_xaxis = TRUE, yaxis_right = FALSE),
      list(name = "_noX", show_legend = FALSE, show_xaxis = FALSE, yaxis_right = FALSE),
      list(name = "_noX_Yright", show_legend = FALSE, show_xaxis = FALSE, yaxis_right = TRUE),
      list(name = "_noLeg_Yright", show_legend = FALSE, show_xaxis = TRUE, yaxis_right = TRUE)
    )
    
    for (variant in variants) {
      plot_filename <- file.path(print_dir, paste0(fname_prefix, variant$name, ".png"))
      
      if (save_transparent) {
        png(plot_filename, width = fig_w, height = fig_h, units = "in", res = 600, bg = "transparent")
      } else {
        png(plot_filename, width = fig_w, height = fig_h, units = "in", res = 600, bg = "white")
      }
      
      create_single_plot(show_legend = variant$show_legend,
                         show_xaxis = variant$show_xaxis,
                         yaxis_right = variant$yaxis_right)
      dev.off()
    }
    
    message("Saved ", length(variants), " variants")
  }
}



# ------------------------- plot_age_range_boxplots ----------------------------- ####

plot_age_range_violins <- function(data_list,
                                   age_ranges,
                                   age_var = "year_ce",
                                   value_vars,
                                   normalize = TRUE,
                                   ylab = "Normalized value",
                                   main_title = NULL,
                                   save_path = NULL,
                                   fname_prefix = "age_violin",
                                   fig_width = 3,
                                   fig_height = 5,
                                   violin_color = "#1473CC",
                                   violin_alpha = 0.6,
                                   add_jitter = TRUE,
                                   jitter_color = "black",
                                   jitter_alpha = 0.3,
                                   jitter_size = 0.5,
                                   jitter_amount = 0.1,
                                   ylims = NULL) {
  
  # Load required library
  if(!require(vioplot)) {
    install.packages("vioplot")
    library(vioplot)
  }
  
  # Helper function for min-max normalization (0 to 1)
  min_max_normalize <- function(x) {
    x_clean <- x[!is.na(x)]
    if(length(x_clean) == 0) return(rep(NA, length(x)))
    min_val <- min(x_clean)
    max_val <- max(x_clean)
    if(min_val == max_val) return(rep(0.5, length(x)))
    return((x - min_val) / (max_val - min_val))
  }
  
  # Create save directory if needed
  if(!is.null(save_path)) {
    dir.create(save_path, recursive = TRUE, showWarnings = FALSE)
  }
  
  # Improved check: determine if this is a list of lists or list of dataframes
  # Check if ALL elements are dataframes (not just the first one)
  all_are_dataframes <- all(sapply(data_list, is.data.frame))
  
  if(all_are_dataframes) {
    # This is a list of dataframes - each should be a separate group
    # DON'T wrap them - they're already structured correctly
    is_single_group <- FALSE
    cat("Detected list of dataframes - treating each as a separate group\n")
  } else {
    # This is a list of lists - each sub-list is a group
    is_single_group <- FALSE
    cat("Detected list of lists - treating each sub-list as a group\n")
  }
  
  # Store all summary statistics
  all_summary_stats <- list()
  
  # Process each top-level group separately
  for(group_name in names(data_list)) {
    
    group_data <- data_list[[group_name]]
    
    cat("\n=== Processing group:", group_name, "===\n")
    
    # For this group, process each age range
    pooled_data <- list()
    
    for(i in 1:nrow(age_ranges)) {
      range_name <- age_ranges$name[i]
      range_min <- age_ranges$min[i]
      range_max <- age_ranges$max[i]
      range_label <- paste0(range_min, "-", range_max)
      
      # Collect all values for this age range from all datasets in this group
      all_values_in_range <- c()
      
      # Handle both cases: group_data is a list of dataframes OR a single dataframe
      if(is.data.frame(group_data)) {
        # Single dataframe case (e.g., fire scar or charcoal for one region)
        datasets_to_process <- list(single = group_data)
      } else {
        # List of dataframes case (e.g., multiple ice cores in one region)
        datasets_to_process <- group_data
      }
      
      for(dataset_name in names(datasets_to_process)) {
        df <- datasets_to_process[[dataset_name]]
        
        # Check if age variable exists
        if(!age_var %in% names(df)) {
          warning(paste("Age variable", age_var, "not found in", dataset_name, "of group", group_name))
          next
        }
        
        # Process each value variable
        for(value_var in value_vars) {
          if(!value_var %in% names(df)) {
            warning(paste("Value variable", value_var, "not found in", dataset_name, "of group", group_name))
            next
          }
          
          # Extract age and value
          age <- df[[age_var]]
          value <- df[[value_var]]
          
          # Remove NAs
          valid_idx <- !is.na(age) & !is.na(value)
          age_valid <- age[valid_idx]
          value_valid <- value[valid_idx]
          
          if(length(age_valid) == 0) next
          
          # Normalize THIS time series separately
          if(normalize) {
            value_normalized <- min_max_normalize(value_valid)
          } else {
            value_normalized <- value_valid
          }
          
          # Find values in this age range from this normalized time series
          in_range <- age_valid >= range_min & age_valid <= range_max
          values_in_range <- value_normalized[in_range]
          
          # Add to the pool for this age range
          if(length(values_in_range) > 0) {
            all_values_in_range <- c(all_values_in_range, values_in_range)
          }
        }
      }
      
      # Store pooled data for this age range
      if(length(all_values_in_range) > 0) {
        pooled_data[[range_name]] <- list(
          values = all_values_in_range,
          age_range = range_name,
          range_label = range_label,
          group_name = group_name
        )
      }
    }
    
    if(length(pooled_data) == 0) {
      warning(paste("No valid data found for group:", group_name))
      next
    }
    
    # Determine y-limits if not provided (use group-specific or global)
    if(is.null(ylims)) {
      all_values <- unlist(lapply(pooled_data, function(x) x$values))
      current_ylims <- c(min(all_values, na.rm = TRUE) - 0.05, 
                         max(all_values, na.rm = TRUE) + 0.05)
    } else {
      current_ylims <- ylims
    }
    
    # Create a separate figure for each age range in this group
    for(range_name in names(pooled_data)) {
      
      item <- pooled_data[[range_name]]
      
      # Get components
      current_age_range <- item$age_range
      current_range_label <- item$range_label
      
      # Prepare data for single violin
      plot_data <- list(item$values)
      plot_label <- current_range_label
      
      # Create plot title
      plot_title <- if(!is.null(main_title)) {
        paste0(main_title, "\n", group_name, "\n", current_range_label)
      } else {
        paste0(group_name, "\n", current_range_label)
      }
      
      # Create filename - always use group_name since we never have single_group anymore
      plot_filename <- paste0(fname_prefix, "_", group_name, "_", current_age_range, ".png")
      
      # Open graphics device if saving
      if(!is.null(save_path)) {
        png(file.path(save_path, plot_filename),
            width = fig_width, height = fig_height, units = "in", res = 300)
      }
      
      # Set up plot margins - increase left margin for y-axis label
      par(mar = c(4, 5, 4, 2), cex = 0.8)
      
      # Create the violin plot without y-axis label
      # Use adjustcolor to apply alpha to the color
      fill_color <- adjustcolor(violin_color, alpha.f = violin_alpha)
      
      vioplot(plot_data,
              names = plot_label,
              col = fill_color,
              border = violin_color,
              main = plot_title,
              ylab = "",  # Remove default y-axis label
              las = 1,
              ylim = current_ylims,
              cex.axis = 0.8,
              cex.names = 0.8,
              cex.main = 0.9)
      
      # Add y-axis label with more spacing (line = 3.5 instead of default ~2)
      mtext(ylab, side = 2, line = 3.5, cex = 0.8)
      
      # Add jittered points if requested
      if(add_jitter) {
        jitter_col <- adjustcolor(jitter_color, alpha.f = jitter_alpha)
        
        # Add jitter to x position (centered at 1 since it's a single violin)
        x_pos <- rep(1, length(plot_data[[1]]))
        x_jittered <- x_pos + runif(length(x_pos), -jitter_amount, jitter_amount)
        
        # Plot points
        points(x_jittered, plot_data[[1]], 
               col = jitter_col, 
               pch = 16, 
               cex = jitter_size)
      }
      
      # Add grid for readability
      grid(nx = NA, ny = NULL, col = "gray90", lty = "dotted")
      
      # Redraw violin plot on top of grid
      vioplot(plot_data,
              names = plot_label,
              col = fill_color,
              border = violin_color,
              las = 1,
              ylim = current_ylims,
              cex.axis = 0.8,
              cex.names = 0.8,
              add = TRUE)
      
      # Re-add jittered points on top
      if(add_jitter) {
        jitter_col <- adjustcolor(jitter_color, alpha.f = jitter_alpha)
        
        x_pos <- rep(1, length(plot_data[[1]]))
        x_jittered <- x_pos + runif(length(x_pos), -jitter_amount, jitter_amount)
        
        points(x_jittered, plot_data[[1]], 
               col = jitter_col, 
               pch = 16, 
               cex = jitter_size)
      }
      
      # Close graphics device if saving
      if(!is.null(save_path)) {
        dev.off()
        message("Saved violin plot to: ", file.path(save_path, plot_filename))
      }
    }
    
    # Store summary statistics for this group
    group_summary <- lapply(names(pooled_data), function(range_name) {
      item <- pooled_data[[range_name]]
      data.frame(
        group = group_name,
        age_range = item$age_range,
        range_label = item$range_label,
        n = length(item$values),
        mean = mean(item$values, na.rm = TRUE),
        median = median(item$values, na.rm = TRUE),
        sd = sd(item$values, na.rm = TRUE),
        min = min(item$values, na.rm = TRUE),
        max = max(item$values, na.rm = TRUE)
      )
    })
    
    all_summary_stats[[group_name]] <- do.call(rbind, group_summary)
  }
  
  # Combine all summary statistics
  summary_df <- do.call(rbind, all_summary_stats)
  rownames(summary_df) <- NULL
  
  return(summary_df)
}

# =============================== Load in data ================================ ####


# ----------------------- Load in the inverted emissions ------------------------- #
ice.model <- read.csv(file.path(data_dir, "Zhang_BB_emission_BB4CMIP_modeled.csv"), header = TRUE)

# Get all column names except year_ce
emission_cols <- setdiff(names(ice.model), "year_ce")

# Create a list of dataframes, one for each emission column
ice.model.list <- lapply(emission_cols, function(col) {
  df <- data.frame(
    year_ce = ice.model$year_ce,
    bb_emissions = ice.model[[col]]
  )
  return(df)
})

# Name the list elements with the original column names
names(ice.model.list) <- emission_cols

ice.data <- read.csv(file.path(data_dir, "Zhang_BB_emission_recon.csv"), header = TRUE)

# Get all column names except year_ce
emission_cols <- setdiff(names(ice.data), "year_ce")

# Create a list of dataframes, one for each emission column
ice.data.list <- lapply(emission_cols, function(col) {
  df <- data.frame(
    year_ce = ice.data$year_ce,
    bb_emissions = ice.data[[col]]
  )
  return(df)
})

# Name the list elements with the original column names
names(ice.data.list) <- emission_cols

# ------------------------ Load in the ice core records -------------------------- #

ice_core.data <- read.csv(file.path(data_dir, '1750_to_2010_global_rBC_deposition.csv'), header = TRUE)
ice_core.meta <- read.csv(file.path(data_dir, '1750_to_2010_global_ice_core_array_information.csv'), header = TRUE)

# Replace -999 with NA
ice_core.data <- ice_core.data %>% 
  mutate(across(everything(), ~na_if(., -999)))

# Separate the ice core records into regions 
# ice_core.boreal_west_sites <- c("McCall.Glacier", "Eclipse") # local only
ice_core.boreal_west_sites <- c("McCall.Glacier", "Eclipse","Humboldt", "NGT.B19", "Tunu2013","Hans.Tausen","NEEM.2011.S1","NasaU", "Summit2010", "D4","ACT2", "ACT11d") #include all of greenland
ice_core.boreal_east_sites <- c("Humboldt", "NGT.B19", "Tunu2013","Hans.Tausen","NEEM.2011.S1","NasaU", "Summit2010", "D4","ACT2", "ACT11d")
ice_core.temperate_west_sites <- c("Upper.Fremont.Glacier")

# Create lists of ice core data by region
ice_core.boreal_west_list <- lapply(ice_core.boreal_west_sites, function(site) {
  if(site %in% names(ice_core.data)) {
    data.frame(
      year_ce = ice_core.data$Year,
      rBC_deposition = ice_core.data[[site]]
    )
  }
})
names(ice_core.boreal_west_list) <- ice_core.boreal_west_sites

ice_core.boreal_east_list <- lapply(ice_core.boreal_east_sites, function(site) {
  if(site %in% names(ice_core.data)) {
    data.frame(
      year_ce = ice_core.data$Year,
      rBC_deposition = ice_core.data[[site]]
    )
  }
})
names(ice_core.boreal_east_list) <- ice_core.boreal_east_sites

ice_core.temperate_west_list <- lapply(ice_core.temperate_west_sites, function(site) {
  if(site %in% names(ice_core.data)) {
    data.frame(
      year_ce = ice_core.data$Year,
      rBC_deposition = ice_core.data[[site]]
    )
  }
})
names(ice_core.temperate_west_list) <- ice_core.temperate_west_sites

# # Apply z-scoring to each region's ice core data
# ice_core.boreal_west_list <- zscore_normalize(ice_core.boreal_west_list)
# ice_core.boreal_east_list <- zscore_normalize(ice_core.boreal_east_list)
# ice_core.temperate_west_list <- zscore_normalize(ice_core.temperate_west_list)

# Apply z-scoring to each region's ice core data
ice_core.boreal_west_list <- zscore_normalize(ice_core.boreal_west_list)
ice_core.boreal_east_list <- zscore_normalize(ice_core.boreal_east_list)
ice_core.temperate_west_list <- zscore_normalize(ice_core.temperate_west_list)

if(iceSmoothType == 'indiv'){
  # Apply smoothing to each region's ice core data
  ice_core.boreal_west_list <- smooth_ice_core_series(ice_core.boreal_west_list,
                                                      hw1 = smooth_params$hw1,
                                                      hw2 = smooth_params$hw2,
                                                      degree = smooth_params$degree,
                                                      family = smooth_params$family)
  
  ice_core.boreal_east_list <- smooth_ice_core_series(ice_core.boreal_east_list,
                                                      hw1 = smooth_params$hw1,
                                                      hw2 = smooth_params$hw2,
                                                      degree = smooth_params$degree,
                                                      family = smooth_params$family)
  
  ice_core.temperate_west_list <- smooth_ice_core_series(ice_core.temperate_west_list,
                                                         hw1 = smooth_params$hw1,
                                                         hw2 = smooth_params$hw2,
                                                         degree = smooth_params$degree,
                                                         family = smooth_params$family)
} else if(iceSmoothType == 'all'){
  # Apply smoothing to POOLED data from each region (fit one curve per region)
  ice_core.boreal_west_list <- smooth_ice_core_collection(ice_core.boreal_west_list,
                                                          hw1 = smooth_params$hw1,
                                                          hw2 = smooth_params$hw2,
                                                          degree = smooth_params$degree,
                                                          family = smooth_params$family)
  
  cat("\n=== Smoothing Boreal East ===\n")
  ice_core.boreal_east_list <- smooth_ice_core_collection(ice_core.boreal_east_list,
                                                          hw1 = smooth_params$hw1,
                                                          hw2 = smooth_params$hw2,
                                                          degree = smooth_params$degree,
                                                          family = smooth_params$family)
  
  cat("\n=== Smoothing Temperate West ===\n")
  ice_core.temperate_west_list <- smooth_ice_core_collection(ice_core.temperate_west_list,
                                                             hw1 = smooth_params$hw1,
                                                             hw2 = smooth_params$hw2,
                                                             degree = smooth_params$degree,
                                                             family = smooth_params$family)
}


# ------------------------- Load in the burn scar data --------------------------- #

scar.boreal_west <- read.csv(file.path(data_dir, "fire_scars_Boreal_NA_west.csv"), header = TRUE)
scar.boreal_east <- read.csv(file.path(data_dir, "fire_scars_Boreal_NA_east.csv"), header = TRUE)
scar.temperate_west <- read.csv(file.path(data_dir, "fire_scars_Temperate_NA_west.csv"), header = TRUE)
scar.temperate_east <- read.csv(file.path(data_dir, "fire_scars_Temperate_NA_east.csv"), header = TRUE)

scar.list <- list(boreal_west = scar.boreal_west,
                  boreal_east = scar.boreal_east,
                  temperate_west = scar.temperate_west,
                  temperate_east = scar.temperate_east)

# ------------------- Load in the charcoal synthesis curves ---------------------- #

if(charSmoothType == 'wide'){
  char.boreal_west = read.csv(file.path(data_dir, 'AERO_NAgfed_Boreal_west_PI_particle_wide_zscore.csv'), header = TRUE)
  char.boreal_east = read.csv(file.path(data_dir, 'AERO_NAgfed_Boreal_east_PI_particle_wide_zscore.csv'), header = TRUE)
  char.temp_west = read.csv(file.path(data_dir, 'AERO_NAgfed_Temperate_west_PI_particle_wide_zscore.csv'), header = TRUE)
  char.temp_east = read.csv(file.path(data_dir, 'AERO_NAgfed_Temperate_east_PI_particle_wide_zscore.csv'), header = TRUE)
} else if(charSmoothType == 'narrow'){
  char.boreal_west = read.csv(file.path(data_dir, 'AERO_NAgfed_Boreal_west_PI_particle_narrow_zscore.csv'), header = TRUE)
  char.boreal_east = read.csv(file.path(data_dir, 'AERO_NAgfed_Boreal_east_PI_particle_narrow_zscore.csv'), header = TRUE)
  char.temp_west = read.csv(file.path(data_dir, 'AERO_NAgfed_Temperate_west_PI_particle_narrow_zscore.csv'), header = TRUE)
  char.temp_east = read.csv(file.path(data_dir, 'AERO_NAgfed_Temperate_east_PI_particle_narrow_zscore.csv'), header = TRUE)
}


# ======================= Create master list of lists ========================= ####

# data.full.list <- list(
#   ice.model = ice.model.list,
#   ice.data = ice.data.list,
#   scar = scar.list
# )

if(iceMode == 'inversion') {
  data.full.list <- list(
    ice.model = ice.model.list,
    ice.data = ice.data.list,
    scar = scar.list
  )
} else if(iceMode == 'raw') {
  data.full.list <- list(
    ice_core.boreal_west = ice_core.boreal_west_list,
    ice_core.boreal_east = ice_core.boreal_east_list,
    ice_core.temperate_west = ice_core.temperate_west_list,
    scar = scar.list
  )
}

# ============================ Set up plot params ============================= ####

# Base plotting parameters (shared across all plots)
pltParams.base <- list()

pltParams.base[["linewidths"]] <- c(1.5, 1.5, 1.5)  # Line width for raw, smooth1, smooth2
pltParams.base[["alpha"]] <- c(0.2)  # Alpha for confidence intervals

# Define x-axis limits (in years CE)
pltParams.base[["xmin"]] <- 1750
pltParams.base[["xmax"]] <- 2050

# Define y-axis limits (adjust based on your data)
pltParams.base[["ymin"]] <- 0
pltParams.base[["ymax"]] <- 100  # Adjust based on your emission data range

# Dataset-specific parameters (colors, labels, etc.)
pltParams.specific <- list(
  
  # Ice core model parameters
  ice.model = list(
    colors = list(
      Boreal.NA.East = data.frame(
        raw = "#1473CC", smooth1 = "#6EB1F1", smooth2 = "#A8D5F7",
        row.names = c("Boreal.NA.East")
      ),
      Boreal.NA.West = data.frame(
        raw = "#0D5A9E", smooth1 = "#5B9DC9", smooth2 = "#8FC1E3",
        row.names = c("Boreal.NA.West")
      ),
      Temperate.NA.East = data.frame(
        raw = "#2E8B57", smooth1 = "#66CDAA", smooth2 = "#98FB98",
        row.names = c("Temperate.NA.East")
      ),
      Temperate.NA.West = data.frame(
        raw = "#1F6B42", smooth1 = "#4FA375", smooth2 = "#7AC69A",
        row.names = c("Temperate.NA.West")
      ),
      Total.NA = data.frame(
        raw = "#4B0082", smooth1 = "#8A2BE2", smooth2 = "#BA55D3",
        row.names = c("Total.NA")
      )
    ),
    single_color = "#F9B5AC",
    ylims = list(
      Boreal.NA.East = c(0, 0.07),
      Boreal.NA.West = c(0, 0.1),
      Temperate.NA.East = c(0, 0.03),
      Temperate.NA.West = c(0, 0.03),
      Total.NA = c(0, 0.15)
    ),
    ylab = "Modeled BC Emissions (Tg C/yr)",
    fname_suffix = "model",
    age_var = "year_ce",
    value_var = "bb_emissions",
    ci_upper_var = NULL,
    ci_lower_var = NULL,
    smooth1_var = NULL,
    smooth2_var = NULL,
    plot_smooth1 = FALSE,
    plot_smooth2 = FALSE
  ),
  
  # Ice core data parameters
  ice.data = list(
    colors = list(
      Boreal.NA.East = data.frame(
        raw = "#A75E39", smooth1 = "#C67C58", smooth2 = "#E5A584",
        row.names = c("Boreal.NA.East")
      ),
      Boreal.NA.West = data.frame(
        raw = "#8B4513", smooth1 = "#A0623C", smooth2 = "#C08866",
        row.names = c("Boreal.NA.West")
      ),
      Temperate.NA.East = data.frame(
        raw = "#D2691E", smooth1 = "#E89A5C", smooth2 = "#F5C18E",
        row.names = c("Temperate.NA.East")
      ),
      Temperate.NA.West = data.frame(
        raw = "#A0522D", smooth1 = "#B87B5A", smooth2 = "#D4A78A",
        row.names = c("Temperate.NA.West")
      ),
      Total.NA = data.frame(
        raw = "#8B0000", smooth1 = "#B22222", smooth2 = "#DC143C",
        row.names = c("Total.NA")
      )
    ),
    # single_color = "#479EEB",
    single_color = "#5D2E8C",
    # single_color = "#D4DFC7",
    ylims = list(
      Boreal.NA.East = c(0, 0.07),
      Boreal.NA.West = c(0, 0.1),
      Temperate.NA.East = c(0, 0.03),
      Temperate.NA.West = c(0, 0.03),
      Total.NA = c(0, 0.15)
    ),
    ylab = "BC Emissions (Tg C/yr)",
    fname_suffix = "data",
    age_var = "year_ce",
    value_var = "bb_emissions",
    ci_upper_var = NULL,
    ci_lower_var = NULL,
    smooth1_var = NULL,
    smooth2_var = NULL,
    plot_smooth1 = FALSE,
    plot_smooth2 = FALSE
  ),
  
  # Raw ice core data parameters
  ice_core.boreal_west = list(
    single_color = list(raw = "#A8D5F7", smooth1 = "#6EB1F1", smooth2 = "#1473CC"),  # Light to dark blue
    ylims = list(default = c(-3, 3)),
    ylab = "rBC deposition (z-score)",
    fname_suffix = "ice_core_raw_zscore_smoothed",
    age_var = "year_ce",
    value_var = "rBC_deposition",
    ci_upper_var = NULL,
    ci_lower_var = NULL,
    smooth1_var = "rBC_deposition_smooth1",
    smooth2_var = "rBC_deposition_smooth2",
    plot_smooth1 = TRUE,
    plot_smooth2 = TRUE
  ),
  
  ice_core.boreal_east = list(
    single_color = list(raw = "#A8D5F7", smooth1 = "#6EB1F1", smooth2 = "#1473CC"),
    ylims = list(default = c(-3, 3)),
    ylab = "rBC deposition (z-score)",
    fname_suffix = "ice_core_raw_zscore_smoothed",
    age_var = "year_ce",
    value_var = "rBC_deposition",
    ci_upper_var = NULL,
    ci_lower_var = NULL,
    smooth1_var = "rBC_deposition_smooth1",
    smooth2_var = "rBC_deposition_smooth2",
    plot_smooth1 = TRUE,
    plot_smooth2 = TRUE
  ),
  
  ice_core.temperate_west = list(
    single_color = list(raw = "#A8D5F7", smooth1 = "#6EB1F1", smooth2 = "#1473CC"),
    ylims = list(default = c(-3, 3)),
    ylab = "rBC deposition (z-score)",
    fname_suffix = "ice_core_raw_zscore_smoothed",
    age_var = "year_ce",
    value_var = "rBC_deposition",
    ci_upper_var = NULL,
    ci_lower_var = NULL,
    smooth1_var = "rBC_deposition_smooth1",
    smooth2_var = "rBC_deposition_smooth2",
    plot_smooth1 = TRUE,
    plot_smooth2 = TRUE
  ),
  
  
  # Fire scar parameters
  scar = list(
    colors = list(
      boreal_west = data.frame(
        raw = "#FF6347", smooth1 = "#FF7F66", smooth2 = "#FFA590",
        row.names = c("boreal_west")
      ),
      boreal_east = data.frame(
        raw = "#FF4500", smooth1 = "#FF6B3D", smooth2 = "#FF9473",
        row.names = c("boreal_east")
      ),
      temperate_west = data.frame(
        raw = "#FFD700", smooth1 = "#FFE34D", smooth2 = "#FFED85",
        row.names = c("temperate_west")
      ),
      temperate_east = data.frame(
        raw = "#FFA500", smooth1 = "#FFBB3D", smooth2 = "#FFD273",
        row.names = c("temperate_east")
      )
    ),
    single_color = if(smooth_scar_data) {
      list(raw = "#91CAA2", smooth1 = "#75BD8B", smooth2 = "#3C7C4F")
    } else {
      "#3C7C4F"
    },
    ylims = list(
      boreal_west = c(0, 20),
      boreal_east = c(0, 20),
      temperate_west = c(0, 20),
      temperate_east = c(0, 20)
    ),
    ylab = "% sites with fire scar",
    fname_suffix = "scar",
    age_var = "year_ce",
    value_var = "percent_scarred",
    ci_upper_var = NULL,
    ci_lower_var = NULL,
    smooth1_var = if(smooth_scar_data) "percent_scarred_smooth1" else NULL,
    smooth2_var = if(smooth_scar_data) "percent_scarred_smooth2" else NULL,
    plot_smooth1 = smooth_scar_data,
    plot_smooth2 = smooth_scar_data
  )
)

if(smooth_scar_data){
  # Apply smoothing to each region
  scar.boreal_west <- smooth_scar_series(scar.boreal_west, 
                                         hw1 = smooth_params$hw1,
                                         hw2 = smooth_params$hw2,
                                         degree = smooth_params$degree,
                                         family = smooth_params$family)
  
  scar.boreal_east <- smooth_scar_series(scar.boreal_east, 
                                         hw1 = smooth_params$hw1,
                                         hw2 = smooth_params$hw2,
                                         degree = smooth_params$degree,
                                         family = smooth_params$family)
  
  scar.temperate_west <- smooth_scar_series(scar.temperate_west, 
                                            hw1 = smooth_params$hw1,
                                            hw2 = smooth_params$hw2,
                                            degree = smooth_params$degree,
                                            family = smooth_params$family)
  
  scar.temperate_east <- smooth_scar_series(scar.temperate_east, 
                                            hw1 = smooth_params$hw1,
                                            hw2 = smooth_params$hw2,
                                            degree = smooth_params$degree,
                                            family = smooth_params$family)
  
  # Update the scar list with smoothed data
  scar.list <- list(boreal_west = scar.boreal_west,
                    boreal_east = scar.boreal_east,
                    temperate_west = scar.temperate_west,
                    temperate_east = scar.temperate_east)
  
  # Update data.full.list
  data.full.list$scar <- scar.list
}

# ============================== Make the plot ================================ ####

# Add a toggle for using single color
use_single_color <- TRUE  # Set to TRUE to use single color per dataset

# Loop through each dataset type (ice.model, ice.data, scar, OR ice_core regions)
for(dataset_type in names(data.full.list)) {
  
  cat("\n=== Processing", dataset_type, "===\n")
  
  # Get the list of regions/datasets for this type
  region_list <- data.full.list[[dataset_type]]
  
  # Get specific parameters for this dataset type
  specific_params <- pltParams.specific[[dataset_type]]
  
  # Check if this is raw ice core data (multiple sites per region)
  if(grepl("ice_core", dataset_type)) {
    
    # For raw ice core mode: plot all sites in the region together
    cat("Creating plot for region:", dataset_type, "with", length(region_list), "sites\n")
    
    # Use the entire list of sites for this region
    data_list <- region_list
    
    # Remove any NULL entries (sites not in data)
    data_list <- data_list[!sapply(data_list, is.null)]
    
    if(length(data_list) == 0) {
      cat("No data available for", dataset_type, ", skipping...\n")
      next
    }
    
    # Combine base parameters with specific parameters
    pltParams <- pltParams.base
    
    # Use default ylims for ice core data
    pltParams[["ymin"]] <- specific_params$ylims$default[1]
    pltParams[["ymax"]] <- specific_params$ylims$default[2]
    
    # Create colors for each site
    # Create colors for each site
    num_sites <- length(data_list)
    if(use_single_color) {
      # Use single color scheme (light for raw, darker for smoothed)
      colors_raw <- rep(specific_params$single_color$raw, num_sites)
      colors_smooth1 <- rep(specific_params$single_color$smooth1, num_sites)
      colors_smooth2 <- rep(specific_params$single_color$smooth2, num_sites)
    } else {
      # Use color gradient for different sites
      colors_raw <- colorRampPalette(c("#A8D5F7", "#87CEEB"))(num_sites)
      colors_smooth1 <- colorRampPalette(c("#6EB1F1", "#4A9FD8"))(num_sites)
      colors_smooth2 <- colorRampPalette(c("#1473CC", "#0D5A9E"))(num_sites)
    }
    
    # Create color dataframe (rows = line types, cols = sites)
    single_color_df <- data.frame(
      raw = colors_raw,
      smooth1 = colors_smooth1,
      smooth2 = colors_smooth2
    )
    single_color_df <- as.data.frame(t(single_color_df))
    colnames(single_color_df) <- names(data_list)
    rownames(single_color_df) <- c("raw", "smooth1", "smooth2")
    pltParams[["colors"]] <- single_color_df
    
    # Set alpha values
    pltParams[["raw_alpha"]] <- 0.25
    pltParams[["smooth_alpha"]] <- 1.0
    
    # Create filename prefix
    fname_prefix <- paste0(dataset_type, "_", specific_params$fname_suffix)
    
    # Create the plots
    plot_timeseries(data_list, 
                    pltParams, 
                    age_var = specific_params$age_var,
                    value_var = specific_params$value_var,
                    ci_upper_var = specific_params$ci_upper_var,
                    ci_lower_var = specific_params$ci_lower_var,
                    smooth1_var = specific_params$smooth1_var,
                    smooth2_var = specific_params$smooth2_var,
                    plot_smooth1 = specific_params$plot_smooth1,
                    plot_smooth2 = specific_params$plot_smooth2,
                    ylab = specific_params$ylab,
                    save_path = output_dir,
                    fname_prefix = fname_prefix,
                    save_transparent = TRUE,
                    save_variants = TRUE)
    
  } else {
    
    # Original behavior for inversion mode and scar data
    # Loop through each region in this dataset type
    for(region_name in names(region_list)) {
      
      cat("Creating plot for:", region_name, "\n")
      
      # Create the data list for plot_timeseries (needs to be a list with named elements)
      data_list <- list(region_list[[region_name]])
      names(data_list) <- region_name
      
      # Combine base parameters with specific parameters
      pltParams <- pltParams.base
      
      # Choose between single color or region-specific colors
      if(use_single_color) {
        # Special handling for scar data with smoothing
        if(dataset_type == "scar" && smooth_scar_data) {
          
          # Create a dataframe with proper structure (rows = line types, cols = region)
          single_color_df <- data.frame(
            region = c(specific_params$single_color$raw,
                       specific_params$single_color$smooth1,
                       specific_params$single_color$smooth2)
          )
          colnames(single_color_df) <- region_name
          rownames(single_color_df) <- c("raw", "smooth1", "smooth2")
          pltParams[["colors"]] <- single_color_df
          
          # Set alpha for raw data
          pltParams[["raw_alpha"]] <- 0.25
          pltParams[["smooth_alpha"]] <- 1.0
          
        } else {
          # Standard single color dataframe
          single_color_value <- if(is.list(specific_params$single_color)) {
            specific_params$single_color$raw
          } else {
            specific_params$single_color
          }
          
          single_color_df <- data.frame(
            col = rep(single_color_value, 3)
          )
          colnames(single_color_df) <- region_name
          rownames(single_color_df) <- c("raw", "smooth1", "smooth2")
          pltParams[["colors"]] <- single_color_df
          
          # Default alpha values
          pltParams[["raw_alpha"]] <- 1.0
          pltParams[["smooth_alpha"]] <- 1.0
        }
      } else {
        # Use region-specific colors
        pltParams[["colors"]] <- specific_params$colors[[region_name]]
        pltParams[["raw_alpha"]] <- 1.0
        pltParams[["smooth_alpha"]] <- 1.0
      }
      
      # Set region-specific y-limits
      pltParams[["ymin"]] <- specific_params$ylims[[region_name]][1]
      pltParams[["ymax"]] <- specific_params$ylims[[region_name]][2]
      
      # Create filename prefix
      if(use_single_color) {
        fname_prefix <- paste0("BB_", specific_params$fname_suffix, "_", region_name, "_singlecolor")
      } else {
        fname_prefix <- paste0("BB_", specific_params$fname_suffix, "_", region_name)
      }
      
      # Add smooth suffix to filename if smoothing is enabled for scar data
      if(dataset_type == "scar" && smooth_scar_data) {
        fname_prefix <- paste0(fname_prefix, "_smoothed")
      }
      
      # Create the plots with dataset-specific variable names
      plot_timeseries(data_list, 
                      pltParams, 
                      age_var = specific_params$age_var,
                      value_var = specific_params$value_var,
                      ci_upper_var = specific_params$ci_upper_var,
                      ci_lower_var = specific_params$ci_lower_var,
                      smooth1_var = specific_params$smooth1_var,
                      smooth2_var = specific_params$smooth2_var,
                      plot_smooth1 = specific_params$plot_smooth1,
                      plot_smooth2 = specific_params$plot_smooth2,
                      ylab = specific_params$ylab,
                      save_path = output_dir,
                      fname_prefix = fname_prefix,
                      save_transparent = TRUE,
                      save_variants = TRUE)
    }
  }
}


#### Create comparison box plots 

# Example usage:
# Define age ranges
age_ranges <- data.frame(
  name = c("Pre-industrial (1750-1780)", "Present day (1997-2010)"),
  min = c(1750, 1997),
  max = c(1780, 2010)
)

# Example usage 3: Multiple lists of lists (e.g., multiple ice core regions)
summary_stats_ice_all <- plot_age_range_violins(
  data_list = list(
    boreal_west = ice_core.boreal_west_list,  # List of 2 ice cores
    boreal_east = ice_core.boreal_east_list,   # List of 10 ice cores
    temp_west = ice_core.temperate_west_list
  ),
  age_ranges = age_ranges,
  age_var = "year_ce",
  value_vars = c("rBC_deposition"),
  normalize = TRUE,
  ylab = "Normalized rBC deposition (0-1)",
  main_title = "Ice Core Data",
  save_path = output_dir,
  fname_prefix = "ice_core_violin",
  fig_width = 2.5,
  fig_height = 5,
  violin_color = "#1473CC",
  violin_alpha = 0.5,
  add_jitter = TRUE,
  jitter_color = "black",
  jitter_alpha = 0.4,
  jitter_size = 0.8,
  jitter_amount = 0.15,
  ylims = c(0, 1)
)

print(summary_stats_ice_all)


summary_stats_scar <- plot_age_range_violins(
  data_list = list(boreal_west = scar.boreal_west,  # Each is a single dataframe
                   boreal_east = scar.boreal_east,
                   temp_west = scar.temperate_west,
                   temp_east = scar.temperate_east),
  age_ranges = age_ranges,
  age_var = "year_ce",
  value_vars = c("percent_scarred"),
  normalize = TRUE,
  ylab = "Normalized % sites scarred (0-1)",
  main_title = "Fire Scar Data",
  save_path = output_dir,
  fname_prefix = "fire_scar_violin",
  fig_width = 2.5,
  fig_height = 5,
  violin_color = "#3C7C4F",
  violin_alpha = 0.5,
  add_jitter = TRUE,
  jitter_color = "darkgreen",
  jitter_alpha = 0.5,
  jitter_size = 1.0,
  jitter_amount = 0.12,
  ylims = c(0, 1)
)

print(summary_stats_scar)

summary_stats_char <- plot_age_range_violins(
  data_list = list(boreal_west = char.boreal_west,  # Each is a single dataframe
                   boreal_east = char.boreal_east,
                   temp_west = char.temp_west,
                   temp_east = char.temp_east),
  age_ranges = age_ranges,
  age_var = "year_ce",
  value_vars = c("fit_zscore"),
  normalize = TRUE,
  ylab = "Normalized charcoal influx z-score (0-1)",
  main_title = "Sedimentary Charcoal Data",
  save_path = output_dir,
  fname_prefix = "char_violin",
  fig_width = 2.5,
  fig_height = 5,
  violin_color = "#7B6868",
  violin_alpha = 0.5,
  add_jitter = TRUE,
  jitter_color = "black",
  jitter_alpha = 0.5,
  jitter_size = 1.0,
  jitter_amount = 0.12,
  ylims = c(0, 1)
)

print(summary_stats_char)




