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
saveOn <- TRUE
smooth_scar_data <- TRUE

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

# ----------------------------- plot_timeseries --------------------------------- #

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

          ci_color <- adjustcolor(color_raw, alpha.f = pltParams[["alpha"]])
          polygon(c(age, rev(age)),
                  c(upper, rev(lower)),
                  col = ci_color, border = NA)
        }
      }

      # Plot the raw line with alpha
      raw_color_with_alpha <- adjustcolor(color_raw, alpha.f = raw_alpha)
      lines(age, value, col = raw_color_with_alpha, lwd = pltParams[["linewidths"]][1], lty = 1)

      # Plot first smooth curve if requested
      if(plot_smooth1 && !is.null(smooth1_var)) {
        if(smooth1_var %in% names(data)) {
          smooth1_value <- data[[smooth1_var]]
          smooth1_color_with_alpha <- adjustcolor(color_smooth1, alpha.f = smooth_alpha)
          lines(age, smooth1_value, col = smooth1_color_with_alpha, lwd = pltParams[["linewidths"]][2], lty = 1)
        }
      }

      # Plot second smooth curve if requested
      if(plot_smooth2 && !is.null(smooth2_var)) {
        if(smooth2_var %in% names(data)) {
          smooth2_value <- data[[smooth2_var]]
          smooth2_color_with_alpha <- adjustcolor(color_smooth2, alpha.f = smooth_alpha)
          lines(age, smooth2_value, col = smooth2_color_with_alpha, lwd = pltParams[["linewidths"]][3], lty = 1)
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



# =============================== Load in data ================================ ####

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

scar.boreal_west <- read.csv(file.path(data_dir, "fire_scars_Boreal_NA_west.csv"), header = TRUE)
scar.boreal_east <- read.csv(file.path(data_dir, "fire_scars_Boreal_NA_east.csv"), header = TRUE)
scar.temperate_west <- read.csv(file.path(data_dir, "fire_scars_Temperate_NA_west.csv"), header = TRUE)
scar.temperate_east <- read.csv(file.path(data_dir, "fire_scars_Temperate_NA_east.csv"), header = TRUE)

scar.list <- list(boreal_west = scar.boreal_west,
                  boreal_east = scar.boreal_east,
                  temperate_west = scar.temperate_west,
                  temperate_east = scar.temperate_east)

# ======================= Create master list of lists ========================= ####

data.full.list <- list(
  ice.model = ice.model.list,
  ice.data = ice.data.list,
  scar = scar.list
)

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
    single_color = "#479EEB",
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

# ========================= Define custom functions =========================== ####



# ============================== Make the plot ================================ ####

# Add a toggle for using single color
use_single_color <- TRUE  # Set to TRUE to use single color per dataset

# Loop through each dataset type (ice.model, ice.data, scar)
for(dataset_type in names(data.full.list)) {
  
  cat("\n=== Processing", dataset_type, "===\n")
  
  # Get the list of regions/datasets for this type
  region_list <- data.full.list[[dataset_type]]
  
  # Get specific parameters for this dataset type
  specific_params <- pltParams.specific[[dataset_type]]
  
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


