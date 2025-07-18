#!/usr/bin/env Rscript

## FRIP_plot.R - Single Table Input Version
## Created by Tovah Markowitz
## June 19, 2020
## Most recent update: October 17, 2024
## Reorganized to accept a single FRiP table file input

# Load required libraries
suppressPackageStartupMessages({
  library(ggplot2)
  library(rjson)
  library(ComplexHeatmap)
  library(tidyr)
  library(circlize)
  library(reshape2)
  library(argparse)
})

options(error = function() { traceback(3) })

# Command line argument parsing
get_command_args <- function() {
  parser <- ArgumentParser(
    description = "Generate FRiP (Fraction of Reads in Peaks) analysis plots and tables from a single input table",
    epilog = "Example: Rscript FRIP_plot.R -t FRiP_data.txt -c config.json -o results/"
  )
  
  parser$add_argument(
    "-t", "--table", 
    type = "character", 
    required = TRUE,
    help = "Path to FRiP table file", 
    metavar = "file.txt"
  )
  
  parser$add_argument(
    "-v", "--verbose", 
    action = "store_true", 
    default = FALSE,
    help = "Print verbose output"
  )
  
  parser$add_argument(
    "-b", "--barplot", 
    type = "character",
    required = TRUE,
    help = "Bar plot output file",
    metavar = "barplot"
  )
  
  parser$add_argument(
    "--hm", 
    type = "character",
    required = TRUE,
    help = "Heatmap output file",
    metavar = "heatmap"
  )

  # parser$add_argument(
  #   "--summary", 
  #   type = "character",
  #   required = TRUE,
  #   help = "Summary output file",
  #   metavar = "summary"
  # )
  
  # parser$add_argument(
  #   "--stats", 
  #   type = "character",
  #   required = TRUE,
  #   help = "Statistics output file",
  #   metavar = "stats"
  # )
  
  parser$add_argument(
    "-s", "--scatter", 
    type = "character",
    required = TRUE,
    help = "Scatter output file",
    metavar = "scatter"
  )

  parser$add_argument(
    "-c", "--config", 
    type = "character",
    required = TRUE,
    help = "Config input file",
    metavar = "config"
  )
  
  # Parse arguments
  args <- parser$parse_args()
  
  # Validate file existence
  if (!file.exists(args$table)) {
    stop(paste("Error: Table file does not exist:", args$table), call. = FALSE)
  }
  
  return(args)
}

# Data processing functions
load_frip_table <- function(table_file, verbose = FALSE) {
  #
  # Load FRiP table file and validate structure
  #
  if (verbose) {
    cat("Loading FRiP table file:", table_file, "\n")
  }
  
  data <- tryCatch({
    read.table(table_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  }, error = function(e) {
    stop(paste("Error reading table file:", e$message), call. = FALSE)
  })
  
  # Check for required columns
  required_cols <- c("bedsample", "bamsample", "FRiP", "bedtool")
  missing_cols <- required_cols[!required_cols %in% colnames(data)]
  if (length(missing_cols) > 0) {
    stop(paste("Missing required columns in data:", paste(missing_cols, collapse = ", ")), 
         call. = FALSE)
  }
  
  # Check for optional columns
  optional_cols <- c("n_basesM", "n_bases")
  has_optional <- optional_cols[optional_cols %in% colnames(data)]
  
  if (verbose) {
    cat("Loaded data:", nrow(data), "rows,", ncol(data), "columns\n")
    cat("Required columns found:", paste(required_cols, collapse = ", "), "\n")
    if (length(has_optional) > 0) {
      cat("Optional columns found:", paste(has_optional, collapse = ", "), "\n")
    }
    cat("Unique bedtools:", paste(unique(data$bedtool), collapse = ", "), "\n")
    cat("Unique bed samples:", length(unique(data$bedsample)), "\n")
    cat("Unique bam samples:", length(unique(data$bamsample)), "\n")
  }
  
  return(data)
}

# Plotting functions
create_barplots <- function(data, group_name, output_file, verbose = FALSE) {
  #
  # Create bar plots for FRiP values by group
  #
  if (verbose) {
    cat("Creating bar plots for group:", group_name, "\n")
  }
  
  p <- ggplot(data, aes(x = bamsample, y = FRiP, fill = bedsample)) +
    geom_bar(position = "dodge", stat = "identity") +
    facet_wrap(. ~ bedtool) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = -15, hjust = 0)) +
    labs(title = paste("FRiP Analysis:", group_name), 
         x = "BAM file", 
         y = "Fraction of Reads in Peaks (FRiP)", 
         fill = "Peak file")

  pdf(output_file, width = 10, height = 6)
  print(p)
  dev.off()
  
  if (verbose) {
    cat("Bar plot saved to:", output_file, "\n")
  }
}

create_scatterplots <- function(data, group_name, output_file, verbose = FALSE) {
  #
  # Create scatter plots for FRiP vs number of bases
  #
  if (verbose) {
    cat("Creating scatter plots for group:", group_name, "\n")
  }
  
  # Check if n_basesM column exists
  if (!"n_basesM" %in% colnames(data)) {
    if ("n_bases" %in% colnames(data)) {
      # Convert n_bases to millions
      data$n_basesM <- data$n_bases / 1000000
      if (verbose) {
        cat("Converted n_bases to n_basesM (millions)\n")
      }
    } else {
      if (verbose) {
        cat("No n_basesM or n_bases column found, skipping scatter plots\n")
      }
      return(NULL)
    }
  }
  
  p <- ggplot(data, aes(x = n_basesM, y = FRiP, shape = bedsample, color = bedtool)) +
    geom_point(size = 2.5) +
    facet_wrap(. ~ bamsample) +
    theme_bw() + 
    scale_x_continuous(trans = "log10") +
    labs(title = paste("FRiP vs Peak Size:", group_name), 
         x = "Number of Bases in Peaks (M)", 
         y = "Fraction of Reads in Peaks (FRiP)",
         shape = "Peak file", 
         color = "Peak calling tool")
  
  # Try to add log ticks, fallback to basic plot if error
  pdf(output_file, width = 12, height = 8)
  tryCatch({
    p_with_ticks <- p + annotation_logticks(sides = "b")
    print(p_with_ticks)
  }, error = function(e) {
    if (verbose) {
      cat("Could not add log ticks, using basic plot\n")
    }
    print(p)
  })
  dev.off()
  
  if (verbose) {
    cat("Scatter plot saved to:", output_file, "\n")
  }
}

create_tool_barplots <- function(data, output_file, verbose = FALSE) {
  #
  # Create bar plots for each bedtool separately
  #
  if (verbose) {
    cat("Creating tool-specific bar plots\n")
  }
  
  bedtools <- unique(data$bedtool)
  
  for (tool in bedtools) {
    tool_data <- data[data$bedtool == tool, ]
    
    if (nrow(tool_data) == 0) {
      if (verbose) {
        cat("No data for tool:", tool, "\n")
      }
      next
    }
    
    p <- ggplot(tool_data, aes(x = bamsample, y = FRiP, fill = groupInfo)) +
      geom_bar(position = "dodge", stat = "identity") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = -15, hjust = 0)) +
      labs(title = paste("FRiP Analysis by Group -", tool),
           x = "BAM file", 
           y = "Fraction of Reads in Peaks (FRiP)", 
           fill = "Group")

    pdf(output_file, width = 10, height = 6)
    print(p)
    dev.off()
    
    if (verbose) {
      cat("Tool bar plot saved to:", output_file, "\n")
    }
  }
}

create_overview_scatterplot <- function(data, output_file, verbose = FALSE) {
  #
  # Create overview scatter plot for all samples
  #
  if (verbose) {
    cat("Creating overview scatter plot\n")
  }
  
  # Check if n_basesM column exists
  if (!"n_basesM" %in% colnames(data)) {
    if ("n_bases" %in% colnames(data)) {
      # Convert n_bases to millions
      data$n_basesM <- data$n_bases / 1000000
      if (verbose) {
        cat("Converted n_bases to n_basesM (millions)\n")
      }
    } else {
      if (verbose) {
        cat("No n_basesM or n_bases column found, skipping overview scatter plot\n")
      }
      return(NULL)
    }
  }
  
  p <- ggplot(data, aes(x = n_basesM, y = FRiP, shape = bedtool, color = groupInfo)) +
    geom_point(size = 2.5) +
    theme_bw() + 
    scale_x_continuous(trans = "log10") +
    annotation_logticks(sides = "b") +
    labs(title = "FRiP Analysis Overview - All Samples", 
         x = "Number of Bases in Peaks (M)", 
         y = "Fraction of Reads in Peaks (FRiP)",
         shape = "Peak calling tool", 
         color = "Group")
  
  pdf(output_file, width = 10, height = 6)
  print(p)
  dev.off()
  
  if (verbose) {
    cat("Overview scatter plot saved to:", output_file, "\n")
  }
}

process_config_json <- function(config_file, verbose = FALSE) {
  #
  # Process config.json to extract group information and sample associations
  # 
  if (verbose) {
    cat("Processing config file:", config_file, "\n")
  }
  
  json_data <- tryCatch({
    fromJSON(file = config_file)
  }, error = function(e) {
    stop(paste("Error reading config file:", e$message), call. = FALSE)
  })
  
  if (is.null(json_data$project$groups)) {
    stop("No groups found in config file", call. = FALSE)
  }
  
  groups_info <- json_data$project$groups
  
  # Handle inputs if they exist
  if (!is.null(json_data$project$peaks$inputs)) {
    inputs <- unlist(json_data$project$peaks$inputs)
    
    # Add input samples to each group
    for (i in 1:length(groups_info)) {
      group_inputs <- unique(unlist(inputs[names(inputs) %in% groups_info[[i]]]))
      if (length(group_inputs) > 1) {
        groups_info[[i]] <- c(groups_info[[i]], as.character(group_inputs))
      } else if (length(group_inputs) == 1 && group_inputs != "") {
        groups_info[[i]] <- c(groups_info[[i]], as.character(group_inputs))
      }
    }
  }
  
  if (verbose) {
    cat("Found", length(groups_info), "groups:\n")
    for (i in 1:length(groups_info)) {
      cat("  ", names(groups_info)[i], ":", length(groups_info[[i]]), "samples\n")
    }
  }
  
  return(groups_info)
}

create_heatmap <- function(data, peakcaller, output_file, verbose = FALSE) {
  #
  # Create heatmap for FRiP values
  #
  if (verbose) {
    cat("Creating heatmap for:", peakcaller, "\n")
  }
  
  # Filter data for specific peak caller
  plot_data <- data[data$bedtool == peakcaller, c('bedsample', 'bamsample', 'FRiP')]
  
  if (nrow(plot_data) == 0) {
    if (verbose) {
      cat("No data found for peak caller:", peakcaller, "\n")
    }
    return(NULL)
  }
  
  # Reshape data for heatmap
  plot_data <- tryCatch({
    tidyr::pivot_wider(plot_data, names_from = bamsample, values_from = FRiP)
  }, error = function(e) {
    warning(paste("Error reshaping data for heatmap:", e$message))
    return(NULL)
  })
  
  if (is.null(plot_data)) {
    return(NULL)
  }
  
  plot_data <- data.frame(plot_data, check.names = FALSE)
  rownames(plot_data) <- plot_data$bedsample
  plot_data <- as.matrix(plot_data[, -1])
  
  # Handle case where matrix is empty or has issues
  if (nrow(plot_data) == 0 || ncol(plot_data) == 0) {
    if (verbose) {
      cat("Empty matrix for peak caller:", peakcaller, "\n")
    }
    return(NULL)
  }
  
  # Create heatmap
  max_val <- max(plot_data, na.rm = TRUE)
  if (is.infinite(max_val) || is.na(max_val)) {
    max_val <- 1
  }
  
  ht <- ComplexHeatmap::Heatmap(
    plot_data, 
    na_col = "grey",
    col = circlize::colorRamp2(c(0, 0.1, max_val), c('red', 'orange', 'blue')), 
    heatmap_legend_param = list(title = 'FRiP'), 
    row_title = 'Peak samples', 
    column_title = 'BAM samples',
    name = paste("FRiP -", peakcaller)
  )
  
  pdf(output_file, width = 8, height = 6)
  tryCatch({
    print(ht)
  }, error = function(e) {
    warning(paste("Error creating heatmap:", e$message))
  })
  dev.off()
  
  if (verbose) {
    cat("Heatmap saved to:", output_file, "\n")
  }
}

create_summary_stats <- function(data, stats_file, summary_file, verbose = FALSE) {
  #
  # Create summary statistics table
  #
  if (verbose) {
    cat("Creating summary statistics\n")
  }
  
  # Overall statistics
  summary_stats <- data.frame(
    Metric = c("Total Samples", "Unique BED Samples", "Unique BAM Samples", 
               "Peak Calling Tools", "Mean FRiP", "Median FRiP", "Min FRiP", 
               "Max FRiP"),
    Value = c(
      nrow(data),
      length(unique(data$bedsample)),
      length(unique(data$bamsample)),
      length(unique(data$bedtool)),
      round(mean(data$FRiP, na.rm = TRUE), 4),
      round(median(data$FRiP, na.rm = TRUE), 4),
      round(min(data$FRiP, na.rm = TRUE), 4),
      round(max(data$FRiP, na.rm = TRUE), 4)
    )
  )
  
  # Per-tool statistics
  tool_stats <- aggregate(FRiP ~ bedtool, data, function(x) {
    c(Count = length(x), Mean = round(mean(x, na.rm = TRUE), 4), 
      Median = round(median(x, na.rm = TRUE), 4), 
      Min = round(min(x, na.rm = TRUE), 4), 
      Max = round(max(x, na.rm = TRUE), 4))
  })
  
  # Save summary statistics
  write.table(summary_stats, summary_file, quote = FALSE, row.names = FALSE, sep = "\t")
  write.table(tool_stats, stats_file, quote = FALSE, row.names = FALSE, sep = "\t")
  
  if (verbose) {
    cat("Summary statistics saved to:", summary_file, "\n")
    cat("Tool statistics saved to:", stats_file, "\n")
  }
}

# Main function
main <- function() {
  # Get command line arguments
  args <- get_command_args()
  
  if (args$verbose) {
    cat("Starting FRiP analysis...\n")
    cat("Input table:", args$table, "\n")
    cat("Config file:", args$config, "\n")
  }
  
  # Load FRiP table
  data <- load_frip_table(args$table, args$verbose)
  
  # Create summary statistics
  # create_summary_stats(data, args$stats, args$summary, args$verbose)

  tool <- basename(args$table)
  tool <- sub("\\..*", "", tool)

  create_heatmap(data, tool, args$hm, args$verbose)
  
  # Process config file and create group plots
  group_list <- process_config_json(args$config, args$verbose)
  
  for (i in 1:length(group_list)) {
    group <- group_list[[i]]
    group_name <- names(group_list)[i]
    
    # Filter data for this group
    group_data <- data[
      (data$bedsample %in% group) & (data$bamsample %in% group), 
    ]
    
    if (nrow(group_data) > 0) {
      create_barplots(group_data, group_name, args$barplot, args$verbose)
      create_scatterplots(group_data, group_name, args$scatterplot, args$verbose)
    } else {
      if (args$verbose) {
        cat("No data found for group:", group_name, "\n")
      }
    }
  }
  
  # Create self-comparison plots (where bedsample == bamsample)
  self_data <- data[data$bedsample == data$bamsample, ]
  
  if (nrow(self_data) > 0) {
    # Add group information
    group_info <- reshape2::melt(group_list)
    names(group_info) <- c("bamsample", "groupInfo")
    self_data_merged <- merge(self_data, group_info, by = "bamsample", all.x = TRUE)
    
    # Handle samples not in any group
    self_data_merged$groupInfo[is.na(self_data_merged$groupInfo)] <- "Ungrouped"
    
    create_tool_barplots(self_data_merged, args$barplot, args$verbose)
    create_overview_scatterplot(self_data_merged, args$scatter, args$verbose)
    
    # Write self-comparison data table
    # write.table(self_data_merged, self_output, quote = FALSE, row.names = FALSE, sep = "\t")
    # if (args$verbose) {
    #   cat("Self-comparison table saved to:", self_output, "\n")
    # }

  } else {
    if (args$verbose) {
      cat("No self-comparison data found\n")
    }
  }
  
  if (args$verbose) {
    cat("FRiP analysis completed successfully!\n")
  }
}

# Run main function if script is executed directly
if (!interactive()) {
  main()
}