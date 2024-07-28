save_plots <- function(qcmetrics,
                       outdir,
                       name,
                       background = NULL,
                       background_alpha = 0.3,
                       width = 7,
                       height = 5,
                       unit = "in",
                       filetype = "png",
                       dpi = 300,
                       colorby = "None",
                       colorbyQCmetric = "None",
                       show_outliers = TRUE,
                       show_which_outliers = c("all", "metric"),
                       selected_samples = c(),
                       probes_plot_cutoff = 0.01,
                       text_size = 13,
                       point_size = 2.5,
                       point_size_outlier = point_size*1.5,
                       outlier_col = "red",
                       alpha = 0.8,
                       theme = theme_classic(),
                       colorscheme = NULL

) {
  plots <- list()
  show_which_outliers <- match.arg(show_which_outliers)

  plots$MU_plot <- .MU_plot(
    qcmetrics,
    threshold = NULL,
    colorby = colorby,
    colorbyQCmetric = colorbyQCmetric,
    outliers = if(show_which_outliers == "metric") getSampleOutliers(qcmetrics) %>% dplyr::filter(MU_outlier) %$% Sample_Name else getSampleOutliers(qcmetrics)$Sample_Name,
    show_outliers = show_outliers,
    selected_samples = selected_samples,
    background = background,
    background_alpha = background_alpha,
    text_size = text_size,
    point_size = point_size,
    point_size_outlier = point_size * 1.5,
    outlier_col = outlier_col,
    alpha = alpha,
    theme = theme,
    colorscheme = colorscheme
  )

  plots$RG_plot <- .RG_plot(
    qcmetrics,
    threshold = NULL,
    colorby = colorby,
    colorbyQCmetric = colorbyQCmetric,
    outliers = if(show_which_outliers == "metric") getSampleOutliers(qcmetrics) %>% dplyr::filter(RG_ratio_outlier | GR_ratio_outlier) %$% Sample_Name else getSampleOutliers(qcmetrics)$Sample_Name,
    show_outliers = show_outliers,
    selected_samples = selected_samples,
    background = background,
    background_alpha = background_alpha,
    text_size = text_size,
    point_size = point_size,
    point_size_outlier = point_size_outlier,
    outlier_col = outlier_col,
    alpha = alpha,
    theme = theme,
    colorscheme = colorscheme
  )

  plots$OP_plot <- .OP_plot(
    qcmetrics,
    threshold = NULL,
    colorby = colorby,
    colorbyQCmetric = colorbyQCmetric,
    outliers = if(show_which_outliers == "metric") getSampleOutliers(qcmetrics) %>% dplyr::filter(OP_outlier) %$% Sample_Name else getSampleOutliers(qcmetrics)$Sample_Name,
    show_outliers = show_outliers,
    selected_samples = selected_samples,
    background = background,
    background_alpha = background_alpha,
    text_size = text_size,
    point_size = point_size,
    point_size_outlier = point_size_outlier,
    outlier_col = outlier_col,
    alpha = alpha,
    theme = theme,
    colorscheme = colorscheme
  )

  plots$HC_plot <- .HC_plot(
    qcmetrics,
    threshold = NULL,
    colorby = colorby,
    colorbyQCmetric = colorbyQCmetric,
    outliers = if(show_which_outliers == "metric") getSampleOutliers(qcmetrics) %>% dplyr::filter(HC_outlier) %$% Sample_Name else getSampleOutliers(qcmetrics)$Sample_Name,
    show_outliers = show_outliers,
    selected_samples = selected_samples,
    background = background,
    background_alpha = background_alpha,
    text_size = text_size,
    point_size = point_size,
    point_size_outlier = point_size_outlier,
    outlier_col = outlier_col,
    alpha = alpha,
    theme = theme,
    colorscheme = colorscheme
  )

  plots$bscon_plot <- .bscon_plot(
    qcmetrics,
    threshold = NULL,
    colorby = colorby,
    colorbyQCmetric = colorbyQCmetric,
    outliers = if(show_which_outliers == "metric") getSampleOutliers(qcmetrics) %>% dplyr::filter(bscon_outlier) %$% Sample_Name else getSampleOutliers(qcmetrics)$Sample_Name,
    show_outliers = show_outliers,
    selected_samples = selected_samples,
    background = NULL,
    text_size = text_size,
    point_size = point_size,
    point_size_outlier = point_size_outlier,
    outlier_col = outlier_col,
    alpha = alpha,
    theme = theme,
    colorscheme = colorscheme
  )

  plots$detP_samples_plot <- .detP_plot(
    qcmetrics,
    threshold = NULL,
    colorby = colorby,
    colorbyQCmetric = colorbyQCmetric,
    outliers = if(show_which_outliers == "metric") getSampleOutliers(qcmetrics) %>% dplyr::filter(detectionP_outlier) %$% Sample_Name else getSampleOutliers(qcmetrics)$Sample_Name,
    show_outliers = show_outliers,
    selected_samples = selected_samples,
    background = NULL,
    text_size = text_size,
    point_size = point_size,
    point_size_outlier = point_size_outlier,
    outlier_col = outlier_col,
    alpha = alpha,
    theme = theme,
    colorscheme = colorscheme
  )

  plots$beadNr_samples_plot <- .beadNr_plot(
    qcmetrics,
    threshold = NULL,
    colorby = colorby,
    colorbyQCmetric = colorbyQCmetric,
    outliers = if(show_which_outliers == "metric")  getSampleOutliers(qcmetrics) %>% dplyr::filter(beadNr_outlier) %$% Sample_Name else getSampleOutliers(qcmetrics)$Sample_Name,
    show_outliers = show_outliers,
    selected_samples = selected_samples,
    background = NULL,
    text_size = text_size,
    point_size = point_size,
    point_size_outlier = point_size_outlier,
    outlier_col = outlier_col,
    alpha = alpha,
    theme = theme,
    colorscheme = colorscheme
  )

  if("IBS_mean" %in% colnames(qcmetrics@Sample_Metrics)) {
    plots$IBS_plot <- .IBS_plot(
      qcmetrics,
      threshold = NULL,
      colorby = colorby,
      colorbyQCmetric = colorbyQCmetric,
      outliers = if(show_which_outliers == "metric") getSampleOutliers(qcmetrics) %>% dplyr::filter(IBS_outlier) %$% Sample_Name else getSampleOutliers(qcmetrics)$Sample_Name,
      show_outliers = show_outliers,
      selected_samples = selected_samples,
      background = background,
      text_size = text_size,
      point_size = point_size,
      point_size_outlier = point_size_outlier,
      outlier_col = outlier_col,
      alpha = alpha,
      theme = theme,
      colorscheme = colorscheme
    )
  }

  plots$XY_plot <- .XY_plot(
    qcmetrics,
    threshold = NULL,
    colorby = colorby,
    colorbyQCmetric = colorbyQCmetric,
    outliers = if(show_which_outliers == "metric") getSampleOutliers(qcmetrics) %>% dplyr::filter(sex_outlier) %$% Sample_Name else getSampleOutliers(qcmetrics)$Sample_Name,
    show_outliers = show_outliers,
    selected_samples = selected_samples,
    background = background,
    background_alpha = background_alpha,
    text_size = text_size,
    point_size = point_size,
    point_size_outlier = point_size_outlier,
    outlier_col = outlier_col,
    alpha = alpha,
    theme = theme,
    colorscheme = colorscheme
  )

  plots$XYdiff_plot <- .XYdiff_plot(
    qcmetrics,
    threshold = NULL,
    colorby = colorby,
    colorbyQCmetric = colorbyQCmetric,
    outliers = if(show_which_outliers == "metric") getSampleOutliers(qcmetrics) %>% dplyr::filter(sex_outlier) %$% Sample_Name else getSampleOutliers(qcmetrics)$Sample_Name,
    show_outliers = show_outliers,
    selected_samples = selected_samples,
    background = background,
    text_size = text_size,
    point_size = point_size,
    point_size_outlier = point_size_outlier,
    outlier_col = outlier_col,
    alpha = alpha,
    theme = theme,
    colorscheme = colorscheme
  )

  if("controlPCA_outlier" %in% colnames(getSampleOutliers(qcmetrics))) {
    controlPCA_outliers <- getSampleOutliers(qcmetrics) %>% dplyr::filter(controlPCA_outlier) %$% Sample_Name
  } else {
    controlPCA_outliers <- NULL
  }
  plots$controlPCA_plot <- .PCA_plot(
    qcmetrics,
    type = "control",
    threshold = NULL,
    colorby = colorby,
    colorbyQCmetric = colorbyQCmetric,
    outliers = if(show_which_outliers == "metric") controlPCA_outliers else getSampleOutliers(qcmetrics)$Sample_Name,
    show_outliers = show_outliers,
    selected_samples = selected_samples,
    background = NULL,
    text_size = text_size,
    point_size = point_size,
    point_size_outlier = point_size_outlier,
    outlier_col = outlier_col,
    alpha = alpha,
    theme = theme,
    colorscheme = colorscheme
  )

  if(!is.null(qcmetrics@PCs)) {
    if("PCA_outlier" %in% colnames(getSampleOutliers(qcmetrics))) {
      PCA_outliers <- getSampleOutliers(qcmetrics) %>% dplyr::filter(PCA_outliers) %$% Sample_Name
    } else {
      PCA_outliers <- NULL
    }

    plots$PCA_plot <- .PCA_plot(
      qcmetrics,
      type = "arraywide",
      threshold = NULL,
      colorby = colorby,
      colorbyQCmetric = colorbyQCmetric,
      outliers = if(show_which_outliers == "metric") PCA_outliers else getSampleOutliers(qcmetrics)$Sample_Name,
      show_outliers = show_outliers,
      selected_samples = selected_samples,
      background = NULL,
      text_size = text_size,
      point_size = point_size,
      point_size_outlier = point_size_outlier,
      outlier_col = outlier_col,
      alpha = alpha,
      theme = theme,
      colorscheme = colorscheme
    )
  }

  plots$detP_probes_plot <- .detP_probes_plot(
    qcmetrics,
    threshold = NULL,
    plot_cutoff = probes_plot_cutoff,
    outliers = if(show_which_outliers == "metric") getProbeOutliers(qcmetrics) %>% dplyr::filter(detP_Percentage) %$% Probe else getProbeOutliers(qcmetrics)$Probe,
    show_outliers = show_outliers,
    text_size = text_size,
    point_size = point_size,
    point_size_outlier = point_size_outlier,
    outlier_col = outlier_col,
    alpha = alpha,
    theme = theme,
    colorscheme = NULL
  )

  plots$beadNr_probes_plot <- .beadNr_probes_plot(
    qcmetrics,
    threshold = NULL,
    plot_cutoff = probes_plot_cutoff,
    outliers = if(show_which_outliers == "metric") getProbeOutliers(qcmetrics) %>% dplyr::filter(BeadNr_Percentage) %$% Probe else getProbeOutliers(qcmetrics)$Probe,
    show_outliers = show_outliers,
    text_size = text_size,
    point_size = point_size,
    point_size_outlier = point_size_outlier,
    outlier_col = outlier_col,
    alpha = alpha,
    theme = theme,
    colorscheme = NULL
  )

  if(!dir.exists(outdir)) {message(sprintf("Folder '%s' did not exist, created it.", outdir)) ; dir.create(outdir)}
  message(sprintf("Saving plots in %s/", outdir))
  for(plotname in names(plots)) {
    ggsave(plots[[plotname]],
           filename = sprintf("%s/%s_%s.%s", outdir, name, plotname, filetype),
           width = width,
           height = height,
           unit = unit,
           dpi = dpi
    )
  }
}

.MU_plot <- function(qcmetrics,
                    threshold = NULL,
                    colorby = "None",
                    colorbyQCmetric = "None",
                    outliers = c(),
                    show_outliers = FALSE,
                    selected_samples = c(),
                    background = NULL,
                    background_alpha = 0.3,
                    text_size = 13,
                    point_size = 2.5,
                    point_size_outlier = point_size * 1.5,
                    outlier_col = "red",
                    alpha = 0.8,
                    theme = theme_classic(),
                    colorscheme = NULL
) {
  plotdata <- qcmetrics@Sample_Metrics %>%
    dplyr::left_join(qcmetrics@samplesheet,by="Sample_Name") %>%
    dplyr::select(-predictedSex) %>%
    dplyr::left_join(getPMS(qcmetrics), by="Sample_Name") %>%
    dplyr::left_join(qcmetrics@sample_qc_outliers$sample_outliers[,!colnames(qcmetrics@sample_qc_outliers$sample_outliers) %in% setdiff(colnames(qcmetrics@samplesheet), "Sample_Name")], by="Sample_Name")
  names <- colnames(qcmetrics@sample_qc_outliers$sample_outliers)
  names <- names[names != "Sample_Name"]
  plotdata[names][is.na(plotdata[names])] <- FALSE

  threshold = if(is.null(threshold)) getThresholds(qcmetrics)$sample_qc_thresholds$MU else threshold

  if((colorby != "None" || colorbyQCmetric != "None") && is.null(colorscheme)) {
    if(colorbyQCmetric != "None") {
      if(is.factor(plotdata[[colorbyQCmetric]]) || is.character(plotdata[[colorbyQCmetric]]) || is.logical(plotdata[[colorbyQCmetric]] )) {
        colorscheme <- scale_color_viridis_d()
      } else {
        colorscheme <- scale_color_viridis_c()
      }
    } else {
      if(is.factor(plotdata[[colorby]]) || is.character(plotdata[[colorby]]) || is.logical(plotdata[[colorby]] )) {
        colorscheme <- scale_color_viridis_d()
      } else {
        colorscheme <- scale_color_viridis_c()
      }
    }
  } else if(is.null(colorscheme)) {
    colorscheme <- scale_color_viridis_c()
  }

  if(!is.null(background)) {
    background <- background@Sample_Metrics

    # Check if any background data is also in present in foreground data
    overlap <- sum(background$Sample_Name %in% plotdata$Sample_Name)
    if(overlap > 0) {
      print(sprintf("%s samples in the foreground data are present in the background data. They will be removed from the background.", overlap))
      background <- background %>% dplyr::filter(!(Sample_Name %in% plotdata$Sample_Name))
    }
    if(nrow(background) == 0) {
      # If there is total overlap between foreground and background, the background option is turned off
      background <- NULL
    } else {
      # Create hexabin backgrounds
      background <- ggplot(background, aes(x = Unmethylated, y = Methylated)) +
        geom_hex(alpha = background_alpha, show.legend=FALSE) +
        scale_fill_gradient2(low ="white", mid = "gray95", high = "#0072B5")
    }
  }

  if(show_outliers) {

    if(!is.null(background)) {
      g <- background + geom_point(data = plotdata,
                                   aes_string(
                                     color = if(colorbyQCmetric != "None") colorbyQCmetric else if(colorby != "None") colorby else NULL
                                   ),
                                   alpha = alpha,
                                   size = point_size) +
        geom_vline(xintercept = threshold, linetype = "dashed") +
        geom_hline(yintercept = threshold, linetype = "dashed") +
        xlim(0, NA) +
        ylim(0, NA) +
        theme +
        theme(
          text = element_text(size = text_size)
        ) +
        colorscheme
    } else {
      g <- ggplot(plotdata, aes_string(x = "Unmethylated", y = "Methylated",
                                       color = if(colorbyQCmetric != "None") colorbyQCmetric else if(colorby != "None") colorby else NULL)) +
        geom_point(alpha = alpha, size = point_size) +
        geom_vline(xintercept = threshold, linetype = "dashed") +
        geom_hline(yintercept = threshold, linetype = "dashed") +
        xlim(0, NA) +
        ylim(0, NA) +
        theme +
        theme(
          text = element_text(size = text_size)
        ) +
        colorscheme
    }

    if(length(selected_samples) > 0) {
      g <- g + geom_point(data = plotdata %>% dplyr::filter(Sample_Name %in% selected_samples), shape = 8, size = point_size_outlier, color = outlier_col)
    } else if(show_outliers && length(outliers) > 0) {
      g <- g + geom_point(data = plotdata %>% dplyr::filter(Sample_Name %in% outliers), shape = 8, size = point_size_outlier, color = outlier_col)
    }
  } else {

    if(!is.null(background)) {
      g <- background + geom_point(data = plotdata, aes_string(
        color = if(colorbyQCmetric != "None") colorbyQCmetric else if(colorby != "None") colorby else NULL),
        size = point_size,
        alpha = alpha
      ) +
        xlim(0, NA) +
        ylim(0, NA) +
        theme +
        theme(
          text = element_text(size = text_size)
        ) +
        colorscheme
    } else {
      g <- ggplot(plotdata, aes_string(x = "Unmethylated", y = "Methylated",
                                       color = if(colorbyQCmetric != "None") colorbyQCmetric else if(colorby != "None") colorby else NULL)) +
        geom_point(alpha=alpha, size = point_size) +
        xlim(0, NA) +
        ylim(0, NA) +
        theme +
        theme(
          text = element_text(size = text_size)
        ) +
        colorscheme

    }

    if(length(selected_samples) > 0) {
      g <- g + geom_point(data = plotdata %>% dplyr::filter(Sample_Name %in% selected_samples), shape = 8, size = point_size_outlier, color = outlier_col)
    }
  }
  g
}


.RG_plot <- function(qcmetrics,
                    threshold = NULL,
                    colorby = "None",
                    colorbyQCmetric = "None",
                    outliers = c(),
                    show_outliers = FALSE,
                    selected_samples = c(),
                    background = NULL,
                    background_alpha = 0.3,
                    text_size = 13,
                    point_size = 2.5,
                    point_size_outlier = point_size * 1.5,
                    outlier_col = "red",
                    alpha = 0.8,
                    theme = theme_classic(),
                    colorscheme = NULL
) {
  plotdata <- qcmetrics@Sample_Metrics %>%
    dplyr::left_join(qcmetrics@samplesheet,by="Sample_Name") %>%
    dplyr::select(-predictedSex) %>%
    dplyr::left_join(getPMS(qcmetrics), by="Sample_Name") %>%
    dplyr::left_join(qcmetrics@sample_qc_outliers$sample_outliers[,!colnames(qcmetrics@sample_qc_outliers$sample_outliers) %in% setdiff(colnames(qcmetrics@samplesheet), "Sample_Name")], by="Sample_Name")
  names <- colnames(qcmetrics@sample_qc_outliers$sample_outliers)
  names <- names[names != "Sample_Name"]
  plotdata[names][is.na(plotdata[names])] <- FALSE

  threshold = if(is.null(threshold)) c(getThresholds(qcmetrics)$sample_qc_thresholds$RG_ratio, getThresholds(qcmetrics)$sample_qc_thresholds$GR_ratio) else threshold

  if((colorby != "None" || colorbyQCmetric != "None") && is.null(colorscheme)) {
    if(colorbyQCmetric != "None") {
      if(is.factor(plotdata[[colorbyQCmetric]]) || is.character(plotdata[[colorbyQCmetric]]) || is.logical(plotdata[[colorbyQCmetric]] )) {
        colorscheme <- scale_color_viridis_d()
      } else {
        colorscheme <- scale_color_viridis_c()
      }
    } else {
      if(is.factor(plotdata[[colorby]]) || is.character(plotdata[[colorby]]) || is.logical(plotdata[[colorby]] )) {
        colorscheme <- scale_color_viridis_d()
      } else {
        colorscheme <- scale_color_viridis_c()
      }
    }
  } else if(is.null(colorscheme)) {
    colorscheme <- scale_color_viridis_c()
  }

  if(!is.null(background)) {
    background <- background@Sample_Metrics

    # Check if any background data is also in present in foreground data
    overlap <- sum(background$Sample_Name %in% plotdata$Sample_Name)
    if(overlap > 0) {
      print(sprintf("%s samples in the foreground data are present in the background data. They will be removed from the background.", overlap))
      background <- background %>% dplyr::filter(!(Sample_Name %in% plotdata$Sample_Name))
    }
    if(nrow(background) == 0) {
      # If there is total overlap between foreground and background, the background option is turned off
      background <- NULL
    } else {
      # Create hexabin backgrounds
      background <- ggplot(background, aes(x = Red, y = Green)) +
        geom_hex(alpha = background_alpha, show.legend=FALSE) +
        scale_fill_gradient2(low ="white", mid = "gray95", high = "#0072B5")
    }
  }
  if(show_outliers) {

    if(!is.null(background)) {
      g <- background + geom_point(data = plotdata,
                                   aes_string(
                                     color = if(colorbyQCmetric != "None") colorbyQCmetric else if(colorby != "None") colorby else NULL
                                   ),
                                   alpha = alpha,
                                   size = point_size) +
        geom_abline(slope = 1/threshold[1], linetype = "dashed") +
        geom_abline(slope = threshold[2], linetype = "dashed") +
        xlim(0, NA) +
        ylim(0, NA) +
        theme +
        theme(
          text = element_text(size = text_size)
        ) +
        colorscheme
    } else {
      g <- ggplot(plotdata, aes_string(x = "Red", y = "Green",
                                       color = if(colorbyQCmetric != "None") colorbyQCmetric else if(colorby != "None") colorby else NULL)) +
        geom_point(alpha = alpha, size = point_size) +
        geom_abline(slope = 1/threshold[1], linetype = "dashed") +
        geom_abline(slope = threshold[2], linetype = "dashed") +
        xlim(0, NA) +
        ylim(0, NA) +
        theme +
        theme(
          text = element_text(size = text_size)
        ) +
        colorscheme
    }

    if(length(selected_samples) > 0) {
      g <- g + geom_point(data = plotdata %>% dplyr::filter(Sample_Name %in% selected_samples), shape = 8, size = point_size_outlier, color = outlier_col)
    } else if(show_outliers && length(outliers) > 0) {
      g <- g + geom_point(data = plotdata %>% dplyr::filter(Sample_Name %in% outliers), shape = 8, size = point_size_outlier, color = outlier_col)
    }
  } else {

    if(!is.null(background)) {
      g <- background + geom_point(data = plotdata, aes_string(
        color = if(colorbyQCmetric != "None") colorbyQCmetric else if(colorby != "None") colorby else NULL),
        size = point_size,
        alpha = alpha
      ) +
        xlim(0, NA) +
        ylim(0, NA) +
        theme +
        theme(
          text = element_text(size = text_size)
        ) +
        colorscheme
    } else {
      g <- ggplot(plotdata, aes_string(x = "Red", y = "Green",
                                       color = if(colorbyQCmetric != "None") colorbyQCmetric else if(colorby != "None") colorby else NULL)) +
        geom_point(alpha=alpha, size = point_size) +
        xlim(0, NA) +
        ylim(0, NA) +
        theme +
        theme(
          text = element_text(size = text_size)
        ) +
        colorscheme
    }

    if(length(selected_samples) > 0) {
      g <- g + geom_point(data = plotdata %>% dplyr::filter(Sample_Name %in% selected_samples), shape = 8, size = point_size_outlier, color = outlier_col)
    }
  }
  g
}

.OP_plot <- function(qcmetrics,
                    threshold = NULL,
                    colorby = "None",
                    colorbyQCmetric = "None",
                    outliers = c(),
                    show_outliers = FALSE,
                    selected_samples = c(),
                    background = NULL,
                    background_alpha = 0.3,
                    text_size = 13,
                    point_size = 2.5,
                    point_size_outlier = point_size * 1.5,
                    outlier_col = "red",
                    alpha = 0.8,
                    theme = theme_classic(),
                    colorscheme = NULL
) {
  plotdata <- qcmetrics@Sample_Metrics %>%
    dplyr::left_join(qcmetrics@samplesheet,by="Sample_Name") %>%
    dplyr::select(-predictedSex) %>%
    dplyr::left_join(getPMS(qcmetrics), by="Sample_Name") %>%
    dplyr::left_join(qcmetrics@sample_qc_outliers$sample_outliers[,!colnames(qcmetrics@sample_qc_outliers$sample_outliers) %in% setdiff(colnames(qcmetrics@samplesheet), "Sample_Name")], by="Sample_Name")
  names <- colnames(qcmetrics@sample_qc_outliers$sample_outliers)
  names <- names[names != "Sample_Name"]
  plotdata[names][is.na(plotdata[names])] <- FALSE


  threshold = if(is.null(threshold)) getThresholds(qcmetrics)$sample_qc_thresholds$OP else threshold

  if((colorby != "None" || colorbyQCmetric != "None") && is.null(colorscheme)) {
    if(colorbyQCmetric != "None") {
      if(is.factor(plotdata[[colorbyQCmetric]]) || is.character(plotdata[[colorbyQCmetric]]) || is.logical(plotdata[[colorbyQCmetric]] )) {
        colorscheme <- scale_color_viridis_d()
      } else {
        colorscheme <- scale_color_viridis_c()
      }
    } else {
      if(is.factor(plotdata[[colorby]]) || is.character(plotdata[[colorby]]) || is.logical(plotdata[[colorby]] )) {
        colorscheme <- scale_color_viridis_d()
      } else {
        colorscheme <- scale_color_viridis_c()
      }
    }
  } else if(is.null(colorscheme)) {
    colorscheme <- scale_color_viridis_c()
  }

  if(!is.null(background)) {
    background <- background@Sample_Metrics

    # Check if any background data is also in present in foreground data
    overlap <- sum(background$Sample_Name %in% plotdata$Sample_Name)
    if(overlap > 0) {
      print(sprintf("%s samples in the foreground data are present in the background data. They will be removed from the background.", overlap))
      background <- background %>% dplyr::filter(!(Sample_Name %in% plotdata$Sample_Name))
    }
    if(nrow(background) == 0) {
      # If there is total overlap between foreground and background, the background option is turned off
      background <- NULL
    } else {
      # Create hexabin backgrounds
      background <- ggplot(background, aes(x = OP_x, y = OP_y)) +
        geom_hex(alpha = background_alpha, show.legend=FALSE) +
        scale_fill_gradient2(low ="white", mid = "gray95", high = "#0072B5")
    }
  }
  if(show_outliers) {

    if(!is.null(background)) {
      g <- background + geom_point(data = plotdata,
                                   aes_string(
                                     color = if(colorbyQCmetric != "None") colorbyQCmetric else if(colorby != "None") colorby else NULL
                                   ),
                                   alpha = alpha,
                                   size = point_size) +
        geom_vline(xintercept = threshold, linetype = "dashed") +
        theme +
        theme(
          text = element_text(size = text_size)
        ) +
        colorscheme
    } else {
      g <- ggplot(plotdata, aes_string(x = "OP_x", y = "OP_y",
                                       color = if(colorbyQCmetric != "None") colorbyQCmetric else if(colorby != "None") colorby else NULL)) +
        geom_point(alpha = alpha, size = point_size) +
        geom_vline(xintercept = threshold, linetype = "dashed") +
        theme +
        theme(
          text = element_text(size = text_size)
        ) +
        colorscheme
    }

    if(length(selected_samples) > 0) {
      g <- g + geom_point(data = plotdata %>% dplyr::filter(Sample_Name %in% selected_samples), shape = 8, size = point_size_outlier, color = outlier_col)
    } else if(show_outliers && length(outliers) > 0) {
      g <- g + geom_point(data = plotdata %>% dplyr::filter(Sample_Name %in% outliers), shape = 8, size = point_size_outlier, color = outlier_col)
    }
  } else {

    if(!is.null(background)) {
      g <- background + geom_point(data = plotdata, aes_string(
        color = if(colorbyQCmetric != "None") colorbyQCmetric else if(colorby != "None") colorby else NULL),
        size = point_size,
        alpha = alpha
      ) +
        theme +
        theme(
          text = element_text(size = text_size)
        ) +
        colorscheme
    } else {
      g <- ggplot(plotdata, aes_string(x = "OP_x", y = "OP_y",
                                       color = if(colorbyQCmetric != "None") colorbyQCmetric else if(colorby != "None") colorby else NULL)) +
        geom_point(alpha=alpha, size = point_size) +
        theme +
        theme(
          text = element_text(size = text_size)
        ) +
        colorscheme
    }

    if(length(selected_samples) > 0) {
      g <- g + geom_point(data = plotdata %>% dplyr::filter(Sample_Name %in% selected_samples), shape = 8, size = point_size_outlier, color = outlier_col)
    }
  }
  g
}

.HC_plot <- function(qcmetrics,
                    threshold = NULL,
                    colorby = "None",
                    colorbyQCmetric = "None",
                    outliers = c(),
                    show_outliers = FALSE,
                    selected_samples = c(),
                    background = NULL,
                    background_alpha = 0.3,
                    text_size = 13,
                    point_size = 2.5,
                    point_size_outlier = point_size * 1.5,
                    outlier_col = "red",
                    alpha = 0.8,
                    theme = theme_classic(),
                    colorscheme = NULL
) {
  plotdata <- qcmetrics@Sample_Metrics %>%
    dplyr::left_join(qcmetrics@samplesheet,by="Sample_Name") %>%
    dplyr::select(-predictedSex) %>%
    dplyr::left_join(getPMS(qcmetrics), by="Sample_Name") %>%
    dplyr::left_join(qcmetrics@sample_qc_outliers$sample_outliers[,!colnames(qcmetrics@sample_qc_outliers$sample_outliers) %in% setdiff(colnames(qcmetrics@samplesheet), "Sample_Name")], by="Sample_Name")
  names <- colnames(qcmetrics@sample_qc_outliers$sample_outliers)
  names <- names[names != "Sample_Name"]
  plotdata[names][is.na(plotdata[names])] <- FALSE

  threshold = if(is.null(threshold)) getThresholds(qcmetrics)$sample_qc_thresholds$HC else threshold

  if((colorby != "None" || colorbyQCmetric != "None") && is.null(colorscheme)) {
    if(colorbyQCmetric != "None") {
      if(is.factor(plotdata[[colorbyQCmetric]]) || is.character(plotdata[[colorbyQCmetric]]) || is.logical(plotdata[[colorbyQCmetric]] )) {
        colorscheme <- scale_color_viridis_d()
      } else {
        colorscheme <- scale_color_viridis_c()
      }
    } else {
      if(is.factor(plotdata[[colorby]]) || is.character(plotdata[[colorby]]) || is.logical(plotdata[[colorby]] )) {
        colorscheme <- scale_color_viridis_d()
      } else {
        colorscheme <- scale_color_viridis_c()
      }
    }
  } else if(is.null(colorscheme)) {
    colorscheme <- scale_color_viridis_c()
  }

  if(!is.null(background)) {
    background <- background@Sample_Metrics

    # Check if any background data is also in present in foreground data
    overlap <- sum(background$Sample_Name %in% plotdata$Sample_Name)
    if(overlap > 0) {
      print(sprintf("%s samples in the foreground data are present in the background data. They will be removed from the background.", overlap))
      background <- background %>% dplyr::filter(!(Sample_Name %in% plotdata$Sample_Name))
    }
    if(nrow(background) == 0) {
      # If there is total overlap between foreground and background, the background option is turned off
      background <- NULL
    } else {
      # Create hexabin backgrounds
      background <- ggplot(background, aes(x = HC_x, y = HC_y)) +
        geom_hex(alpha = background_alpha, show.legend=FALSE) +
        scale_fill_gradient2(low ="white", mid = "gray95", high = "#0072B5")
    }
  }
  if(show_outliers) {

    if(!is.null(background)) {
      g <- background + geom_point(data = plotdata,
                                   aes_string(
                                     color = if(colorbyQCmetric != "None") colorbyQCmetric else if(colorby != "None") colorby else NULL
                                   ),
                                   alpha = alpha,
                                   size = point_size) +
        geom_vline(xintercept = threshold, linetype = "dashed") +
        theme +
        theme(
          text = element_text(size = text_size)
        ) +
        colorscheme
    } else {
      g <- ggplot(plotdata, aes_string(x = "HC_x", y = "HC_y",
                                       color = if(colorbyQCmetric != "None") colorbyQCmetric else if(colorby != "None") colorby else NULL)) +
        geom_point(alpha = alpha, size = point_size) +
        geom_vline(xintercept = threshold, linetype = "dashed") +
        theme +
        theme(
          text = element_text(size = text_size)
        ) +
        colorscheme
    }

    if(length(selected_samples) > 0) {
      g <- g + geom_point(data = plotdata %>% dplyr::filter(Sample_Name %in% selected_samples), shape = 8, size = point_size_outlier, color = outlier_col)
    } else if(show_outliers && length(outliers) > 0) {
      g <- g + geom_point(data = plotdata %>% dplyr::filter(Sample_Name %in% outliers), shape = 8, size = point_size_outlier, color = outlier_col)
    }
  } else {

    if(!is.null(background)) {
      g <- background + geom_point(data = plotdata, aes_string(
        color = if(colorbyQCmetric != "None") colorbyQCmetric else if(colorby != "None") colorby else NULL),
        size = point_size,
        alpha = alpha
      ) +
        theme +
        theme(
          text = element_text(size = text_size)
        ) +
        colorscheme
    } else {
      g <- ggplot(plotdata, aes_string(x = "HC_x", y = "HC_y",
                                       color = if(colorbyQCmetric != "None") colorbyQCmetric else if(colorby != "None") colorby else NULL)) +
        geom_point(alpha=alpha, size = point_size) +
        theme +
        theme(
          text = element_text(size = text_size)
        ) +
        colorscheme
    }

    if(length(selected_samples) > 0) {
      g <- g + geom_point(data = plotdata %>% dplyr::filter(Sample_Name %in% selected_samples), shape = 8, size = point_size_outlier, color = outlier_col)
    }
  }
  g
}


.bscon_plot <- function(qcmetrics,
                       threshold = NULL,
                       colorby = "None",
                       colorbyQCmetric = "None",
                       outliers = c(),
                       show_outliers = FALSE,
                       selected_samples = c(),
                       background = NULL,
                       text_size = 13,
                       point_size = 2.5,
                       point_size_outlier = point_size * 1.5,
                       outlier_col = "red",
                       alpha = 0.8,
                       theme = theme_classic(),
                       colorscheme = NULL
) {
  plotdata <- qcmetrics@Sample_Metrics %>%
    dplyr::left_join(qcmetrics@samplesheet,by="Sample_Name") %>%
    dplyr::select(-predictedSex) %>%
    dplyr::left_join(getPMS(qcmetrics), by="Sample_Name") %>%
    dplyr::left_join(qcmetrics@sample_qc_outliers$sample_outliers[,!colnames(qcmetrics@sample_qc_outliers$sample_outliers) %in% setdiff(colnames(qcmetrics@samplesheet), "Sample_Name")], by="Sample_Name")
  names <- colnames(qcmetrics@sample_qc_outliers$sample_outliers)
  names <- names[names != "Sample_Name"]
  plotdata[names][is.na(plotdata[names])] <- FALSE

  threshold = if(is.null(threshold)) getThresholds(qcmetrics)$sample_qc_thresholds$bscon else threshold

  if((colorby != "None" || colorbyQCmetric != "None") && is.null(colorscheme)) {
    if(colorbyQCmetric != "None") {
      if(is.factor(plotdata[[colorbyQCmetric]]) || is.character(plotdata[[colorbyQCmetric]]) || is.logical(plotdata[[colorbyQCmetric]] )) {
        colorscheme <- scale_color_viridis_d()
      } else {
        colorscheme <- scale_color_viridis_c()
      }
    } else {
      if(is.factor(plotdata[[colorby]]) || is.character(plotdata[[colorby]]) || is.logical(plotdata[[colorby]] )) {
        colorscheme <- scale_color_viridis_d()
      } else {
        colorscheme <- scale_color_viridis_c()
      }
    }
  } else if(is.null(colorscheme)) {
    colorscheme <- scale_color_viridis_c()
  }

  if(!is.null(background)) {
    print("background not implemented for bscon")
  }

  if(show_outliers) {

    g <- ggplot(plotdata, aes_string(x = "bscon", y = "Sample_Name",
                                     color = if(colorbyQCmetric != "None") colorbyQCmetric else if(colorby != "None") colorby else NULL)) +
      geom_point(alpha = alpha, size = point_size) +
      geom_vline(xintercept = threshold, linetype = "dashed") +
      theme +
      theme(
        text = element_text(size = text_size),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank()
      ) +
      ylab("Samples") +
      xlab("bscon percentage") +
      colorscheme


    if(length(selected_samples) > 0) {
      g <- g + geom_point(data = plotdata %>% dplyr::filter(Sample_Name %in% selected_samples), shape = 8, size = point_size_outlier, color = outlier_col)
    } else if(show_outliers && length(outliers) > 0) {
      g <- g + geom_point(data = plotdata %>% dplyr::filter(Sample_Name %in% outliers), shape = 8, size = point_size_outlier, color = outlier_col)
    }
  } else {

    g <- ggplot(plotdata, aes_string(x = "bscon", y = "Sample_Name",
                                     color = if(colorbyQCmetric != "None") colorbyQCmetric else if(colorby != "None") colorby else NULL)) +
      geom_point(alpha=alpha, size = point_size) +
      theme +
      theme(
        text = element_text(size = text_size),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank()
      ) +
      ylab("Samples") +
      xlab("bscon percentage") +
      colorscheme

    if(length(selected_samples) > 0) {
      g <- g + geom_point(data = plotdata %>% dplyr::filter(Sample_Name %in% selected_samples), shape = 8, size = point_size_outlier, color = outlier_col)
    }
  }
  g
}

.detP_plot <- function(qcmetrics,
                      threshold = NULL,
                      colorby = "None",
                      colorbyQCmetric = "None",
                      outliers = c(),
                      show_outliers = FALSE,
                      selected_samples = c(),
                      background = NULL,
                      text_size = 13,
                      point_size = 2.5,
                      point_size_outlier = point_size * 1.5,
                      outlier_col = "red",
                      alpha = 0.8,
                      theme = theme_classic(),
                      colorscheme = NULL
) {
  plotdata <- qcmetrics@Sample_Metrics %>%
    dplyr::left_join(qcmetrics@samplesheet,by="Sample_Name") %>%
    dplyr::select(-predictedSex) %>%
    dplyr::left_join(getPMS(qcmetrics), by="Sample_Name") %>%
    dplyr::left_join(qcmetrics@sample_qc_outliers$sample_outliers[,!colnames(qcmetrics@sample_qc_outliers$sample_outliers) %in% setdiff(colnames(qcmetrics@samplesheet), "Sample_Name")], by="Sample_Name")
  names <- colnames(qcmetrics@sample_qc_outliers$sample_outliers)
  names <- names[names != "Sample_Name"]
  plotdata[names][is.na(plotdata[names])] <- FALSE

  threshold = if(is.null(threshold)) getThresholds(qcmetrics)$sample_qc_thresholds$detP else threshold
  threshold2 = getThresholds(qcmetrics)$detp_threshold

  if((colorby != "None" || colorbyQCmetric != "None") && is.null(colorscheme)) {
    if(colorbyQCmetric != "None") {
      if(is.factor(plotdata[[colorbyQCmetric]]) || is.character(plotdata[[colorbyQCmetric]]) || is.logical(plotdata[[colorbyQCmetric]] )) {
        colorscheme <- scale_color_viridis_d()
      } else {
        colorscheme <- scale_color_viridis_c()
      }
    } else {
      if(is.factor(plotdata[[colorby]]) || is.character(plotdata[[colorby]]) || is.logical(plotdata[[colorby]] )) {
        colorscheme <- scale_color_viridis_d()
      } else {
        colorscheme <- scale_color_viridis_c()
      }
    }
  } else if(is.null(colorscheme)) {
    colorscheme <- scale_color_viridis_c()
  }

  if(!is.null(background)) {
    print("background not implemented for detP")
  }

  if(show_outliers) {

    g <- ggplot(plotdata, aes_string(x = "Sample_Name", y = "detP_Percentage",
                                     color = if(colorbyQCmetric != "None") colorbyQCmetric else if(colorby != "None") colorby else NULL)) +
      geom_point(alpha = alpha, size = point_size) +
      geom_hline(yintercept = threshold, linetype = "dashed") +
      theme +
      theme(
        text = element_text(size = text_size),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()
      ) +
      ylab(sprintf("Proportion of probes with detP > %s", threshold)) +
      xlab("Samples") +
      colorscheme


    if(length(selected_samples) > 0) {
      g <- g + geom_point(data = plotdata %>% dplyr::filter(Sample_Name %in% selected_samples), shape = 8, size = point_size_outlier, color = outlier_col)
    } else if(show_outliers && length(outliers) > 0) {
      g <- g + geom_point(data = plotdata %>% dplyr::filter(Sample_Name %in% outliers), shape = 8, size = point_size_outlier, color = outlier_col)
    }
  } else {

    g <- ggplot(plotdata, aes_string(x = "Sample_Name", y = "detP_Percentage",
                                     color = if(colorbyQCmetric != "None") colorbyQCmetric else if(colorby != "None") colorby else NULL)) +
      geom_point(alpha=alpha, size = point_size) +
      theme +
      theme(
        text = element_text(size = text_size),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()
      ) +
      ylab(sprintf("Proportion of probes with detP > %s", threshold)) +
      xlab("Samples") +
      colorscheme


    if(length(selected_samples) > 0) {
      g <- g + geom_point(data = plotdata %>% dplyr::filter(Sample_Name %in% selected_samples), shape = 8, size = point_size_outlier, color = outlier_col)
    }
  }
  g
}

.beadNr_plot <- function(qcmetrics,
                        threshold = NULL,
                        colorby = "None",
                        colorbyQCmetric = "None",
                        outliers = c(),
                        show_outliers = FALSE,
                        selected_samples = c(),
                        background = NULL,
                        text_size = 13,
                        point_size = 2.5,
                        point_size_outlier = point_size * 1.5,
                        outlier_col = "red",
                        alpha = 0.8,
                        theme = theme_classic(),
                        colorscheme = NULL
) {
  plotdata <- qcmetrics@Sample_Metrics %>%
    dplyr::left_join(qcmetrics@samplesheet,by="Sample_Name") %>%
    dplyr::select(-predictedSex) %>%
    dplyr::left_join(getPMS(qcmetrics), by="Sample_Name") %>%
    dplyr::left_join(qcmetrics@sample_qc_outliers$sample_outliers[,!colnames(qcmetrics@sample_qc_outliers$sample_outliers) %in% setdiff(colnames(qcmetrics@samplesheet), "Sample_Name")], by="Sample_Name")
  names <- colnames(qcmetrics@sample_qc_outliers$sample_outliers)
  names <- names[names != "Sample_Name"]
  plotdata[names][is.na(plotdata[names])] <- FALSE

  threshold = if(is.null(threshold)) getThresholds(qcmetrics)$sample_qc_thresholds$beadNr else threshold
  threshold2 =  getThresholds(qcmetrics)$beadnr_threshold

  if((colorby != "None" || colorbyQCmetric != "None") && is.null(colorscheme)) {
    if(colorbyQCmetric != "None") {
      if(is.factor(plotdata[[colorbyQCmetric]]) || is.character(plotdata[[colorbyQCmetric]]) || is.logical(plotdata[[colorbyQCmetric]] )) {
        colorscheme <- scale_color_viridis_d()
      } else {
        colorscheme <- scale_color_viridis_c()
      }
    } else {
      if(is.factor(plotdata[[colorby]]) || is.character(plotdata[[colorby]]) || is.logical(plotdata[[colorby]] )) {
        colorscheme <- scale_color_viridis_d()
      } else {
        colorscheme <- scale_color_viridis_c()
      }
    }
  } else if(is.null(colorscheme)) {
    colorscheme <- scale_color_viridis_c()
  }

  if(!is.null(background)) {
    print("background not implemented for beadNr")
  }

  if(show_outliers) {

    g <- ggplot(plotdata, aes_string(x = "Sample_Name", y = "BeadNr_Percentage",
                                     color = if(colorbyQCmetric != "None") colorbyQCmetric else if(colorby != "None") colorby else NULL)) +
      geom_point(alpha = alpha, size = point_size) +
      geom_hline(yintercept = threshold, linetype = "dashed") +
      theme +
      theme(
        text = element_text(size = text_size),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()
      ) +
      ylab(sprintf("Proportion of probes with beadnr < %s", threshold2)) +
      xlab("Samples") +
      colorscheme

    if(length(selected_samples) > 0) {
      g <- g + geom_point(data = plotdata %>% dplyr::filter(Sample_Name %in% selected_samples), shape = 8, size = point_size_outlier, color = outlier_col)
    } else if(show_outliers && length(outliers) > 0) {
      g <- g + geom_point(data = plotdata %>% dplyr::filter(Sample_Name %in% outliers), shape = 8, size = point_size_outlier, color = outlier_col)
    }
  } else {

    g <- ggplot(plotdata, aes_string(x = "Sample_Name", y = "BeadNr_Percentage",
                                     color = if(colorbyQCmetric != "None") colorbyQCmetric else if(colorby != "None") colorby else NULL)) +
      geom_point(alpha=alpha, size = point_size) +
      theme +
      theme(
        text = element_text(size = text_size),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()
      ) +
      ylab(sprintf("Proportion of probes with beadnr < %s", threshold2)) +
      xlab("Samples") +
      colorscheme

    if(length(selected_samples) > 0) {
      g <- g + geom_point(data = plotdata %>% dplyr::filter(Sample_Name %in% selected_samples), shape = 8, size = point_size_outlier, color = outlier_col)
    }
  }
  g
}

.IBS_plot <- function(qcmetrics,
                     threshold = NULL,
                     colorby = "None",
                     colorbyQCmetric = "None",
                     outliers = c(),
                     show_outliers = FALSE,
                     selected_samples = c(),
                     background = NULL,
                     text_size = 13,
                     point_size = 2.5,
                     point_size_outlier = point_size * 1.5,
                     outlier_col = "red",
                     alpha = 0.8,
                     theme = theme_classic(),
                     colorscheme = NULL
) {
  plotdata <- qcmetrics@Sample_Metrics %>%
    dplyr::left_join(qcmetrics@samplesheet,by="Sample_Name") %>%
    dplyr::select(-predictedSex) %>%
    dplyr::left_join(getPMS(qcmetrics), by="Sample_Name") %>%
    dplyr::left_join(qcmetrics@sample_qc_outliers$sample_outliers[,!colnames(qcmetrics@sample_qc_outliers$sample_outliers) %in% setdiff(colnames(qcmetrics@samplesheet), "Sample_Name")], by="Sample_Name")
  names <- colnames(qcmetrics@sample_qc_outliers$sample_outliers)
  names <- names[names != "Sample_Name"]
  plotdata[names][is.na(plotdata[names])] <- FALSE

  threshold = if(is.null(threshold)) c(getThresholds(qcmetrics)$sample_qc_thresholds$IBS_mean, getThresholds(qcmetrics)$sample_qc_thresholds$IBS_var) else threshold

  if((colorby != "None" || colorbyQCmetric != "None") && is.null(colorscheme)) {
    if(colorbyQCmetric != "None") {
      if(is.factor(plotdata[[colorbyQCmetric]]) || is.character(plotdata[[colorbyQCmetric]]) || is.logical(plotdata[[colorbyQCmetric]] )) {
        colorscheme <- scale_color_viridis_d()
      } else {
        colorscheme <- scale_color_viridis_c()
      }
    } else {
      if(is.factor(plotdata[[colorby]]) || is.character(plotdata[[colorby]]) || is.logical(plotdata[[colorby]] )) {
        colorscheme <- scale_color_viridis_d()
      } else {
        colorscheme <- scale_color_viridis_c()
      }
    }
  } else if(is.null(colorscheme)) {
    colorscheme <- scale_color_viridis_c()
  }

  if(show_outliers) {

    g <- ggplot(plotdata, aes_string(x = "IBS_mean", y = "IBS_var",
                                     color = if(colorbyQCmetric != "None") colorbyQCmetric else if(colorby != "None") colorby else NULL)) +
      geom_point(alpha = alpha, size = point_size) +
      geom_vline(xintercept = threshold[1], linetype = "dashed") +
      geom_hline(yintercept = threshold[2], linetype = "dashed") +
      theme +
      theme(
        text = element_text(size = text_size)
      ) +
      colorscheme

    if(length(selected_samples) > 0) {
      g <- g + geom_point(data = plotdata %>% dplyr::filter(Sample_Name %in% selected_samples), shape = 8, size = point_size_outlier, color = outlier_col)
    } else if(show_outliers && length(outliers) > 0) {
      g <- g + geom_point(data = plotdata %>% dplyr::filter(Sample_Name %in% outliers), shape = 8, size = point_size_outlier, color = outlier_col)
    }
  } else {


    g <- ggplot(plotdata, aes_string(x = "IBS_mean", y = "IBS_var",
                                     color = if(colorbyQCmetric != "None") colorbyQCmetric else if(colorby != "None") colorby else NULL)) +
      geom_point(alpha=alpha, size = point_size) +
      theme +
      theme(
        text = element_text(size = text_size)
      ) +
      colorscheme

    if(length(selected_samples) > 0) {
      g <- g + geom_point(data = plotdata %>% dplyr::filter(Sample_Name %in% selected_samples), shape = 8, size = point_size_outlier, color = outlier_col)
    }
  }
  g
}


.XY_plot <- function(qcmetrics,
                    threshold = NULL,
                    colorby = "None",
                    colorbyQCmetric = "None",
                    outliers = c(),
                    show_outliers = FALSE,
                    selected_samples = c(),
                    background = NULL,
                    background_alpha = 0.3,
                    text_size = 13,
                    point_size = 2.5,
                    point_size_outlier = point_size * 1.5,
                    outlier_col = "red",
                    alpha = 0.8,
                    theme = theme_classic(),
                    colorscheme = NULL
) {
  plotdata <- qcmetrics@Sample_Metrics %>%
    dplyr::left_join(qcmetrics@samplesheet,by="Sample_Name") %>%
    dplyr::select(-predictedSex) %>%
    dplyr::left_join(getPMS(qcmetrics), by="Sample_Name") %>%
    dplyr::left_join(qcmetrics@sample_qc_outliers$sample_outliers[,!colnames(qcmetrics@sample_qc_outliers$sample_outliers) %in% setdiff(colnames(qcmetrics@samplesheet), "Sample_Name")], by="Sample_Name")
  names <- colnames(qcmetrics@sample_qc_outliers$sample_outliers)
  names <- names[names != "Sample_Name"]
  plotdata[names][is.na(plotdata[names])] <- FALSE

  threshold = if(is.null(threshold)) getThresholds(qcmetrics)$sample_qc_thresholds$XY_diff else threshold

  if((colorby != "None" || colorbyQCmetric != "None") && is.null(colorscheme)) {
    if(colorbyQCmetric != "None") {
      if(is.factor(plotdata[[colorbyQCmetric]]) || is.character(plotdata[[colorbyQCmetric]]) || is.logical(plotdata[[colorbyQCmetric]] )) {
        colorscheme <- scale_color_viridis_d()
      } else {
        colorscheme <- scale_color_viridis_c()
      }
    } else {
      if(is.factor(plotdata[[colorby]]) || is.character(plotdata[[colorby]]) || is.logical(plotdata[[colorby]] )) {
        colorscheme <- scale_color_viridis_d()
      } else {
        colorscheme <- scale_color_viridis_c()
      }
    }
  } else if(is.null(colorscheme)) {
    colorscheme <- scale_color_viridis_c()
  }

  if(!is.null(background)) {
    background <- background@Sample_Metrics

    # Check if any background data is also in present in foreground data
    overlap <- sum(background$Sample_Name %in% plotdata$Sample_Name)
    if(overlap > 0) {
      print(sprintf("%s samples in the foreground data are present in the background data. They will be removed from the background.", overlap))
      background <- background %>% dplyr::filter(!(Sample_Name %in% plotdata$Sample_Name))
    }
    if(nrow(background) == 0) {
      # If there is total overlap between foreground and background, the background option is turned off
      background <- NULL
    } else {
      # Create hexabin backgrounds
      background <- ggplot(background, aes(x = xMed, y = yMed)) +
        geom_hex(alpha = background_alpha, show.legend=FALSE) +
        scale_fill_gradient2(low ="white", mid = "gray95", high = "#0072B5")
    }
  }

  if(show_outliers) {

    if(!is.null(background)) {
      g <- background + geom_point(data = plotdata,
                                   aes_string(
                                     color = if(colorbyQCmetric != "None") colorbyQCmetric else if(colorby != "None") colorby else NULL,
                                     shape = "Sex"
                                   ),
                                   alpha = alpha,
                                   size = point_size)  +
      theme +
      theme(
          text = element_text(size = text_size)
        ) +
        colorscheme
    } else {
      g <- ggplot(plotdata, aes_string(x = "xMed", y = "yMed",
                                       color = if(colorbyQCmetric != "None") colorbyQCmetric else if(colorby != "None") colorby else NULL,
                                       shape = "Sex"
                                       )) +
        geom_point(alpha = alpha, size = point_size) +
        theme +
        theme(
          text = element_text(size = text_size)
        ) +
        colorscheme
    }

    if(length(selected_samples) > 0) {
      g <- g + geom_point(data = plotdata %>% dplyr::filter(Sample_Name %in% selected_samples), shape = 8, size = point_size_outlier, color = outlier_col)
    } else if(show_outliers && length(outliers) > 0) {
      g <- g + geom_point(data = plotdata %>% dplyr::filter(Sample_Name %in% outliers), shape = 8, size = point_size_outlier, color = outlier_col)
    }
  } else {

    if(!is.null(background)) {
      g <- background + geom_point(data = plotdata, aes_string(
        color = if(colorbyQCmetric != "None") colorbyQCmetric else if(colorby != "None") colorby else NULL,
        shape = "Sex"
        ),
        size = point_size,
        alpha = alpha
      ) +
        theme +
        theme(
          text = element_text(size = text_size)
        ) +
        colorscheme
    } else {
      g <- ggplot(plotdata, aes_string(x = "xMed", y = "yMed",
                                       color = if(colorbyQCmetric != "None") colorbyQCmetric else if(colorby != "None") colorby else NULL,
                                       shape = "Sex"
                                       )) +
        geom_point(alpha=alpha, size = point_size) +
        theme +
        theme(
          text = element_text(size = text_size)
        ) +
        colorscheme
    }

    if(length(selected_samples) > 0) {
      g <- g + geom_point(data = plotdata %>% dplyr::filter(Sample_Name %in% selected_samples), shape = 8, size = point_size_outlier, color = outlier_col)
    }
  }
  g
}


.XYdiff_plot <- function(qcmetrics,
                        threshold = NULL,
                        colorby = "None",
                        colorbyQCmetric = "None",
                        outliers = c(),
                        show_outliers = FALSE,
                        selected_samples = c(),
                        background = NULL,
                        text_size = 13,
                        point_size = 2.5,
                        point_size_outlier = point_size * 1.5,
                        alpha = 0.8,
                        outlier_col = "red",
                        theme = theme_classic(),
                        colorscheme = NULL
) {
  plotdata <- qcmetrics@Sample_Metrics %>%
    dplyr::left_join(qcmetrics@samplesheet,by="Sample_Name") %>%
    dplyr::select(-predictedSex) %>%
    dplyr::left_join(getPMS(qcmetrics), by="Sample_Name") %>%
    dplyr::left_join(qcmetrics@sample_qc_outliers$sample_outliers[,!colnames(qcmetrics@sample_qc_outliers$sample_outliers) %in% setdiff(colnames(qcmetrics@samplesheet), "Sample_Name")], by="Sample_Name")
  names <- colnames(qcmetrics@sample_qc_outliers$sample_outliers)
  names <- names[names != "Sample_Name"]
  plotdata[names][is.na(plotdata[names])] <- FALSE

  threshold = if(is.null(threshold)) getThresholds(qcmetrics)$sample_qc_thresholds$XY_diff else threshold

  if((colorby != "None" || colorbyQCmetric != "None") && is.null(colorscheme)) {
    if(colorbyQCmetric != "None") {
      if(is.factor(plotdata[[colorbyQCmetric]]) || is.character(plotdata[[colorbyQCmetric]]) || is.logical(plotdata[[colorbyQCmetric]] )) {
        colorscheme <- scale_color_viridis_d()
      } else {
        colorscheme <- scale_color_viridis_c()
      }
    } else {
      if(is.factor(plotdata[[colorby]]) || is.character(plotdata[[colorby]]) || is.logical(plotdata[[colorby]] )) {
        colorscheme <- scale_color_viridis_d()
      } else {
        colorscheme <- scale_color_viridis_c()
      }
    }
  } else if(is.null(colorscheme)) {
    colorscheme <- scale_color_viridis_c()
  }

  if(show_outliers) {

    g <- ggplot(plotdata, aes_string(x = "XYdiff", y = "Sample_Name",
                                     color = if(colorbyQCmetric != "None") colorbyQCmetric else if(colorby != "None") colorby else NULL,
                                     shape = "Sex"
                                     )) +
      geom_point(alpha = alpha, size = point_size) +
      geom_vline(xintercept = threshold, linetype = "dashed") +
      theme +
      theme(
        text = element_text(size = text_size),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank()
      ) +
      ylab("Samples") +
      colorscheme

    if(length(selected_samples) > 0) {
      g <- g + geom_point(data = plotdata %>% dplyr::filter(Sample_Name %in% selected_samples), shape = 8, size = point_size_outlier, color = outlier_col)
    } else if(show_outliers && length(outliers) > 0) {
      g <- g + geom_point(data = plotdata %>% dplyr::filter(Sample_Name %in% outliers), shape = 8, size = point_size_outlier, color = outlier_col)
    }
  } else {

    g <- ggplot(plotdata, aes_string(x = "XYdiff", y = "Sample_Name",
                                     color = if(colorbyQCmetric != "None") colorbyQCmetric else if(colorby != "None") colorby else NULL,
                                     shape = "Sex"
                                     )) +
      geom_point(alpha=alpha, size = point_size) +
      theme +
      theme(
        text = element_text(size = text_size),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank()
      ) +
      ylab("Samples") +
      colorscheme

    if(length(selected_samples) > 0) {
      g <- g + geom_point(data = plotdata %>% dplyr::filter(Sample_Name %in% selected_samples), shape = 8, size = point_size_outlier, color = outlier_col)
    }
  }
  g
}


## PCA plots (option: controlPCA, array-wide)
.PCA_plot <- function(qcmetrics,
                      type = c("arraywide", "control"),
                      threshold = NULL,
                      colorby = "None",
                      colorbyQCmetric = "None",
                      outliers = c(),
                      show_outliers = FALSE,
                      selected_samples = c(),
                      background = NULL,
                      text_size = 13,
                      point_size = 2.5,
                      point_size_outlier = point_size * 1.5,
                      outlier_col = "red",
                      alpha = 0.8,
                      theme = theme_classic(),
                      colorscheme = NULL
) {
  type <- match.arg(type)
  plotdata <- qcmetrics@Sample_Metrics %>%
    dplyr::left_join(qcmetrics@samplesheet,by="Sample_Name") %>%
    dplyr::select(-predictedSex) %>%
    dplyr::left_join(getPMS(qcmetrics), by="Sample_Name") %>%
    dplyr::left_join(qcmetrics@sample_qc_outliers$sample_outliers[,!colnames(qcmetrics@sample_qc_outliers$sample_outliers) %in% setdiff(colnames(qcmetrics@samplesheet), "Sample_Name")], by="Sample_Name")

  if(type == "control") {
    plotdata <- plotdata %>% dplyr::left_join(qcmetrics@controlPCs, by="Sample_Name")
  } else if (type == "arraywide") {
    plotdata <- plotdata %>% dplyr::left_join(qcmetrics@PCs, by="Sample_Name")
  }
  names <- colnames(qcmetrics@sample_qc_outliers$sample_outliers)
  names <- names[names != "Sample_Name"]
  plotdata[names][is.na(plotdata[names])] <- FALSE

  threshold <- if(type == "control" && is.null(threshold) && !is.null(getThresholds(qcmetrics)$sample_qc_thresholds$controlPCA_sd) && !is.na(getThresholds(qcmetrics)$sample_qc_thresholds$controlPCA_sd)) getThresholds(qcmetrics)$sample_qc_thresholds$controlPCA_sd else threshold
  threshold <- if(type == "arraywide" && is.null(threshold) && !is.null(getThresholds(qcmetrics)$sample_qc_thresholds$PCA_sd) && !is.na(getThresholds(qcmetrics)$sample_qc_thresholds$PCA_sd)) getThresholds(qcmetrics)$sample_qc_thresholds$PCA_sd else threshold

  if(!is.null(threshold) && type == "control") {
    controlPCs_melt <- qcmetrics@controlPCs %>% tidyr::gather(key = "PC_control", value = "value", dplyr::ends_with("_control"))
    dat <- controlPCs_melt %>% dplyr::mutate(PC_control = factor(PC_control, levels = paste0("PC", 1:10, "_control"))) %>%
      dplyr::group_by(PC_control) %>%
      dplyr::mutate(mean = mean(value), sd = sqrt(var(value)),
                    sd_min = mean - as.numeric(threshold) * sd, sd_plus = mean + as.numeric(threshold) * sd,
                    outlier = ifelse(value > sd_plus | value < sd_min, TRUE, FALSE))

    PCs_melt_sel <- dat[dat$PC_control %in% c("PC1_control", "PC2_control"),]

  } else if(!is.null(threshold) && type == "arraywide") {
    PCs_melt <- qcmetrics@PCs %>% tidyr::gather(key = "PC", value = "value", dplyr::starts_with("PC"))
    dat <- PCs_melt %>% dplyr::mutate(PC = factor(PC, levels = paste0("PC", 1:10))) %>%
      dplyr::group_by(PC) %>%
      dplyr::mutate(mean = mean(value), sd = sqrt(var(value)),
                    sd_min = mean - as.numeric(threshold) * sd, sd_plus = mean + as.numeric(threshold) * sd,
                    outlier = ifelse(value > sd_plus | value < sd_min, TRUE, FALSE))

    PCs_melt_sel <- dat[dat$PC %in% c("PC1", "PC2"),]
  }

  if((colorby != "None" || colorbyQCmetric != "None") && is.null(colorscheme)) {
    if(colorbyQCmetric != "None") {
      if(is.factor(plotdata[[colorbyQCmetric]]) || is.character(plotdata[[colorbyQCmetric]]) || is.logical(plotdata[[colorbyQCmetric]] )) {
        colorscheme <- scale_color_viridis_d()
      } else {
        colorscheme <- scale_color_viridis_c()
      }
    } else {
      if(is.factor(plotdata[[colorby]]) || is.character(plotdata[[colorby]]) || is.logical(plotdata[[colorby]] )) {
        colorscheme <- scale_color_viridis_d()
      } else {
        colorscheme <- scale_color_viridis_c()
      }
    }
  } else if(is.null(colorscheme)) {
    colorscheme <- scale_color_viridis_c()
  }

  if(!is.null(outliers) && !is.null(threshold)) {
    if(type == "control") {
      g <- ggplot(plotdata, aes_string(x = "PC1_control", y = "PC2_control",
                                       color = if(colorbyQCmetric != "None") colorbyQCmetric else if(colorby != "None") colorby else NULL)) +
        geom_point(alpha = alpha, size = point_size) +
        geom_vline(xintercept = dat %>% dplyr::filter(PC_control == "PC1_control") %$% sd_min[1], linetype = "dashed") +
        geom_vline(xintercept = dat %>% dplyr::filter(PC_control == "PC1_control") %$% sd_plus[1], linetype = "dashed") +
        geom_hline(yintercept = dat %>% dplyr::filter(PC_control == "PC2_control") %$% sd_min[1], linetype = "dashed") +
        geom_hline(yintercept = dat %>% dplyr::filter(PC_control == "PC2_control") %$% sd_plus[1], linetype = "dashed") +
        theme +
        theme(
          text = element_text(size = text_size)
        ) +
        colorscheme
    } else {
      g <- ggplot(plotdata, aes_string(x = "PC1", y = "PC2",
                                       color = if(colorbyQCmetric != "None") colorbyQCmetric else if(colorby != "None") colorby else NULL)) +
        geom_point(alpha = alpha, size = point_size) +
        geom_vline(xintercept = dat %>% dplyr::filter(PC == "PC1") %$% sd_min[1], linetype = "dashed") +
        geom_vline(xintercept = dat %>% dplyr::filter(PC == "PC1") %$% sd_plus[1], linetype = "dashed") +
        geom_hline(yintercept = dat %>% dplyr::filter(PC == "PC2") %$% sd_min[1], linetype = "dashed") +
        geom_hline(yintercept = dat %>% dplyr::filter(PC == "PC2") %$% sd_plus[1], linetype = "dashed") +
        theme +
        theme(
          text = element_text(size = text_size)
        ) +
        colorscheme
    }

  } else {

    if(type == "control") {
      g <- ggplot(plotdata, aes_string(x = "PC1_control", y = "PC2_control",
                                       color = if(colorbyQCmetric != "None") colorbyQCmetric else if(colorby != "None") colorby else NULL)) +
        geom_point(alpha=alpha, size = point_size) +
        theme +
        theme(
          text = element_text(size = text_size)
        ) +
        colorscheme
    } else {
      g <- ggplot(plotdata, aes_string(x = "PC1", y = "PC2",
                                       color = if(colorbyQCmetric != "None") colorbyQCmetric else if(colorby != "None") colorby else NULL)) +
        geom_point(alpha=alpha, size = point_size) +
        theme +
        theme(
          text = element_text(size = text_size)
        ) +
        colorscheme
    }
  }
  if(length(selected_samples) > 0) {
    g <- g + geom_point(data = plotdata %>% dplyr::filter(Sample_Name %in% selected_samples), shape = 8, size = point_size_outlier, color = outlier_col)
  } else if(show_outliers && length(outliers) > 0) {
    g <- g + geom_point(data = plotdata %>% dplyr::filter(Sample_Name %in% outliers), shape = 8, size = point_size_outlier, color = outlier_col)
  }
  g
}

## Probe QC
.beadNr_probes_plot <- function(qcmetrics,
                              threshold = NULL,
                              plot_cutoff = 0.01,
                              outliers = c(),
                              show_outliers = FALSE,
                              text_size = 13,
                              point_size = 2.5,
                              point_size_outlier = point_size * 1.5,
                              outlier_col = "red",
                              alpha = 0.8,
                              theme = theme_classic(),
                              colorscheme = NULL) {

  threshold <- getThresholds(qcmetrics)$probe_qc_thresholds$beadNr
  threshold2 <- getThresholds(qcmetrics)$beadnr_threshold

  if(show_outliers) {
    g <- qcmetrics@Probe_Metrics %>% dplyr::filter(BeadNr_Percentage >= plot_cutoff) %>%
      ggplot(aes(x = Probe, y = BeadNr_Percentage)) +
      geom_point(alpha=alpha, size = point_size) +
      theme +
      theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size = text_size),
        axis.text = element_text(size = text_size)) +
      geom_hline(yintercept = threshold, linetype = "dashed", color = "red") +
      ylab(sprintf("Proportion of samples with beadNr < %s", threshold2)) +
      xlab("Probes")
  } else {
    g <- qcmetrics@Probe_Metrics %>% dplyr::filter(BeadNr_Percentage >= plot_cutoff) %>%
      ggplot(aes(x = Probe, y = BeadNr_Percentage)) +
      geom_point(alpha=alpha, size = point_size) +
      theme +
      theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size = text_size),
        axis.text = element_text(size = text_size)) +
      ylab(sprintf("Proportion of samples with beadNr < %s", threshold2)) +
      xlab("Probes")
  }

  if(show_outliers && length(outliers) > 0) {
    g <- g + geom_point(data = qcmetrics@Probe_Metrics %>%
                          dplyr::filter(Probe %in% outliers), aes(x = Probe, y = BeadNr_Percentage), shape = 8, size = point_size_outlier,color = outlier_col)
  }
  g
}



.detP_probes_plot <- function(qcmetrics,
                              threshold = NULL,
                              plot_cutoff = 0.01,
                              outliers = c(),
                              show_outliers = FALSE,
                              text_size = 13,
                              point_size = 2.5,
                              point_size_outlier = point_size * 1.5,
                              outlier_col = "red",
                              alpha = 0.8,
                              theme = theme_classic(),
                              colorscheme = NULL) {

  threshold <- getThresholds(qcmetrics)$probe_qc_thresholds$detP
  threshold2 <- getThresholds(qcmetrics)$detp_threshold

  if(show_outliers) {
    g <- qcmetrics@Probe_Metrics %>% dplyr::filter(detP_Percentage >= plot_cutoff) %>%
      ggplot(aes(x = Probe, y = detP_Percentage)) +
      geom_point(alpha=alpha, size = point_size) +
      theme +
      theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size = text_size),
        axis.text = element_text(size = text_size)) +
      geom_hline(yintercept = threshold, linetype = "dashed", color = "red") +
      ylab(sprintf("Proportion of samples with detP > %s", threshold2)) +
      xlab("Probes")
  } else {
    g <- qcmetrics@Probe_Metrics %>% dplyr::filter(detP_Percentage >= plot_cutoff) %>%
      ggplot(aes(x = Probe, y = detP_Percentage)) +
      geom_point(alpha=alpha, size = point_size) +
      theme +
      theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size = text_size),
        axis.text = element_text(size = text_size)) +
      ylab(sprintf("Proportion of samples with detP > %s", threshold2)) +
      xlab("Probes")
  }

  if(show_outliers && length(outliers) > 0) {
    g <- g + geom_point(data = qcmetrics@Probe_Metrics %>%
                          dplyr::filter(Probe %in% outliers), aes(x = Probe, y = detP_Percentage), shape = 8, size = point_size_outlier, color = outlier_col)
  }
  g
}
