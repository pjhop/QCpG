shiny_PCA_plot <- function(
  scatterPCcolumns,
  melt_PCs,
  define_outliers,
  webGL,
  plotdata,
  colorbyQCmetric,
  scatterColorByColumn,
  selected_samples,
  outliers,
  PCA_outliers,
  show_outliers,
  colorscheme
  ) {
  if(define_outliers == 1) {

    g <- ggplot(plotdata, aes_string(x = scatterPCcolumns[1], y = scatterPCcolumns[2],
                                       color = if(colorbyQCmetric != "None") colorbyQCmetric else if(scatterColorByColumn != "None") scatterColorByColumn else NULL,
                                       key = "Sample_Name")) +
      geom_point(alpha=0.8) +
      theme_minimal() +
      colorscheme

    if(PCA_outliers == 1) {
      g <- g + geom_vline(xintercept = melt_PCs %>% dplyr::filter(PC == scatterPCcolumns[1]) %$% sd_min[1], linetype = "dashed") +
        geom_vline(xintercept = melt_PCs %>% dplyr::filter(PC == scatterPCcolumns[1]) %$% sd_plus[1], linetype = "dashed") +
        geom_hline(yintercept = melt_PCs %>% dplyr::filter(PC == scatterPCcolumns[2]) %$% sd_min[1], linetype = "dashed") +
        geom_hline(yintercept = melt_PCs %>% dplyr::filter(PC == scatterPCcolumns[2]) %$% sd_plus[1], linetype = "dashed")
    }
    if(length(selected_samples) > 0) {
      g <- g + geom_point(data = plotdata %>% dplyr::filter(Sample_Name %in% selected_samples),
                          shape = 23,
                          size = 2.2,
                          color = "darkturquoise",
                          fill = "darkturquoise")
      ggply <- ggplotly(g)
      ggply <- style(ggply,
                     hoverinfo = "none",
                     traces = length(ggply$x$data))

    } else if(show_outliers == 1 && length(outliers) > 0) {
      g <- g + geom_point(data = plotdata %>% dplyr::filter(Sample_Name %in% outliers),
                          shape = 8,
                          size = 2.2,
                          color = "red")
      ggply <- ggplotly(g)
      ggply <- style(ggply,
                     hoverinfo = "none",
                     traces = length(ggply$x$data))
    } else {
      ggply <- ggplotly(g)
    }
    if(webGL) toWebGL(style(ggply, hoveron = NULL)) else ggply

  } else {

    g <- ggplot(plotdata,
                aes_string(x = scatterPCcolumns[1],
                           y = scatterPCcolumns[2],
                           color = if(scatterColorByColumn == "None") NULL else scatterColorByColumn,
                           key = "Sample_Name")) +
      geom_point(alpha=0.8) +
      theme_minimal() +
      colorscheme

    if(length(selected_samples) > 0) {
      g <- g + geom_point(data = plotdata %>% dplyr::filter(Sample_Name %in% selected_samples),
                          shape = 23,
                          size = 2.2,
                          color = "darkturquoise",
                          fill = "darkturquoise")
      ggply <- ggplotly(g)
      ggply <- style(ggply,
                     hoverinfo = "none",
                     traces = length(ggply$x$data))
    } else {
      ggply <- ggplotly(g)
    }
    if(webGL) toWebGL(style(ggply, hoveron = NULL)) else ggply
  }
}

shiny_controlPCA_plot <- function(
                     scattercontrolPCcolumns,
                     melt_controlPCs,
                     define_outliers, webGL,
                     plotdata, colorbyQCmetric, scatterColorByColumn,
                     selected_samples,
                     outliers,
                     controlPCA_outliers,
                     show_outliers,
                     colorscheme
                     ) {
  if(define_outliers == 1) {
    g <- ggplot(plotdata, aes_string(x = scattercontrolPCcolumns[1],
                                     y = scattercontrolPCcolumns[2],
                                     color = if(colorbyQCmetric != "None") colorbyQCmetric else if(scatterColorByColumn != "None") scatterColorByColumn else NULL,
                                     key = "Sample_Name")) +
      geom_point(alpha=0.8) +
      theme_minimal() +
      colorscheme

    if(controlPCA_outliers == 1) {
      g <- g + geom_vline(xintercept = melt_controlPCs %>% dplyr::filter(PC_control == scattercontrolPCcolumns[1]) %$% sd_min[1], linetype = "dashed") +
        geom_vline(xintercept = melt_controlPCs %>% dplyr::filter(PC_control == scattercontrolPCcolumns[1]) %$% sd_plus[1], linetype = "dashed") +
        geom_hline(yintercept = melt_controlPCs %>% dplyr::filter(PC_control == scattercontrolPCcolumns[2]) %$% sd_min[1], linetype = "dashed") +
        geom_hline(yintercept = melt_controlPCs %>% dplyr::filter(PC_control == scattercontrolPCcolumns[2]) %$% sd_plus[1], linetype = "dashed")
    }

    if(length(selected_samples) > 0) {
      g <- g + geom_point(data = plotdata %>% dplyr::filter(Sample_Name %in% selected_samples),
                          shape = 23,
                          size = 2.2,
                          color = "darkturquoise",
                          fill = "darkturquoise")
      ggply <- ggplotly(g)
      ggply <- style(ggply,
                     hoverinfo = "none",
                     traces = length(ggply$x$data))

    } else if(show_outliers == 1 && length(outliers) > 0) {
      g <- g + geom_point(data = plotdata %>% dplyr::filter(Sample_Name %in% outliers),
                          shape = 8,
                          size = 2.2,
                          color = "red")
      ggply <- ggplotly(g)
      ggply <- style(ggply,
                     hoverinfo = "none",
                     traces = length(ggply$x$data))
    } else {
      ggply <- ggplotly(g)
    }
    if(webGL) toWebGL(style(ggply, hoveron = NULL)) else ggply

  } else {
    g <- ggplot(plotdata,
                aes_string(x = scattercontrolPCcolumns[1],
                           y = scattercontrolPCcolumns[2],
                           color = if(scatterColorByColumn == "None") NULL else scatterColorByColumn, key = "Sample_Name")) +
      geom_point(alpha=0.8) +
      theme_minimal() +
      colorscheme

    if(length(selected_samples) > 0) {
      g <- g + geom_point(data = plotdata %>% dplyr::filter(Sample_Name %in% selected_samples),
                          shape = 23,
                          size = 2.2,
                          color = "darkturquoise",
                          fill = "darkturquoise"
                          )
      ggply <- ggplotly(g)
      ggply <- style(ggply,
                     hoverinfo = "none",
                     traces = length(ggply$x$data))
    } else {
      ggply <- ggplotly(g)
    }
    if(webGL) toWebGL(style(ggply, hoveron = NULL)) else ggply
  }
}

shiny_XYplot <- function(define_outliers, background, show_background, webGL,
                   plotdata, colorbyQCmetric, scatterColorByColumn,
                   selected_samples, outliers,
                   show_outliers, bg, colorscheme) {
  data <- plotdata %>%
    dplyr::mutate(text = sprintf("Reported Sex: %s<br>Predicted Sex: %s", Sex, predictedSex))

  if(define_outliers == 1) {

    if(!is.null(background) && show_background && !webGL) {
      g <- bg + geom_point(data = data,
                              aes_string(key = "Sample_Name",
                                         text = "text",
                                         color = if(colorbyQCmetric != "None") colorbyQCmetric else if(scatterColorByColumn != "None") scatterColorByColumn else NULL),
                              alpha = 0.8) +
        theme_minimal() +
        colorscheme
    } else {
      g <- ggplot(data,
                  aes_string(x = "xMed",
                             y = "yMed",
                             key = "Sample_Name",
                             text = "text",
                             color = if(colorbyQCmetric != "None") colorbyQCmetric else if(scatterColorByColumn != "None") scatterColorByColumn else NULL)) +
        geom_point(alpha = 0.8) +
        theme_minimal() +
        colorscheme
    }

    if(length(selected_samples) > 0) {
      g <- g + geom_point(data = data %>% dplyr::filter(Sample_Name %in% selected_samples),
                          shape = 23,
                          size = 2.2,
                          color = "darkturquoise",
                          fill = "darkturquoise")
      ggply <- ggplotly(g)
      traces_bg <- which(purrr::map_chr(ggply$x$data, "mode") != "markers")
      ggply <- style(ggply,
                     hoverinfo = "none",
                     traces = if(!is.null(background) && show_background && !webGL) c(traces_bg,length(ggply$x$data)) else length(ggply$x$data))

    } else if(show_outliers == 1 && length(outliers) > 0) {
      g <- g + geom_point(data = data %>% dplyr::filter(Sample_Name %in% outliers),
                          shape = 8,
                          size = 2.2,
                          color = "red")
      ggply <- ggplotly(g)
      traces_bg <- which(purrr::map_chr(ggply$x$data, "mode") != "markers")
      ggply <- style(ggply,
                     hoverinfo = "none",
                     traces = if(!is.null(background) && show_background && !webGL) c(traces_bg,length(ggply$x$data)) else length(ggply$x$data))
    } else {
      ggply <- ggplotly(g)
      if(!is.null(background) && show_background && !webGL) {
        traces_bg <- which(purrr::map_chr(ggply$x$data, "mode") != "markers")
        ggply <- style(ggply, hoverinfo = "none", traces = traces_bg)
      }
    }
    if(webGL) toWebGL(style(ggply, hoveron = NULL)) else ggply

  } else {
    if(!is.null(background) && show_background && !webGL) {
      g <- bg + geom_point(data = data,
                              aes_string(key = "Sample_Name",
                                         text = "text",
                                         color = if(scatterColorByColumn == "None") NULL else scatterColorByColumn),
                              alpha = 0.8) +
        theme_minimal() +
        colorscheme
    } else {
      g <- ggplot(data, aes_string(x = "xMed",
                                   y = "yMed",
                                   key = "Sample_Name",
                                   text = "text",
                                   color = if(scatterColorByColumn == "None") NULL else scatterColorByColumn)) +
        geom_point(alpha=0.8) +
        theme_minimal() +
        colorscheme

    }

    if(length(selected_samples) > 0) {
      g <- g + geom_point(data = data %>% dplyr::filter(Sample_Name %in% selected_samples),
                          shape = 23,
                          size = 2.2,
                          color = "darkturquoise",
                          fill = "darkturquoise")
      ggply <- ggplotly(g)
      traces_bg <- which(purrr::map_chr(ggply$x$data, "mode") != "markers")
      ggply <- style(ggply,
                     hoverinfo = "none",
                     traces = if(!is.null(background) && show_background && !webGL) c(traces_bg,length(ggply$x$data)) else length(ggply$x$data))
    } else {
      ggply <- ggplotly(g)
      if(!is.null(background) && show_background && !webGL) {
        traces_bg <- which(purrr::map_chr(ggply$x$data, "mode") != "markers")
        ggply <- style(ggply, hoverinfo = "none", traces = traces_bg)
      }
    }
    if(webGL) toWebGL(style(ggply, hoveron = NULL)) else ggply

  }
}

shiny_XYdiff_plot <- function(XY_diff,define_outliers, webGL,
                   plotdata, colorbyQCmetric, scatterColorByColumn,
                   selected_samples, outliers,
                   show_outliers, colorscheme) {
  data <- plotdata %>%
    dplyr::mutate(text = sprintf("Reported Sex: %s<br>Predicted Sex: %s", Sex, predictedSex))

  if(define_outliers == 1) {

    g <- ggplot(data,
                aes_string(x = "XYdiff",
                           y = "Sample_Name",
                           text = "text",
                           key = "Sample_Name",
                           color = if(colorbyQCmetric != "None") colorbyQCmetric else if(scatterColorByColumn != "None") scatterColorByColumn else NULL
    )) +
      geom_jitter(alpha=0.8) +
      theme_minimal2() +
      theme(panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.text.y = element_blank()) +
      geom_vline(xintercept = XY_diff,
                 linetype = "dashed") +
      ylab("Samples") +
      colorscheme

    if(length(selected_samples) > 0) {
      g <- g + geom_point(data = data %>% dplyr::filter(Sample_Name %in% selected_samples),
                          shape = 23,
                          size = 2.2,
                          color = "darkturquoise",
                          fill = "darkturquoise")
      ggply <- ggplotly(g)
      ggply <- style(ggply,
                     hoverinfo = "none",
                     traces = length(ggply$x$data))

    } else if(show_outliers == 1 && length(outliers) > 0) {
      g <- g + geom_point(data = data %>% dplyr::filter(Sample_Name %in% outliers),
                          shape = 8,
                          size = 2.2,
                          color = "red")
      ggply <- ggplotly(g)
      ggply <- style(ggply,
                     hoverinfo = "none",
                     traces = length(ggply$x$data))
    } else {
      ggply <- ggplotly(g)
    }
    if(webGL) toWebGL(style(ggply, hoveron = NULL)) else ggply

  } else {
    g <- ggplot(data,
                aes_string(x = "XYdiff",
                           y = "Sample_Name",
                           text = "text",
                           key = "Sample_Name",
                           color = if(scatterColorByColumn == "None") NULL else scatterColorByColumn
    )) +
      geom_jitter(alpha=0.8) +
      theme_minimal2() +
      theme(panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.text.y = element_blank()) +
      ylab("Samples") +
      colorscheme

    if(length(selected_samples) > 0) {
      g <- g + geom_point(data = data %>% dplyr::filter(Sample_Name %in% selected_samples),
                          shape = 23,
                          size = 2.2,
                          color = "darkturquoise",
                          fill = "darkturquoise")
      ggply <- ggplotly(g)
      ggply <- style(ggply,
                     hoverinfo = "none",
                     traces = length(ggply$x$data))
    } else {
      ggply <- ggplotly(g)
    }
    if(webGL) toWebGL(style(ggply, hoveron = NULL)) else ggply

  }
}

shiny_plot_general <- function(type, define_outliers, background, show_background, webGL,
                    plotdata, colorbyQCmetric, scatterColorByColumn,
                    thresholds, selected_samples, outliers,
                    show_outliers, bg, colorscheme) {

if(define_outliers == 1) {

  if(!is.null(background) && show_background && !webGL) {
    if(type == "OP") {
      g <- bg + geom_point(data = plotdata,
                           aes_string(key = "Sample_Name",
                                      color = if(colorbyQCmetric != "None") colorbyQCmetric else if(scatterColorByColumn != "None") scatterColorByColumn else NULL),
                           alpha = 0.8) +
        geom_vline(xintercept = thresholds[1],
                   linetype = "dashed") +
        theme_minimal() +
        colorscheme
    } else if (type == "MU") {
      g <- bg + geom_point(data = plotdata,
                           aes_string(key = "Sample_Name",
                                      color = if(colorbyQCmetric != "None") colorbyQCmetric else if(scatterColorByColumn != "None") scatterColorByColumn else NULL,
                                      alpha = 0.8)) +
        xlim(0, NA) +
        ylim(0, NA) +
        geom_vline(xintercept = thresholds[1],
                   linetype = "dashed") +
        geom_hline(yintercept = thresholds[1],
                   linetype = "dashed") +
        theme_minimal() +
        colorscheme
    } else if (type == "RG") {
      g <- bg + geom_point(data = plotdata,
                           aes_string(key = "Sample_Name",
                                      color = if(colorbyQCmetric != "None") colorbyQCmetric else if(scatterColorByColumn != "None") scatterColorByColumn else NULL), alpha = 0.8) +
        geom_abline(slope = 1/thresholds[1],
                    linetype = "dashed") +
        geom_abline(slope = thresholds[2],
                    linetype = "dashed") +
        xlim(0, NA) +
        ylim(0, NA) +
        theme_minimal() +
        colorscheme
    } else if (type == "HC") {
      g <- bg + geom_point(data = plotdata,
                           aes_string(key = "Sample_Name",
                                      color = if(colorbyQCmetric != "None") colorbyQCmetric else if(scatterColorByColumn != "None") scatterColorByColumn else NULL),
                              alpha = 0.8) +
        geom_vline(xintercept = thresholds[1],
                   linetype = "dashed") +
        theme_minimal() +
        colorscheme
    }

  } else {
    if(type == "OP") {
      g <- ggplot(plotdata,
                  aes_string(x = "OP_x",
                             y = "OP_y",
                             key = "Sample_Name",
                             color = if(colorbyQCmetric != "None") colorbyQCmetric else if(scatterColorByColumn != "None") scatterColorByColumn else NULL)) +
        geom_point(alpha = 0.8) +
        geom_vline(xintercept = thresholds[1],
                   linetype = "dashed") +
        theme_minimal() +
        colorscheme
    } else if (type == "MU") {
      g <- ggplot(plotdata,
                  aes_string(x = "Unmethylated",
                             y = "Methylated",
                             key = "Sample_Name",
                             color = if(colorbyQCmetric != "None") colorbyQCmetric else if(scatterColorByColumn != "None") scatterColorByColumn else NULL)) +
        geom_point(alpha = 0.8) +
        xlim(0, NA) +
        ylim(0, NA) +
        geom_vline(xintercept = thresholds[1],
                   linetype = "dashed") +
        geom_hline(yintercept = thresholds[1],
                   linetype = "dashed") +
        theme_minimal() +
        colorscheme
    } else if (type == "RG") {
      g <- ggplot(plotdata,
                  aes_string(x = "Red",
                             y = "Green",
                             key = "Sample_Name",
                             color = if(colorbyQCmetric != "None") colorbyQCmetric else if(scatterColorByColumn != "None") scatterColorByColumn else NULL)) +
        geom_point(alpha = 0.8) +
        geom_abline(slope = 1/thresholds[1],
                    linetype = "dashed") +
        geom_abline(slope = thresholds[2],
                    linetype = "dashed") +
        xlim(0, NA) +
        ylim(0, NA) +
        theme_minimal() +
        colorscheme
    } else if (type == "HC") {
      g <- ggplot(plotdata,
                  aes_string(x = "HC_x",
                             y = "HC_y",
                             key = "Sample_Name",
                             color = if(colorbyQCmetric != "None") colorbyQCmetric else if(scatterColorByColumn != "None") scatterColorByColumn else NULL)) +
        geom_point(alpha = 0.8) +
        geom_vline(xintercept = thresholds[1],
                   linetype = "dashed") +
        theme_minimal() +
        colorscheme
    }
  }

  if(length(selected_samples) > 0) {
    g <- g + geom_point(data = plotdata %>% dplyr::filter(Sample_Name %in% selected_samples),
                        shape = 23,
                        size = 2.2,
                        color = "darkturquoise",
                        fill = "darkturquoise")
    ggply <- ggplotly(g)
    traces_bg <- which(purrr::map_chr(ggply$x$data, "mode") != "markers")
    ggply <- style(ggply,
                   hoverinfo = "none",
                   traces = if(!is.null(background) && show_background && !webGL) c(traces_bg,length(ggply$x$data)) else length(ggply$x$data))
  } else if(show_outliers == 1 && length(outliers) > 0) {
    g <- g + geom_point(data = plotdata %>% dplyr::filter(Sample_Name %in% outliers),
                        shape = 8,
                        size = 2.2,
                        color = "red")
    ggply <- ggplotly(g)
    traces_bg <- which(purrr::map_chr(ggply$x$data, "mode") != "markers")
    ggply <- style(ggply,
                   hoverinfo = "none",
                   traces = if(!is.null(background) && show_background && !webGL) c(traces_bg,length(ggply$x$data)) else length(ggply$x$data))
  } else {
    ggply <- ggplotly(g)
    if(!is.null(background) && show_background && !webGL) {
      traces_bg <- which(purrr::map_chr(ggply$x$data, "mode") != "markers")
      ggply <- style(ggply, hoverinfo = "none", traces = traces_bg)
    }
  }
  if(webGL) toWebGL(style(ggply, hoveron = NULL)) else ggply

} else {

  if(!is.null(background) && show_background && !webGL) {
    g <- bg + geom_point(data = plotdata,
                         aes_string(key = "Sample_Name",
                                    color = if(scatterColorByColumn == "None") NULL else scatterColorByColumn),
                            alpha = 0.8) +
      theme_minimal() +
      colorscheme
  } else {
    if(type == "OP") {
      g <- ggplot(plotdata,
                  aes_string(x = "OP_x",
                             y = "OP_y",
                             key = "Sample_Name",
                             color = if(scatterColorByColumn == "None") NULL else scatterColorByColumn)) +
        geom_point(alpha=0.8) +
        theme_minimal() +
        colorscheme
    } else if (type == "MU") {
      g <- ggplot(plotdata,
                  aes_string(x = "Unmethylated",
                             y = "Methylated",
                             key = "Sample_Name",
                             color = if(scatterColorByColumn == "None") NULL else scatterColorByColumn)) +
        geom_point(alpha=0.8) +
        xlim(0, NA) +
        ylim(0, NA) +
        theme_minimal() +
        colorscheme
    } else if (type == "RG") {
      g <- ggplot(plotdata,
                  aes_string(x = "Red",
                             y = "Green",
                             key = "Sample_Name",
                             color = if(scatterColorByColumn == "None") NULL else scatterColorByColumn)) +
        geom_point(alpha=0.8) +
        xlim(0, NA) +
        ylim(0, NA) +
        theme_minimal() +
        colorscheme
    } else if (type == "HC") {
      g <- ggplot(plotdata,
                  aes_string(x = "HC_x",
                             y = "HC_y",
                             key = "Sample_Name",
                             color = if(scatterColorByColumn == "None") NULL else scatterColorByColumn)) +
        geom_point(alpha=0.8) +
        theme_minimal() +
        colorscheme
    }
  }

  if(length(selected_samples) > 0) {
    g <- g + geom_point(data = plotdata %>% dplyr::filter(Sample_Name %in% selected_samples),
                        shape = 23,
                        size = 2.2,
                        color = "darkturquoise",
                        fill = "darkturquoise")
    ggply <- ggplotly(g)
    traces_bg <- which(purrr::map_chr(ggply$x$data, "mode") != "markers")
    ggply <- style(ggply,
                   hoverinfo = "none",
                   traces = if(!is.null(background) && show_background && !webGL) c(traces_bg,length(ggply$x$data)) else length(ggply$x$data))
  } else {
    ggply <- ggplotly(g)
    if(!is.null(background) && show_background && !webGL) {
      traces_bg <- which(purrr::map_chr(ggply$x$data, "mode") != "markers")
      ggply <- style(ggply, hoverinfo = "none", traces = traces_bg)
    }
  }
  if(webGL) toWebGL(style(ggply, hoveron = NULL)) else ggply

}}

shiny_bscon_plot <- function(bscon_threshold, define_outliers, webGL,
                       plotdata, colorbyQCmetric, scatterColorByColumn,
                       thresholds, selected_samples, outliers,
                       show_outliers, colorscheme) {
  if(define_outliers == 1) {

    g <- ggplot(plotdata, aes_string(x = "bscon",
                                     y = "Sample_Name",
                                     color = if(colorbyQCmetric != "None") colorbyQCmetric else if(scatterColorByColumn != "None") scatterColorByColumn else NULL,
                                      key = "Sample_Name")) +
      geom_point(alpha=0.8) +
      geom_vline(xintercept = bscon_threshold,
                 linetype = "dashed") +
      theme_minimal2() +
      theme(panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.text.y = element_blank()) +
      ylab("Samples") +
      xlab("Percentage") +
      colorscheme

    if(length(selected_samples) > 0) {
      g <- g + geom_point(data = plotdata %>% dplyr::filter(Sample_Name %in% selected_samples),
                          shape = 23,
                          size = 2.2,
                          color = "darkturquoise",
                          fill = "darkturquoise")
      ggply <- ggplotly(g)
      ggply <- style(ggply, hoverinfo = "none", traces = length(ggply$x$data))
    } else if(show_outliers == 1 && length(outliers) > 0) {
      g <- g + geom_point(data = plotdata %>% dplyr::filter(Sample_Name %in% outliers),
                          size = 2.2,
                          shape = 8,
                          colour = "red")
      ggply <- ggplotly(g)
      ggply <- style(ggply, hoverinfo = "none", traces = length(ggply$x$data))
    } else {
      ggply <- ggplotly(g)
    }
    if(webGL) toWebGL(style(ggply, hoveron = NULL)) else ggply

  } else {

    g <- ggplot(plotdata, aes_string(x = "bscon",
                                     y = "Sample_Name",
                                     color = if(scatterColorByColumn == "None") NULL else scatterColorByColumn,
                                     key = "Sample_Name")) +
      geom_point(alpha=0.8) +
      theme_minimal2() +
      theme(panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.text.y = element_blank()) +
      ylab("Samples") +
      xlab("Percentage") +
      colorscheme

    if(length(selected_samples) > 0) {
      g <- g + geom_point(data = plotdata %>% dplyr::filter(Sample_Name %in% selected_samples),
                          shape = 23,
                          size = 2.2,
                          color = "darkturquoise",
                          fill = "darkturquoise")
      ggply <- ggplotly(g)
      ggply <- style(ggply,
                     hoverinfo = "none",
                     traces = length(ggply$x$data))
    } else {
      ggply <- ggplotly(g)
    }
    if(webGL) toWebGL(style(ggply, hoveron = NULL)) else ggply
  }
}

# threshold = input$beadNr_threshold
# threshold2 = beadnr_threshold

# threshold = input$detP_threshold
# threshold2 = detp_threshold
shiny_detp_beadnr_plot <- function(
  type, define_outliers, webGL,
  plotdata, colorbyQCmetric, scatterColorByColumn,
  threshold, threshold2, selected_samples, outliers,
  show_outliers, colorscheme
) {
  if(define_outliers == 1) {
    if(type == "beadNr") {
      g <- ggplot(plotdata,
                  aes_string(x = "Sample_Name",
                             y = "BeadNr_Percentage",
                             color = if(colorbyQCmetric != "None") colorbyQCmetric else if(scatterColorByColumn != "None") scatterColorByColumn else NULL,
                             key = "Sample_Name"
      )) +
        geom_point(alpha=0.8) +
        theme_minimal2() +
        theme(panel.grid.major.x = element_blank(),
              panel.grid.minor.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.text.x = element_blank()) +
        geom_hline(yintercept = threshold,
                   linetype = "dashed") +
        ylab(sprintf("Proportion of probes with beadnr < %s", threshold2)) +
        xlab("Samples") +
        colorscheme
    } else if (type == "detP") {
      g <- ggplot(plotdata,
                  aes_string(x = "Sample_Name",
                             y = "detP_Percentage",
                             color = if(colorbyQCmetric != "None") colorbyQCmetric else if(scatterColorByColumn != "None") scatterColorByColumn else NULL,
                             key = "Sample_Name"
      )) +
        geom_point(alpha=0.8) +
        theme_minimal2() +
        theme(panel.grid.major.x = element_blank(),
              panel.grid.minor.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.text.x = element_blank()) +
        geom_hline(yintercept = threshold,
                   linetype = "dashed") +
        ylab(sprintf("Proportion of probes with detP > %s", threshold2)) +
        xlab("Samples") +
        colorscheme
    }

    if(length(selected_samples) > 0) {
      g <- g + geom_point(data = plotdata %>% dplyr::filter(Sample_Name %in% selected_samples),
                          shape = 23,
                          size = 2.2,
                          color = "darkturquoise",
                          fill = "darkturquoise")
      ggply <- ggplotly(g)
      ggply <- style(ggply, hoverinfo = "none", traces = length(ggply$x$data))
    } else if(show_outliers == 1 && length(outliers) > 0) {
      g <- g + geom_point(data = plotdata %>% dplyr::filter(Sample_Name %in% outliers),
                          size = 2.2,
                          shape = 8,
                          colour = "red")
      ggply <- ggplotly(g)
      ggply <- style(ggply, hoverinfo = "none", traces = length(ggply$x$data))
    } else {
      ggply <- ggplotly(g)
    }
    if(webGL) toWebGL(style(ggply, hoveron = NULL)) else ggply

  } else {

    if(type == "beadNr") {
      g <- ggplot(plotdata,
                  aes_string(x = "Sample_Name",
                             y = "BeadNr_Percentage",
                             color = if(scatterColorByColumn == "None") NULL else scatterColorByColumn,
                             key = "Sample_Name"
      )) +
        geom_point(alpha=0.8) +
        theme_minimal2() +
        theme(panel.grid.major.x = element_blank(),
              panel.grid.minor.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.text.x = element_blank()) +
        ylab(sprintf("Proportion of probes with beadnr < %s", threshold2)) +
        xlab("Samples") +
        colorscheme
    } else if (type == "detP") {
      g <- ggplot(plotdata,
                  aes_string(x = "Sample_Name",
                             y = "detP_Percentage",
                             color = if(scatterColorByColumn == "None") NULL else scatterColorByColumn,
                             key = "Sample_Name"
      )) +
        geom_point(alpha=0.8) +
        theme_minimal2() +
        theme(panel.grid.major.x = element_blank(),
              panel.grid.minor.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.text.x = element_blank()) +
        ylab(sprintf("Proportion of probes with detP > %s", threshold2)) +
        xlab("Samples") +
        colorscheme
    }

    if(length(selected_samples) > 0) {
      g <- g + geom_point(data = plotdata %>% dplyr::filter(Sample_Name %in% selected_samples),
                          shape = 23,
                          size = 2.2,
                          color = "darkturquoise",
                          fill = "darkturquoise")
      ggply <- ggplotly(g)
      ggply <- style(ggply, hoverinfo = "none", traces = length(ggply$x$data))
    } else {
      ggply <- ggplotly(g)
    }
    if(webGL) toWebGL(style(ggply, hoveron = NULL)) else ggply
  }
}

shiny_IBS_plot <- function(
  IBS_mean, IBS_var, define_outliers, webGL,
  plotdata, colorbyQCmetric, scatterColorByColumn,
  thresholds, selected_samples, outliers,
  show_outliers, colorscheme
)
  if(define_outliers == 1) {

    g <- ggplot(plotdata,
                aes_string(x = "IBS_mean",
                           y = "IBS_var",
                           color = if(colorbyQCmetric != "None") colorbyQCmetric else if(scatterColorByColumn != "None") scatterColorByColumn else NULL,
                           key = "Sample_Name")) +
      geom_point(alpha=0.8)  +
      geom_vline(xintercept = IBS_mean,
                 linetype = "dashed") +
      geom_hline(yintercept = IBS_var,
                 linetype = "dashed") +
      xlab("IBS Mean") +
      ylab("IBS Var") +
      theme_minimal() +
      colorscheme

    if(length(selected_samples) > 0) {
      g <- g + geom_point(data = plotdata %>% dplyr::filter(Sample_Name %in% selected_samples),
                          shape = 23,
                          size = 2.2,
                          color = "darkturquoise",
                          fill = "darkturquoise")
      ggply <- ggplotly(g)
      ggply <- style(ggply, hoverinfo = "none", traces = length(ggply$x$data))
    } else if(show_outliers == 1 && length(outliers) > 0) {
      g <- g + geom_point(data = plotdata %>% dplyr::filter(Sample_Name %in% outliers),
                          size = 2.2,
                          shape = 8,
                          colour = "red")
      ggply <- ggplotly(g)
      ggply <- style(ggply, hoverinfo = "none", traces = length(ggply$x$data))
    } else {
      ggply <- ggplotly(g)
    }
    if(webGL) toWebGL(style(ggply, hoveron = NULL)) else ggply

  } else {

    g <- ggplot(plotdata,
                aes_string(x = "IBS_mean",
                           y = "IBS_var",
                           color = if(scatterColorByColumn == "None") NULL else scatterColorByColumn,
                           key = "Sample_Name")) +
      geom_point(alpha=0.8)  +
      xlab("IBS Mean") +
      ylab("IBS Var") +
      theme_minimal() +
      colorscheme

    if(length(selected_samples) > 0) {
      g <- g + geom_point(data = plotdata %>% dplyr::filter(Sample_Name %in% selected_samples),
                          shape = 23,
                          size = 2.2,
                          color = "darkturquoise",
                          fill = "darkturquoise")
      ggply <- ggplotly(g)
      ggply <- style(ggply, hoverinfo = "none", traces = length(ggply$x$data))
    } else {
      ggply <- ggplotly(g)
    }
    if(webGL) toWebGL(style(ggply, hoveron = NULL)) else ggply
  }

shiny_relatedness_geno_plot <- function(
  qcmetrics,
  related_samples_geno,
  ibs_relatedness_mean, ibs_relatedness_var, define_outliers, webGL,
  plotdata, colorbyQCmetric, scatterColorByColumn,
  thresholds, selected_samples, outliers,
  show_outliers,
  show_pairs
) {
  qcmetrics@ibs$ibs_geno_relatedness <- qcmetrics@ibs$ibs_geno_relatedness %>% dplyr::mutate(key = sprintf("Sample_Name: %s\nSample_Name_geno: %s", Sample_Name, Sample_Name_geno))

  if(!qcmetrics@ibs$save_all_relationships) {
    g <- ggplot(qcmetrics@ibs$hexagons_relatedness_geno %>% dplyr::rename(IBS_mean=x, IBS_var=y),
                aes(x = IBS_mean,
                    y = IBS_var)) +
      geom_hex(mapping = aes(fill = count),
               stat = "identity",
               colour = NA,
               show.legend=FALSE)  +
      xlab("IBS Mean") +
      ylab("IBS Var") +
      theme_minimal() +
      xlim(NA, 2)  +
      ylim(0, NA)
    g <- g + geom_point(data = qcmetrics@ibs$ibs_geno_relatedness %>% dplyr::filter(IBS_mean > qcmetrics@relatedness_thresholds$IBS_mean & IBS_var < qcmetrics@relatedness_thresholds$IBS_var),
                        aes(text = key),
                        show.legend=FALSE)
  } else {
    g <- qcmetrics@ibs$ibs_geno_relatedness %>% dplyr::filter(IBS_mean < ibs_relatedness_mean | IBS_var > ibs_relatedness_var) %>%
      ggplot(aes(x = IBS_mean,
                 y = IBS_var)) +
      geom_hex(show.legend=FALSE) +
      xlab("IBS Mean") +
      ylab("IBS Var") +
      theme_minimal() +
      xlim(NA, 2)  +
      ylim(0, NA)
    g <- g + geom_point(data = qcmetrics@ibs$ibs_geno_relatedness %>% dplyr::filter(IBS_mean > ibs_relatedness_mean & IBS_var < ibs_relatedness_var),
                        aes(text = key),
                        show.legend=FALSE)
  }

  g <- g + geom_vline(xintercept = ibs_relatedness_mean,
                      linetype = "dashed") +
    geom_hline(yintercept = ibs_relatedness_var,
               linetype = "dashed")

  if(show_pairs) {
    pairs <- related_samples_geno
    if(nrow(pairs) > 0) {
      samples <- unique(c(pairs$Sample_Name, pairs$Sample_Name_geno))
      prs <-  t(apply(pairs[,c("Sample_Name", "Sample_Name_geno")], 1, sort))
      prs <- data.frame(prs[!duplicated(prs),])
      prs$pair_id <- paste0("pair", 1:nrow(prs))
      prs$pair <- paste0(prs$X1, prs$X2, sep = "_")
      tmp <- qcmetrics@ibs$ibs_geno_relatedness %>%
        dplyr::filter(Sample_Name %in% samples & Sample_Name_geno %in% samples)
      tmp.samples <- data.frame(t(apply(tmp[,c("Sample_Name", "Sample_Name_geno")], 1, sort)))
      tmp$pair <- paste0(tmp.samples$X1, tmp.samples$X2, sep = "_")
      tmp <- tmp %>% dplyr::left_join(prs[,c("pair", "pair_id")], by = "pair")
      tmp <- tmp %>% dplyr::filter(!is.na(pair_id))
      g <- g + geom_point(data = tmp, aes(text = key, x = IBS_mean, y = IBS_var, color = pair_id))
    }

  }

  if(show_outliers & define_outliers) {
    if(!qcmetrics@ibs$save_all_relationships) {
      tmp.data <- qcmetrics@ibs$ibs_geno_relatedness %>% dplyr::filter(IBS_mean > qcmetrics@relatedness_thresholds$IBS_mean & IBS_var < qcmetrics@relatedness_thresholds$IBS_var)
    } else {
      tmp.data <- qcmetrics@ibs$ibs_geno_relatedness %>% dplyr::filter(IBS_mean > ibs_relatedness_mean & IBS_var < ibs_relatedness_var)
    }
    tmp.data <- tmp.data %>% dplyr::filter(Sample_Name %in% outliers)
    g <- g + geom_point(data = tmp.data,
                        aes(text = key, x = IBS_mean, y = IBS_var),
                        shape = 8,
                        size = 2.2,
                        color = "red")
    ggply <- ggplotly(g)
    ggply <- style(ggply, hoverinfo = "none", traces = length(ggply$x$data))
  } else if(length(selected_samples) > 0) {
    if(!qcmetrics@ibs$save_all_relationships) {
      g <- g + geom_point(data = qcmetrics@ibs$ibs_geno_relatedness %>%
                            dplyr::filter(IBS_mean > qcmetrics@relatedness_thresholds$IBS_mean & IBS_var < qcmetrics@relatedness_thresholds$IBS_var) %>%
                            dplyr::filter(Sample_Name %in% selected_samples),
                          shape = 23,
                          size = 2.2,
                          color = "darkturquoise",
                          fill = "darkturquoise")
    } else {
      g <- g + geom_point(data = qcmetrics@ibs$ibs_geno_relatedness %>%
                            dplyr::filter(IBS_mean > ibs_relatedness_mean & IBS_var < ibs_relatedness_var) %>%
                            dplyr::filter(Sample_Name %in% selected_samples),
                          shape = 23,
                          size = 2.2,
                          color = "darkturquoise",
                          fill = "darkturquoise")
    }

    ggply <- ggplotly(g)
    ggply <- style(ggply, hoverinfo = "none", traces = length(ggply$x$data))
  } else {
    ggply <- ggplotly(g)
  }
  if(webGL) toWebGL(ggply) else ggply

}

shiny_relatedness_dnam_plot <- function(qcmetrics,
                                  related_samples_geno,
                                  ibs_relatedness_mean, ibs_relatedness_var, define_outliers, webGL,
                                  plotdata, colorbyQCmetric, scatterColorByColumn,
                                  thresholds, selected_samples, outliers,
                                  show_outliers
                                  ) {
  qcmetrics@ibs$ibs_dnam <- qcmetrics@ibs$ibs_dnam %>% dplyr::mutate(key = sprintf("Sample_Name.x: %s\nSample_Name.y: %s", Sample_Name.x, Sample_Name.y))

  if(!qcmetrics@ibs$save_all_relationships) {
    g <- ggplot(qcmetrics@ibs$hexagons_relatedness_dnam %>% dplyr::rename(IBS_mean=x, IBS_var=y),
                aes(x = IBS_mean, y = IBS_var)) +
      geom_hex(mapping = aes(fill = count),
               stat = "identity",
               show.legend=FALSE)  +
      xlab("IBS Mean") +
      ylab("IBS Var") +
      theme_minimal() +
      coord_cartesian(xlim = c(NA, 2.05), ylim = c(0, NA))

    g <- g + geom_point(data = qcmetrics@ibs$ibs_dnam %>% dplyr::filter(IBS_mean > qcmetrics@relatedness_thresholds$IBS_mean & IBS_var < qcmetrics@relatedness_thresholds$IBS_var),
                        aes(text = key),
                        show.legend=FALSE)
  } else {
    g <- qcmetrics@ibs$ibs_dnam %>% dplyr::filter(IBS_mean < ibs_relatedness_mean | IBS_var > ibs_relatedness_var) %>% ggplot(aes(x = IBS_mean, y = IBS_var)) +
      geom_hex(show.legend=FALSE) +
      xlab("IBS Mean") +
      ylab("IBS Var") +
      theme_minimal() +
      xlim(NA, 2.05) +
      ylim(0, NA)
    g <- g + geom_point(data = qcmetrics@ibs$ibs_dnam %>% dplyr::filter(IBS_mean > ibs_relatedness_mean & IBS_var < ibs_relatedness_var),
                        aes(text = key),
                        show.legend=FALSE)
  }

  g <- g + geom_vline(xintercept = ibs_relatedness_mean,
                      linetype = "dashed") +
    geom_hline(yintercept = ibs_relatedness_var,
               linetype = "dashed")


  if(define_outliers & show_outliers) {
    if(!qcmetrics@ibs$save_all_relationships) {
      tmp.data <- qcmetrics@ibs$ibs_dnam %>%
        dplyr::filter(IBS_mean > qcmetrics@relatedness_thresholds$IBS_mean & IBS_var < qcmetrics@relatedness_thresholds$IBS_var)
    } else {
      tmp.data <- qcmetrics@ibs$ibs_dnam %>%
        dplyr::filter(IBS_mean > ibs_relatedness_mean & IBS_var < ibs_relatedness_var)
    }

    #tmp.data <- tmp.data %>% filter(Sample_Name.x %in% outliers() | Sample_Name.y %in% outliers())
    # Update show whether one of the samples or both are outliers
    tmp.data <- tmp.data %>% dplyr::filter(Sample_Name.x %in% outliers | Sample_Name.y %in% outliers) %>%
      dplyr::mutate(Outlier = dplyr::case_when(
        (Sample_Name.x %in% outliers) & (Sample_Name.y %in% outliers) ~ "Both",
        (Sample_Name.x %in% outliers) | (Sample_Name.y %in% outliers) ~ "One sample",
        TRUE ~ "None"
      ))
    g <- g + geom_point(data = tmp.data,
                        aes(text = key, x = IBS_mean, y = IBS_var, color = Outlier))
    ggply <- ggplotly(g)
  } else if(length(selected_samples) > 0) {
    if(!qcmetrics@ibs$save_all_relationships) {
      g <- g + geom_point(data = qcmetrics@ibs$ibs_dnam %>% dplyr::filter(IBS_mean > qcmetrics@relatedness_thresholds$IBS_mean & IBS_var < qcmetrics@relatedness_thresholds$IBS_var) %>% dplyr::filter(Sample_Name.x %in% selected_samples),
                          shape = 23,
                          size = 2.2,
                          color = "darkturquoise",
                          fill = "darkturquoise")
    } else {
      g <- g + geom_point(data = qcmetrics@ibs$ibs_dnam %>% dplyr::filter(IBS_mean > ibs_relatedness_mean & IBS_var < ibs_relatedness_var) %>% dplyr::filter(Sample_Name.x %in% selected_samples),
                          shape = 23,
                          size = 2.2,
                          color = "darkturquoise",
                          fill = "darkturquoise")
    }
    ggply <- ggplotly(g)
    ggply <- style(ggply, hoverinfo = "none", traces = length(ggply$x$data))
  } else {
    ggply <- ggplotly(g)
  }

  if(webGL) toWebGL(ggply) else ggply
}

shiny_CTF_plot <- function(qcmetrics, define_outliers, webGL,
                     plotdata, colorbyQCmetric, scatterColorByColumn,
                     thresholds, selected_samples, outliers,
                     show_outliers, colorscheme) {
  if(define_outliers == 1) {
    data <- plotdata %>%
      tidyr::gather(key = "Cell_Type", value = "Proportion",
                    colnames(qcmetrics@cellcounts)[colnames(qcmetrics@cellcounts) != "Sample_Name"]) %>%
      dplyr::group_by(Cell_Type) %>%
      dplyr::mutate(quant = Proportion > quantile(Proportion, 0.95) | Proportion < quantile(Proportion, 0.05)) %>%
      dplyr::ungroup()

    g <- ggplot(data ,
                aes_string(x = "Cell_Type",
                           y = "Proportion",
                           color = if(colorbyQCmetric != "None") colorbyQCmetric else if(scatterColorByColumn != "None") scatterColorByColumn else NULL)) +
      geom_boxplot() +
      geom_point(data = data  %>% dplyr::filter(quant), aes_string(key = "Sample_Name"),
                 alpha = 0.8) +
      theme_minimal() +
      colorscheme

    if(length(selected_samples) > 0) {
      g <- g + geom_point(data = data  %>% dplyr::filter(Sample_Name %in% selected_samples),
                          shape = 23,
                          size = 2.2,
                          color = "darkturquoise",
                          fill = "darkturquoise")
      ggply <- ggplotly(g)
      ggply <- style(ggply,
                     hoverinfo = "none",
                     traces = length(ggply$x$data))
    } else if(show_outliers == 1) {
      g <- g + geom_point(data = data  %>% dplyr::filter(Sample_Name %in% outliers),
                          shape = 23,
                          size = 2.2,
                          color = "darkturquoise",
                          fill = "darkturquoise")
      ggply <- ggplotly(g)
      ggply <- style(ggply, hoverinfo = "none", traces = length(ggply$x$data))
    } else {
      ggply <- ggplotly(g)
    }
    ggply

  } else {
    data <- plotdata  %>%
      tidyr::gather(key = "Cell_Type",
                    value = "Proportion",
                    colnames(qcmetrics@cellcounts)[colnames(qcmetrics@cellcounts) != "Sample_Name"]) %>%
      dplyr::group_by(Cell_Type) %>%
      dplyr::mutate(quant = Proportion > quantile(Proportion, 0.95) | Proportion < quantile(Proportion, 0.05)) %>%
      dplyr::ungroup()

    g <- ggplot(data ,
                aes_string(x = "Cell_Type",
                           y = "Proportion",
                           color = if(scatterColorByColumn == "None") NULL else scatterColorByColumn)) +
      geom_boxplot() +
      geom_point(data = data  %>% dplyr::filter(quant),
                 aes_string(key = "Sample_Name"),
                 alpha = 0.8) +
      theme_minimal() +
      colorscheme


    if(length(selected_samples) > 0) {
      g <- g + geom_point(data = data  %>% dplyr::filter(Sample_Name %in% selected_samples),
                          shape = 23,
                          size = 2.2,
                          color = "darkturquoise",
                          fill = "darkturquoise")
      ggply <- ggplotly(g)
      ggply <- style(ggply,
                     hoverinfo = "none",
                     traces = length(ggply$x$data))
    } else {
      ggply <- ggplotly(g)
    }
    ggply

  }
}

shiny_age_plot <- function(qcmetrics, samplesheet, define_outliers, webGL,
                     plotdata, colorbyQCmetric, scatterColorByColumn,
                     thresholds, selected_samples, outliers,
                     show_outliers, select_age_predictor, colorscheme) {
  if("Age" %in% colnames(samplesheet) && !is.null(qcmetrics@pms)) {
    if(define_outliers == 1) {
      g <- ggplot(plotdata,
                  aes_string(x = "Age",
                             y = if(select_age_predictor == "Horvath") {
                               "Predicted_Age"} else if (select_age_predictor == "Zhang_ElasticNet") {
                                 "Predicted_Age_Zhang_en" } else {"Predicted_Age_Zhang_blup"} ,
                             key = "Sample_Name",
                             color = if(colorbyQCmetric != "None") colorbyQCmetric else if(scatterColorByColumn != "None") scatterColorByColumn else NULL)) +
        geom_point(alpha = 0.8) +
        theme_minimal() +
        colorscheme

      if(length(selected_samples) > 0) {
        g <- g + geom_point(data = plotdata %>% dplyr::filter(Sample_Name %in% selected_samples),
                            shape = 23,
                            size = 2.2,
                            color = "darkturquoise",
                            fill = "darkturquoise")
        ggply <- ggplotly(g)
        ggply <- style(ggply,
                       hoverinfo = "none",
                       traces = length(ggply$x$data))
      } else if(show_outliers == 1) {
        g <- g + geom_point(data = plotdata %>% dplyr::filter(Sample_Name %in% outliers),
                            shape = 8,
                            size = 2.2,
                            color = "red")
        ggply <- ggplotly(g)
        ggply <- style(ggply,
                       hoverinfo = "none",
                       traces = length(ggply$x$data))
      } else {
        ggply <- ggplotly(g)
      }
      ggply

    } else {
      g <- ggplot(plotdata,
                  aes_string(x = "Age",
                             y = if(select_age_predictor == "Horvath") {
                               "Predicted_Age"} else if (select_age_predictor == "Zhang_ElasticNet") {
                                 "Predicted_Age_Zhang_en" } else {"Predicted_Age_Zhang_blup"},
                             key = "Sample_Name",
                             color = if(scatterColorByColumn == "None") NULL else scatterColorByColumn)) +
        geom_point(alpha=0.8) +
        theme_minimal() +
        colorscheme

      if(length(selected_samples) > 0) {
        g <- g + geom_point(data = plotdata %>% dplyr::filter(Sample_Name %in% selected_samples),
                            shape = 23,
                            size = 2.2,
                            color = "darkturquoise",
                            fill = "darkturquoise")
        ggply <- ggplotly(g)
        ggply <- style(ggply,
                       hoverinfo = "none",
                       traces = length(ggply$x$data))
      } else {
        ggply <- ggplotly(g)
      }
      ggply

    }} else if(!is.null(qcmetrics@pms)) {
      g <- ggplot(data = plotdata,
                  aes_string(x = if(select_age_predictor == "Horvath") {
                    "Predicted_Age"} else if (select_age_predictor == "Zhang_ElasticNet") {
                      "Predicted_Age_Zhang_en" } else {"Predicted_Age_Zhang_blup"})) +
        geom_histogram() +
        theme_minimal() +
        colorscheme
      ggplotly(g)
    }
}

shiny_smoking_plot <- function(qcmetrics, samplesheet, define_outliers, webGL,
                         plotdata, colorbyQCmetric, scatterColorByColumn,
                         thresholds, selected_samples, outliers,
                         show_outliers, colorscheme) {
  if("Smoking_Status" %in% colnames(samplesheet)) {
    if(define_outliers == 1) {
      data <- plotdata %>%
        dplyr::group_by(Smoking_Status) %>%
        dplyr::mutate(quant = Smoking_PMS > quantile(Smoking_PMS, 0.95) | Smoking_PMS < quantile(Smoking_PMS, 0.05)) %>%
        dplyr::ungroup()

      if(scatterColorByColumn == "None" & colorbyQCmetric == "None") {
        g <- ggplot(data, aes(x = Smoking_Status, y = Smoking_PMS)) +
          geom_boxplot() +
          geom_point(data = data %>% dplyr::filter(quant),
                     aes(key = Sample_Name),
                     alpha = 0.8) +
          theme_minimal() +
          colorscheme

        if(length(selected_samples) > 0) {
          g <- g + geom_point(data = data %>% dplyr::filter(Sample_Name %in% selected_samples),
                              shape = 23,
                              size = 2.2,
                              color = "darkturquoise",
                              fill = "darkturquoise")
          ggply <- ggplotly(g)
          ggply <- style(ggply, hoverinfo = "none", traces = length(ggply$x$data))
        } else if(show_outliers == 1) {
          g <- g + geom_point(data = data %>% dplyr::filter(Sample_Name %in% outliers),
                              shape = 8,
                              size = 2.2,
                              color = "red")
          ggply <- ggplotly(g)
          ggply <- style(ggply, hoverinfo = "none", traces = length(ggply$x$data))
        } else {
          ggply <- ggplotly(g)
        }
        ggply

      } else {

        g <- ggplot(data, aes_string(x = "Smoking_Status", y = "Smoking_PMS",
                                     color = if(colorbyQCmetric != "None") colorbyQCmetric else scatterColorByColumn)) +
          geom_boxplot() +
          geom_point(data = data %>% dplyr::filter(quant),
                     aes(key = "Sample_Name"),
                     alpha = 0.8) +
          theme_minimal() +
          colorscheme

        if(length(selected_samples) > 0) {
          g <- g + geom_point(data = data %>% dplyr::filter(Sample_Name %in% selected_samples),
                              shape = 23,
                              size = 2.2,
                              color = "darkturquoise",
                              fill = "darkturquoise")
          ggply <- ggplotly(g)
          ggply <- style(ggply,
                         hoverinfo = "none",
                         traces = length(ggply$x$data))
        } else if(show_outliers == 1) {
          g <- g + geom_point(data = data %>% dplyr::filter(Sample_Name %in% outliers),
                              shape = 8,
                              size = 2.2,
                              color = "red")
          ggply <- ggplotly(g)
          ggply <- style(ggply,
                         hoverinfo = "none",
                         traces = length(ggply$x$data))
        } else {
          ggply <- ggplotly(g)
        }
        ggply

      }
    } else {
      data <-  plotdata %>%
        dplyr::group_by(Smoking_Status) %>%
        dplyr::mutate(quant = Smoking_PMS > quantile(Smoking_PMS, 0.95) | Smoking_PMS < quantile(Smoking_PMS, 0.05)) %>%
        dplyr::ungroup()

      if(scatterColorByColumn == "None") {
        g <- ggplot(data,
                    aes(x = Smoking_Status, y = Smoking_PMS)) +
          geom_boxplot() +
          geom_point(data = data %>% dplyr::filter(quant),
                     aes(key = Sample_Name),
                     alpha = 0.8) +
          theme_minimal() +
          colorscheme

        if(length(selected_samples) > 0) {
          g <- g + geom_point(data = data %>% dplyr::filter(Sample_Name %in% selected_samples),
                              shape = 23,
                              size = 2.2,
                              color = "darkturquoise",
                              fill = "darkturquoise")
          ggply <- ggplotly(g)
          ggply <- style(ggply,
                         hoverinfo = "none",
                         traces = length(ggply$x$data))
        } else {
          ggply <- ggplotly(g)
        }

        ggply
      } else {

        g <- ggplot(data,
                    aes_string(x = "Smoking_Status",
                               y = "Smoking_PMS",
                               color = scatterColorByColumn)) +
          geom_boxplot() +
          geom_point(data = data %>% dplyr::filter(quant),
                     aes_string(key = "Sample_Name"),
                     alpha = 0.8) +
          theme_minimal() +
          colorscheme

        if(length(selected_samples) > 0) {
          g <- g + geom_point(data = data %>% dplyr::filter(Sample_Name %in% selected_samples),
                              shape = 23,
                              size = 2.2,
                              color = "darkturquoise",
                              fill = "darkturquoise")
          ggply <- ggplotly(g)
          ggply <- style(ggply,
                         hoverinfo = "none",
                         traces = length(ggply$x$data))
        } else {
          ggply <- ggplotly(g)
        }
        ggply
      }
    }} else {
      g <- ggplot(data = plotdata,
                  aes(x = Smoking_PMS)) +
        geom_histogram() +
        theme_minimal() +
        colorscheme
      ggplotly(g)
    }
}

shiny_plot_PMS <- function(type, qcmetrics, samplesheet, define_outliers, webGL,
                         plotdata, colorbyQCmetric, scatterColorByColumn,
                         thresholds, selected_samples, outliers,
                         show_outliers, colorscheme) {
  if(type %in% colnames(samplesheet)) {
    if(define_outliers == 1) {

      if(type == "BMI") {
        g <- ggplot(plotdata,
                    aes_string(x = "BMI_PMS",
                               y = "BMI",
                               key = "Sample_Name",
                               color = if(colorbyQCmetric != "None") colorbyQCmetric else if(scatterColorByColumn != "None") scatterColorByColumn else NULL)) +
          geom_point(alpha = 0.8) +
          theme_minimal() +
          colorscheme
      } else if (type == "Alcohol") {
        g <- ggplot(plotdata,
                    aes_string(x = "Alcohol_PMS", y = "Alcohol", key = "Sample_Name",
                    color = if(colorbyQCmetric != "None") colorbyQCmetric else if(scatterColorByColumn != "None") scatterColorByColumn else NULL)) +
          geom_point(alpha = 0.8) +
          theme_minimal() +
          colorscheme
      }

      if(length(selected_samples) > 0) {
        g <- g + geom_point(data = plotdata %>% dplyr::filter(Sample_Name %in% selected_samples),
                            shape = 23,
                            size = 2.2,
                            color = "darkturquoise",
                            fill = "darkturquoise")
        ggply <- ggplotly(g)
        ggply <- style(ggply, hoverinfo = "none", traces = length(ggply$x$data))
      } else if(show_outliers == 1 && length(outliers) > 0) {
        g <- g + geom_point(data = plotdata %>% dplyr::filter(Sample_Name %in% outliers),
                            shape = 8,
                            size = 2.2,
                            color = "red")
        ggply <- ggplotly(g)
        ggply <- style(ggply,
                       hoverinfo = "none",
                       traces = length(ggply$x$data))
      } else {
        ggply <- ggplotly(g)
      }
      if(webGL) toWebGL(style(ggply, hoveron = NULL)) else ggply

    } else {
      if(type == "BMI") {
        g <- ggplot(plotdata,
                    aes_string(x = "BMI_PMS",
                               y = "BMI",
                               key = "Sample_Name",
                               color = if(scatterColorByColumn == "None") NULL else scatterColorByColumn)) +
          geom_point(alpha=0.8) +
          theme_minimal() +
          colorscheme
      } else if (type == "Alcohol") {
        g <- ggplot(plotdata,
                    aes_string(x = "Alcohol_PMS",
                               y = "Alcohol",
                               key = "Sample_Name",
                               color = if(scatterColorByColumn == "None") NULL else scatterColorByColumn)) +
          geom_point(alpha=0.8) +
          theme_minimal() +
          colorscheme
      }

      if(length(selected_samples) > 0) {
        g <- g + geom_point(data = plotdata %>% dplyr::filter(Sample_Name %in% selected_samples),
                            shape = 23,
                            size = 2.2,
                            color = "darkturquoise",
                            fill = "darkturquoise")
        ggply <- ggplotly(g)
        ggply <- style(ggply,
                       hoverinfo = "none",
                       traces = length(ggply$x$data))
      } else {
        ggply <- ggplotly(g)
      }
      if(webGL) toWebGL(style(ggply, hoveron = NULL)) else ggply
    }
  } else {
    if (type == "BMI") {
      g <- ggplot(data = plotdata,
                  aes(x = BMI_PMS)) +
        geom_histogram() +
        theme_minimal()
    } else if (type == "Alcohol") {
      g <- ggplot(data = plotdata,
                  aes(x = Alcohol_PMS)) +
        geom_histogram() +
        theme_minimal()
    }

    ggplotly(g)
  }
}

#
# MU_plot <- function(define_outliers, background, show_background, webGL,
#                     plotdata, colorbyQCmetric, scatterColorByColumn,
#                     MU_threshold, selected_samples, outliers,
#                     show_outliers, bg_mu
# ) {
#   if(define_outliers == 1) {
#
#     if(!is.null(background) && show_background && !webGL) {
#       g <- bg_mu + geom_point(data = plotdata, aes_string(key =
#                                                             "Sample_Name",
#                                                           color =
#                                                             if(colorbyQCmetric != "None") colorbyQCmetric else
#                                                               if(scatterColorByColumn != "None") scatterColorByColumn else NULL,
#                                                           alpha = 0.8)) +
#         geom_vline(xintercept = MU_threshold, linetype = "dashed") +
#         geom_hline(yintercept = MU_threshold, linetype = "dashed") +
#         theme_minimal()
#     } else {
#       g <- ggplot(plotdata, aes_string(x = "Unmethylated", y =
#                                          "Methylated", key = "Sample_Name",
#                                        color = if(colorbyQCmetric !=
#                                                   "None") colorbyQCmetric else if(scatterColorByColumn != "None")
#                                                     scatterColorByColumn else NULL)) +
#         geom_point(alpha = 0.8) +
#         geom_vline(xintercept = MU_threshold, linetype = "dashed") +
#         geom_hline(yintercept = MU_threshold, linetype = "dashed") +
#         theme_minimal()
#     }
#
#     if(length(selected_samples) > 0) {
#       g <- g + geom_point(data = plotdata %>% dplyr::filter(Sample_Name
#                                                             %in% selected_samples), shape = 8, size = 2.2, color = "red")
#       ggply <- ggplotly(g)
#       traces_bg_MU <- which(purrr::map_chr(ggply$x$data, "mode") !=
#                               "markers")
#       ggply <- style(ggply, hoverinfo = "none", traces =
#                        if(!is.null(background) && show_background && !webGL)
#                          c(traces_bg_MU,length(ggply$x$data)) else length(ggply$x$data))
#     } else if(show_outliers == 1 && length(outliers > 0)) {
#       g <- g + geom_point(data = plotdata %>% dplyr::filter(Sample_Name
#                                                             %in% outliers), shape = 8, size = 2.2, color = "red")
#       ggply <- ggplotly(g)
#       traces_bg_MU <- which(purrr::map_chr(ggply$x$data, "mode") !=
#                               "markers")
#       ggply <- style(ggply, hoverinfo = "none", traces =
#                        if(!is.null(background) && show_background && !webGL)
#                          c(traces_bg_MU,length(ggply$x$data)) else length(ggply$x$data))
#     } else {
#       ggply <- ggplotly(g)
#       if(!is.null(background) && show_background %% !webGL) {
#         traces_bg_MU <- which(purrr::map_chr(ggply$x$data, "mode") !=
#                                 "markers")
#         ggply <- style(ggply, hoverinfo = "none", traces = traces_bg_MU)
#       }
#     }
#     if(webGL) toWebGL(style(ggply, hoveron = NULL)) else ggply
#
#   } else {
#
#     if(!is.null(background) && show_background && !webGL) {
#       g <- bg_mu + geom_point(data = plotdata, aes_string(key =
#                                                             "Sample_Name",
#                                                           color =
#                                                             if(scatterColorByColumn == "None") NULL else scatterColorByColumn, alpha
#                                                           = 0.8)) +
#         theme_minimal()
#     } else {
#       g <- ggplot(plotdata, aes_string(x = "Unmethylated", y =
#                                          "Methylated", key = "Sample_Name",
#                                        color = if(scatterColorByColumn
#                                                   == "None") NULL else scatterColorByColumn)) +
#         geom_point(alpha=0.8) +  theme_minimal()
#     }
#
#     if(length(selected_samples) > 0) {
#       g <- g + geom_point(data = plotdata %>% dplyr::filter(Sample_Name
#                                                             %in% selected_samples), shape = 8, size = 2.2, color = "red")
#       ggply <- ggplotly(g)
#       traces_bg_MU <- which(purrr::map_chr(ggply$x$data, "mode") !=
#                               "markers")
#       ggply <- style(ggply, hoverinfo = "none", traces =
#                        if(!is.null(background) && show_background && !webGL)
#                          c(traces_bg_MU,length(ggply$x$data)) else length(ggply$x$data))
#     } else {
#       ggply <- ggplotly(g)
#       if(!is.null(background) && show_background && !webGL) {
#         traces_bg_MU <- which(purrr::map_chr(ggply$x$data, "mode") !=
#                                 "markers")
#         ggply <- style(ggply, hoverinfo = "none", traces = traces_bg_MU)
#       }
#     }
#     if(webGL) toWebGL(style(ggply, hoveron = NULL)) else ggply
#
#   } }
#
