## Check with dfOrNull representations
setClassUnion("dfOrNULL", c("data.frame", "NULL"))
setClassUnion("listOrNULL", c("list", "NULL"))

#' S4 class to represent QCmetrics
#'
#' Returned by the \code{\link{getQCmetrics}} function. The \code{\link{QCapp}} method interactively visualizes the QCmetrics object.
#' @slot Sample_Metrics Contains all sample QC metrics
#' @slot Probe_Metrics Contains all probe QC metrics
#' @slot samplesheet User-provided samplesheet
#' @slot ibs List containing IBS (identity-by-state) info.
#' @slot PCs Array-wide principal components
#' @slot controlPCs Control probe principal components
#' @slot heatmap_pvals Associations between array-wide PCs and technical/biological variables
#' @slot heatmap_pvals_controlPCA Associations between control-probe PCs and technical/biological variables
#' @slot cellcounts data.frame containing predicted white blood cell proportions
#' @slot controlmatrix matrix containing control probe info
#' @slot pms data.frame containg polymethylation scores for age, smoking, alcohol and BMI
#' @slot probe_qc_method Which probe QC method is used (before or after sample QC)
#' @slot sample_qc_outliers Sample outliers and thresholds
#' @slot probe_qc_outliers Probe outliers and thresholds
#' @slot detp_threshold Detection p-value threshold used in the `getQCmetrics()` function
#' @slot beadnr_threshold Beadnr threshold used in the `getQCmetrics()` function
#' @slot nprobes_PCA number of probes used for array-wide PCA, either "all" or a number.
#' @slot relatedness_thresholds thresholds used for defining related samples
#' @slot selected_samples Samples selected within the \code{\link{QCapp}}.
#' @slot metadata metadata
#' @export
setClass("QCmetrics",
         representation(Sample_Metrics = "data.frame",
                        Probe_Metrics = "listOrNULL",
                        samplesheet = "data.frame",
                        ibs = "list",
                        PCs = "dfOrNULL",
                        controlPCs = "data.frame",
                        heatmap_pvals = "dfOrNULL",
                        heatmap_pvals_controlPCA = "dfOrNULL",
                        cellcounts = "data.frame",
                        controlmatrix = "matrix",
                        pms = "data.frame",
                        probe_qc_method = "character",
                        sample_qc_outliers = "list",
                        probe_qc_outliers = "list",
                        detp_threshold = "numeric",
                        beadnr_threshold = "numeric",
                        nprobes_PCA = "character",
                        relatedness_thresholds = "list",
                        selected_samples = "character",
                        metadata = "list"
         )
)


#' @describeIn QCmetrics show summary info
#' @import methods
#' @param object a \linkS4class{QCmetrics} object
#' @export
setMethod("show", "QCmetrics",
          function(object) cat("class: QCmetrics",
                               sprintf("array: %s", object@metadata$array),
                               sprintf("dim: %s %s", BiocGenerics::nrow(object), BiocGenerics::ncol(object)),
                               sep = "\n"))

#' @describeIn QCmetrics number of samples
#' @param x An object of class \linkS4class{QCmetrics}
#' @export
setMethod("ncol", "QCmetrics",
          function(x) nrow(x@samplesheet))

#' @describeIn QCmetrics number of probes
#' @param x An object of class \linkS4class{QCmetrics}
#' @export
setMethod("nrow", "QCmetrics",
          function(x) nrow(x@Probe_Metrics))

#' @describeIn QCmetrics number of probes and number of samples
#' @param x An object of class \linkS4class{QCmetrics}
#' @export
setMethod("dim", "QCmetrics",
          function(x) c(BiocGenerics::nrow(x), BiocGenerics::ncol(x))
          )


#' @describeIn QCmetrics returns sample IDs
#' @param x An object of class \linkS4class{QCmetrics}
#' @export
setMethod("colnames", "QCmetrics",
          function(x) x@Sample_Metrics$Sample_Name)

#' @describeIn QCmetrics returns probe IDs
#' @param x An object of class \linkS4class{QCmetrics}
#' @export
setMethod("rownames", "QCmetrics",
          function(x) x@Probe_Metrics$Probe)

#' Return poly-methylation scores (PMS) for a \linkS4class{QCmetrics} object
#'
#' returns predicted phenotypes for \linkS4class{QCmetrics} objects
#' (three age predictors, sex, smoking score, alcohol score, BMI, and cell type fractions)
#' @param qcmetrics a \linkS4class{QCmetrics} object
#' @return a data.frame containing the poly-methylation scores
#' (three age predictors, sex, smoking score, alcohol score, BMI, and cell type fractions)
#' @export
setGeneric("getPMS",
           function(qcmetrics)
             standardGeneric("getPMS")
  )

#' @rdname getPMS
setMethod("getPMS", "QCmetrics",
          function(qcmetrics) {
            predictions <- qcmetrics@pms %>%
              dplyr::left_join(qcmetrics@Sample_Metrics[,c("Sample_Name", "predictedSex")], by = "Sample_Name") %>%
              dplyr::left_join(qcmetrics@cellcounts, by = "Sample_Name")
            predictions
          }
        )

#' Return array-wide PCs for \linkS4class{QCmetrics} objects
#'
#' Return array-wide PCs for \linkS4class{QCmetrics} objects
#' @param qcmetrics a \linkS4class{QCmetrics} object
#' @return a data.frame containing the array-wide PCs
#' @export
setGeneric("getPCs",
           function(qcmetrics)
             standardGeneric("getPCs")
)

#' @rdname getPCs
setMethod("getPCs", "QCmetrics",
          function(qcmetrics) {
            if(is.null(qcmetrics@PCs)) {
              message("Array-wide PCs are not available for this QCmetrics object, set `calculate_PCs = TRUE` in `getQCmetrics()` to return array-wide PCs.")
              return(NULL)
            } else {
              qcmetrics@PCs
            }
          }
)

#' Return control-probe PCs for \linkS4class{QCmetrics} objects
#'
#' Return control-probe PCs for \linkS4class{QCmetrics} objects
#' @param qcmetrics a \linkS4class{QCmetrics} object
#' @return a data.frame containing the array-wide PCs
#' @export
setGeneric("getControlPCs",
           function(qcmetrics)
             standardGeneric("getControlPCs")
)

#' @rdname getPCs
setMethod("getControlPCs", "QCmetrics",
          function(qcmetrics) {
              qcmetrics@controlPCs
          }
)

#' Return samplesheet in a \linkS4class{QCmetrics} object
#'
#' Return samplesheet in a \linkS4class{QCmetrics} object
#' @return a data.frame containing the sample info
#' @param qcmetrics a \linkS4class{QCmetrics} object
#' @export
setGeneric("getSamplesheet",
           function(qcmetrics)
             standardGeneric("getSamplesheet")
)

#' @rdname getSamplesheet
#' @param qcmetrics a \linkS4class{QCmetrics} object
setMethod("getSamplesheet", "QCmetrics",
          function(qcmetrics) {
            qcmetrics@samplesheet
          }
)

#' Combine two or more \linkS4class{QCmetrics} objects into one object.
#'
#' Combine two or more \linkS4class{QCmetrics} objects into one object.
#' This function will first check whether the metric-files were:
#' (1) measured using the same array
#' (2) generated using the same package version
#' (3) Used the same probe_qc_method
#' (4) The same detP and beadNr thresholds were used.
#' By default an error will be raised when any of the above checks fail.
#' @param x \linkS4class{QCmetrics} object
#' @param y \linkS4class{QCmetrics} object
#' @param ... \linkS4class{QCmetrics} objects
#' @param force force combining when the metric-files do not have the same package version,
#' probe_qc_method,detection P-value/beadnr threshold or are measured using different arrays.
#' @param include_probe_metrics Should probe metrics be included after combining? Defaults to TRUE.
#' @param npcs Number of control probe PCs to retain.
#' @param name Name of the experiment.
#' @param PCs Optionally, a dataframe containing array-wide PCs to add to the combined \linkS4class{QCmetrics}object.
#' @param heatmap_variables test for an association between PCs and these variables.
#' @param heatmap_variables_fct if `heatmap_variables` are specified, which of them should be treated as factors.
#' @return a \linkS4class{QCmetrics} object
#' @export
#'
setGeneric("combineQCmetrics",
           function(x, y, ...)
             standardGeneric("combineQCmetrics")
  )

#' @rdname combineQCmetrics
setMethod("combineQCmetrics", signature = c("QCmetrics", "QCmetrics"),
          function(x, y, ..., force = FALSE, include_probe_metrics = TRUE, npcs = 30, name = NULL, PCs = NULL, heatmap_variables = NULL, heatmap_variables_fct = NULL) {
            combine_metrics(x, y, ..., force = force, include_probe_metrics = include_probe_metrics, npcs = npcs, name = name,
                            PCs = PCs, heatmap_variables = heatmap_variables, heatmap_variables_fct = heatmap_variables_fct)
          }
  )

#' @rdname combineQCmetrics
setMethod("combineQCmetrics", signature = c("list", "missing"),
          function(x, force = FALSE, include_probe_metrics = TRUE, npcs = 30, name = NULL, PCs = NULL, heatmap_variables = NULL, heatmap_variables_fct = NULL) {
            combine_metrics(QCmetricslist = x, force = force, include_probe_metrics = include_probe_metrics, npcs = npcs, name = name,
                            PCs = PCs, heatmap_variables = heatmap_variables, heatmap_variables_fct = heatmap_variables_fct)
          }
)



#' Add PCs to a \linkS4class{QCmetrics} object
#'
#' Add PCs to a \linkS4class{QCmetrics} object
#' @param qcmetrics \linkS4class{QCmetrics} object
#' @param PCs data.frame containing PCs and a `Sample_Name` column
#' @param heatmap_variables test for an association between PCs and these variables.
#' @param heatmap_variables_fct if `heatmap_variables` are specified, which of them should be treated as factors.
#' @return a  \linkS4class{QCmetrics} object
#' @export
#'
setGeneric("addPCs",
           function(qcmetrics, PCs, heatmap_variables = NULL, heatmap_variables_fct = NULL)
             standardGeneric("addPCs")
  )

#' @rdname addPCs
setMethod("addPCs","QCmetrics",
           function(qcmetrics, PCs, heatmap_variables = NULL, heatmap_variables_fct = NULL) {
             stopifnot("Sample_Name" %in% colnames(PCs))
             # Check whether qcmetrics file and PCs contain the same samples
             test <- intersect(qcmetrics@Sample_Metrics$Sample_Name, PCs$Sample_Name)
             test <- length(test) == nrow(PCs) && length(test) == nrow(qcmetrics@Sample_Metrics)
             if(!test) {
               stop("QCmetrics object and PCs do not contain the same samples!")
             }
             # Check whether all columns (except Sample_Name) are named PC*
             col <- colnames(PCs)[colnames(PCs) != "Sample_Name"]
             test <- all(col %in% paste0("PC", 1:length(col)))
             if(!test) {
               stop("PC columns should be named PC1, PC2, etc.")
             }

             # Add
             PCs <- PCs %>% dplyr::select(Sample_Name, dplyr::everything())
             qcmetrics@PCs <- PCs

             # Calculate heatmap
             message("Recalculating associations for heatmap..")
             heatmap_pvals <- .get_heatmap_pvals(samplesheet = qcmetrics@samplesheet,
                                                 PCs = qcmetrics@PCs,
                                                 cellcounts = qcmetrics@cellcounts,
                                                 pms = qcmetrics@pms,
                                                 variables = heatmap_variables,
                                                 factor = heatmap_variables_fct,
                                                 logfile = NULL)
             qcmetrics@heatmap_pvals <- heatmap_pvals
             # Return
             qcmetrics
           }

)

#' Retrieve metadata for a \linkS4class{QCmetrics} object

#' @title retrieve metadata for a \linkS4class{QCmetrics} object
#' @param x \linkS4class{QCmetrics} object
#' @importFrom S4Vectors metadata
#' @export
#'
setMethod("metadata", "QCmetrics",
         function(x) {
           cat(sprintf("Experiment Name: %s\nExperimenter: %s\nArray: %s\nPackage Version: %s\nCreated at: %s",
                       if(!is.null(x@metadata$experiment_name)) x@metadata$experiment_name else "-",
                       if(!is.null(x@metadata$experimenter)) x@metadata$experimenter else "-",
                       x@metadata$array,
                       as.character(x@metadata$package_version),
                       x@metadata$analysis_date
                       ))
         }
  )

#' Get sample outliers from a \linkS4class{QCmetrics} object
#'
#' Get sample outliers from a \linkS4class{QCmetrics} object
#' @param qcmetrics A \linkS4class{QCmetrics} object
#' @export
#'
setGeneric("getSampleOutliers",
           function(qcmetrics)
             standardGeneric("getSampleOutliers")
)

#' @rdname getSampleOutliers
setMethod("getSampleOutliers", "QCmetrics",
          function(qcmetrics) {
            qcmetrics@sample_qc_outliers$sample_outliers %>% dplyr::filter(Total)
          }
)

#' Get related pair(s) of samples from a \linkS4class{QCmetrics} object
#'
#' Get related pair(s) of samples from a \linkS4class{QCmetrics} object
#' @param qcmetrics \linkS4class{QCmetrics} object
#' @export
#'
setGeneric("getRelatedPairs",
           function(qcmetrics)
             standardGeneric("getRelatedPairs")
)

#' @rdname getRelatedPairs
setMethod("getRelatedPairs", "QCmetrics",
          function(qcmetrics) {
            if(!is.null(qcmetrics@ibs$ibs_dnam)) {
              qcmetrics@ibs$ibs_dnam %>%
                dplyr::filter(IBS_mean > qcmetrics@relatedness_thresholds$IBS_mean,
                              IBS_var < qcmetrics@relatedness_thresholds$IBS_var) %>%
                dplyr::select(Sample_Name.x, Sample_Name.y, IBS_mean, IBS_var)
            } else {
              message("No relatedness info contained in this object!")
            }

          }
)

#' Get probe outliers from a \linkS4class{QCmetrics} object
#'
#' Get probe outliers from a \linkS4class{QCmetrics} object
#' @param qcmetrics \linkS4class{QCmetrics} object
#' @export
#'
setGeneric("getProbeOutliers",
           function(qcmetrics)
             standardGeneric("getProbeOutliers")
)

#' @rdname getProbeOutliers
setMethod("getProbeOutliers", "QCmetrics",
          function(qcmetrics) {
            qcmetrics@probe_qc_outliers$probe_outliers
          }
)


#' Retrieve thresholds from a \linkS4class{QCmetrics} object
#'
#' Retrieve thresholds from a \linkS4class{QCmetrics} object
#' @param qcmetrics \linkS4class{QCmetrics} object
#' @export
#'
setGeneric("getThresholds",
           function(qcmetrics)
             standardGeneric("getThresholds")
)

#' @rdname getThresholds
setMethod("getThresholds", "QCmetrics",
          function(qcmetrics) {
            list(sample_qc_thresholds = qcmetrics@sample_qc_outliers$thresholds,
                 probe_qc_thresholds = qcmetrics@probe_qc_outliers$thresholds,
                 relatedness_thresholds = qcmetrics@relatedness_thresholds,
                 detp_threshold = qcmetrics@detp_threshold,
                 beadnr_threshold = qcmetrics@beadnr_threshold
                 )
          }
)

#' Extract all relatedness info from a QCmetrics object
#'
#' Extract all relatedness info from a QCmetrics object
#' @param qcmetrics QCmetrics object
#' @export
#'

setGeneric("relatedness",
           function(qcmetrics)
             standardGeneric("relatedness")
)

#' @rdname relatedness
setMethod("relatedness", "QCmetrics",
          function(qcmetrics) {
            qcmetrics@ibs$ibs_dnam
          }
)

#' Launch QC app
#'
#' Visualize a \linkS4class{QCmetrics} object generated using the \code{\link{getQCmetrics}} function.
#' @param qcmetrics A \linkS4class{QCmetrics} object.
#' @param samplesheet A samplesheet (data.frame), optional argument.
#' By default, the samplesheet stored within the \linkS4class{QCmetrics} object will be used.
#' @param show_all_relationships Show all relationships in relatedness plots?
#' For large data this can slow down the app considerably given that the number of relations is N^2.
#' If `FALSE`, individual-level data is shown only for relations passing the relatedness threshold,
#' all other relations are summarized into hexbins.
#' @param background Object of class \linkS4class{QCmetrics} or class \linkS4class{QCmetricsBG}.
#' These metrics will be shown as background in the relevant plots.
#' @param variables Optional, which variables in samplesheet to use.
#' @import shiny
#' @import ggplot2
#' @importFrom plotly ggplotly style renderPlotly plotlyOutput event_data toWebGL
#' @examples
#' library(minfiData)
#' ## Load example data from minfiData package (6 samples)
#' library(minfiData)
#' library(dplyr)
#' baseDir <- system.file("extdata", package="minfiData")
#' samplesheet <- read.metharray.sheet(baseDir)
#'
#' # Make sure that the samplesheet contains a 'Sample_Name'
#' # column with the sample IDs and a 'Sex' column!
#' samplesheet <- samplesheet %>%
#'                    dplyr::mutate(Sample_Name =  paste(samplesheet$Slide,
#'                                  samplesheet$Array, sep = "_")) %>%
#'                    dplyr::rename(Sex = sex)
#'
#' ## Get QC metrics
#' metrics <- getQCmetrics(samplesheet,
#'   name = "example minfiData",
#'   calculate_PCs = TRUE,
#'   n_probes_PCA = 20000,
#'   heatmap_variables = c("Array", "Slide"),
#'   heatmap_variables_fct = c("Array", "Slide"))
#'
#' \dontrun{
#' QCapp(metrics)
#' }
#'
#' # By default the samplesheet used in the getQCmetrics() function will be used.
#' # However, users can provide an (updated) samplesheet.
#' # For example: rename the 'age' column to 'Age' so that it will be recognized and used
#' # to create a scatterplot between Age and Predicted Age
#' samplesheet <- samplesheet %>% dplyr::rename(Age = age)
#' \dontrun{
#' QCapp(metrics, samplesheet = samplesheet)
#' }
#' @export
#'
setGeneric("QCapp",
           function(qcmetrics, background = NULL, samplesheet = NULL, show_all_relationships = TRUE, variables = "all")
             standardGeneric("QCapp")
)

#' @rdname QCapp
setMethod("QCapp", "QCmetrics",
          function(qcmetrics, background = NULL, samplesheet = NULL, show_all_relationships = TRUE, variables = "all") {
            run_app(qcmetrics = qcmetrics, background = background, samplesheet = samplesheet, show_all_relationships = show_all_relationships, variables = variables)
          }
  )

#' Update a QCmetrics object with a cleaned/normalized beta-matrix
#'
#' Update a QCmetrics object with a cleaned/normalized beta-matrix
#' @param qcmetrics QCmetrics object
#' @param beta beta-matrix or a filepath pointing to a txt-file where the first column contains the probe identifiers, and the remaining columns are samples.
#' @param beta_format if a filepath is provided to the 'beta' parameter (instead of an R matrix), this parameter describes the format ('rds' or 'txt'). This is inferred from the file extension if this parameter is empty.
#' @param chunk_size chunk size, applicable if 'read_chunks'=TRUE and a filepath to a txt-file is provided in the 'beta' parameter.
#' @param read_chunks TRUE/FALSE, read the data in chunks? (applicable if 'beta' points toward a txt-file.)
#' @param delimiter Defaults to `\\t`, only applicable for txt-files.
#' @param samplesheet samplesheet
#' @param calculate_PCs should PCA be performed? (TRUE/FALSE). Default = TRUE
#' @param n_probes_PCA number of most variable probes to maintain for PCA, defaults to 'NULL' (all probes are used).
#' @param n_estimate_PCA number of samples used to estimate the most variable probes (defaults to 'chunk_size').
#' @param nrPCs nr of PCs to retain. Default = 10
#' @param epidish_method Which method to using in the EpiDish cellcount prediction (RPC, CBS, or CP)
#' @param heatmap_variables If heatmap = y, this parameter can be used to add extra variables (in a character vector) that should be tested for an association with PCs.
#' @param heatmap_variables_fct If extra variables are supplied for heatmap (heatmap_variables option), which of them should be treated as a factor?
#' @param verbose Should the function be verbose? (TRUE/FALSE)
#' @export
#'

setGeneric("updateQCmetrics",
           function(qcmetrics, beta, samplesheet = NULL, beta_format = NULL,
                    chunk_size = 500, read_chunks = TRUE,
                    delimiter = "\t",
                    calculate_PCs = FALSE,
                    n_probes_PCA = NULL,
                    n_estimate_PCA = chunk_size,
                    nrPCs = 10,
                    epidish_method = c("RPC", "CBS", "CP"), heatmap_variables = NULL,
                    heatmap_variables_fct = NULL, verbose = TRUE)
             standardGeneric("updateQCmetrics")
)

setMethod("updateQCmetrics", "QCmetrics",
          function(qcmetrics, beta, samplesheet = NULL, beta_format = NULL,
                   chunk_size = 500, read_chunks = TRUE,
                   delimiter = "\t", calculate_PCs = FALSE,
                   n_probes_PCA = NULL,
                   n_estimate_PCA = chunk_size,
                   nrPCs = 10, epidish_method = c("RPC", "CBS", "CP"),
                   heatmap_variables = NULL, heatmap_variables_fct = NULL, verbose = TRUE) {
                     update_function(qcmetrics = qcmetrics, beta = beta, samplesheet = samplesheet,
                                     beta_format = beta_format, chunk_size = chunk_size,
                                     read_chunks = read_chunks, delimiter = delimiter,
                                     calculate_PCs = calculate_PCs, nrPCs = nrPCs,
                                     n_probes_PCA = n_probes_PCA, n_estimate_PCA = n_estimate_PCA,
                                     epidish_method = epidish_method, heatmap_variables = heatmap_variables,
                                     heatmap_variables_fct = heatmap_variables_fct, verbose = verbose)
          }
)

#' Save QC plots contained in a \linkS4class{QCmetrics} object in high quality.
#'
#' Save QC plots contained in a \linkS4class{QCmetrics} object in high quality.
#' @param qcmetrics \linkS4class{QCmetrics} object
#' @param outdir path to the directory where the plots should be save (will be created if it does not exist)
#' @param name filename
#' @param background A \linkS4class{QCmetrics} or \linkS4class{QCmetricsBG} object that will be used as background.
#' Defaults to `NULL`.
#' @param background_alpha transparency, defaults to 0.3.
#' @param width Plot size in units ("in", "cm", or "mm"). Defaults to 7.
#' @param height Plot size in units ("in", "cm", or "mm"). Defaults to 5.
#' @param unit "in", "cm", or "mm". Defaults to "in"
#' @param filetype Device to use ("png", "eps", "ps", "tex" (pictex), "pdf", "jpeg", "tiff", "png", "bmp", "svg" or "wmf). Defaults to "png".
#' @param dpi Plot resolution. Defaults to 300.
#' @param colorby Column to color the plot by. Defaults to "None".
#' @param colorbyQCmetric Color by outlier type (MU_outlier etc.)
#' @param show_outliers Show outliers? (TRUE/FALSE). Based on outliers as defined in the \linkS4class{QCmetrics} object.
#' Defaults to `TRUE`.
#' @param show_which_outliers "all" or "metric".
#' The former will highlight all outliers (based on all metrics),
#' the latter will only the outliers based on the metric plotted.
#' @param selected_samples Highlight user-provided vector of Sample Names.
#' @param probes_plot_cutoff Limit the number of probes plotted in the probe QC plots.
#' Defaults to 0.01 (only showing probes with percentage > 0.01).
#' @param text_size Text size, defaults to 13.
#' @param point_size Point size, defaults to 2.5.
#' @param point_size_outlier Point size for highlighted outliers, defaults to
#' @param outlier_col Color used for highlighting outliers, defaults to "red".
#' @param alpha transparency, defaults to 0.8.
#' @param theme ggplot theme to use, defaults to "theme_classic()"
#' @param colorscheme which colorscheme to use, defaults to ggplot2::scale_color_viridis_c()/ggplot2::scale_color_viridis_d()
#' @export
#'
#' @examples
#' ## Load example data from minfiData package (6 samples)
#' library(minfiData)
#' library(dplyr)
#' baseDir <- system.file("extdata", package="minfiData")
#' samplesheet <- read.metharray.sheet(baseDir)

#' # Make sure that the samplesheet contains a 'Sample_Name'
#' # column with the sample IDs and a 'Sex' column!
#' samplesheet <- samplesheet %>%
#'                    dplyr::mutate(Sample_Name =  paste(samplesheet$Slide,
#'                                  samplesheet$Array, sep = "_")) %>%
#'                    dplyr::rename(Sex = sex)
#'
#' ## Get QC metrics
#' metrics <- getQCmetrics(samplesheet,
#'   name = "example minfiData",
#'   calculate_PCs = TRUE,
#'   n_probes_PCA = 20000,
#'   heatmap_variables = c("Array", "Slide"),
#'   heatmap_variables_fct = c("Array", "Slide"))
#'
#' ## Save plots
#' \dontrun{
#' savePlots(metrics,
#'           outdir = "./",
#'           name = "example_minfidata")
#'}
#'
#' ## color by 'Slide' variable and show metric-specific outliers:
#' \dontrun{
#' savePlots(metrics,
#'           outdir = "./",
#'           colorby = "Slide",
#'           show_which_outliers = "metric",
#'           name = "example_minfidata")
#'}
#'
#' ## Change dimensions and use pdf-format
#' \dontrun{
#' savePlots(metrics,
#'           outdir = "./",
#'           colorby = "Slide",
#'           show_which_outliers = "metric",
#'           name = "example_minfidata",
#'           width = 10, height = 8, filetype="pdf"
#'           )
#'}


setGeneric("savePlots",
           function(qcmetrics,
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
                    colorscheme = NULL)
             standardGeneric("savePlots")
)


setMethod("savePlots", "QCmetrics",
          function(qcmetrics,
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
                   colorscheme = NULL) {
            save_plots(qcmetrics = qcmetrics,
                       outdir = outdir,
                       name = name,
                       background = background,
                       background_alpha = background_alpha,
                       width = width,
                       height = width,
                       unit = unit,
                       filetype = filetype,
                       dpi = dpi,
                       colorby = colorby,
                       colorbyQCmetric = colorbyQCmetric,
                       show_outliers = show_outliers,
                       show_which_outliers = show_which_outliers,
                       selected_samples = selected_samples,
                       probes_plot_cutoff = probes_plot_cutoff,
                       text_size = text_size,
                       point_size = point_size,
                       point_size_outlier = point_size_outlier,
                       outlier_col = outlier_col,
                       alpha = alpha,
                       theme = theme,
                       colorscheme = colorscheme)
          }
)

#' getOutlierOverview
#'
#' Get an overview of the number of samples failed on each QC metric
#' @param qcmetrics QCmetrics object
#' @export
setGeneric("getOutlierOverview",
           function(qcmetrics)
           standardGeneric("getOutlierOverview")
)

setMethod("getOutlierOverview", "QCmetrics",
          function(qcmetrics) {
            outliers <- getSampleOutliers(qcmetrics)
            outlier_vars <- colnames(outliers)[!colnames(outliers) %in% c("Sample_Name", "Total")]
            outlier_table_temp <- outliers[,outlier_vars]

            add_unique_outliers <- function(var) {
              outliers_vars_current <- outlier_table_temp[,colnames(outlier_table_temp) == var]
              outliers_vars_other <- outlier_table_temp[,colnames(outlier_table_temp) != var]
              other_outliers <- apply(outliers_vars_other, 1, function(x) any(x))
              unique <- outliers_vars_current & !other_outliers
              unique
            }
            uniques <- lapply(outlier_vars, FUN = add_unique_outliers)
            uniques <- do.call(cbind, uniques)
            colnames(uniques) <- paste0(outlier_vars, "_unique")
            outliers <- cbind(outliers, uniques)

            getSampleOutliers(qcmetrics) %>%
              dplyr::select(-c("Sample_Name", "Total")) %>%
              tidyr::gather(key = "Metric", value = "Outliers", dplyr::everything()) %>%

              dplyr::group_by(Metric) %>%
              dplyr::summarize(Nr_Outliers = sum(Outliers), .groups = "drop") %>%
              dplyr::mutate(Metric = stringr::str_replace(Metric, "_outlier$", ""))

            outliers_total <- getSampleOutliers(qcmetrics) %>%
              tidyr::gather(key = "Metric", value = "Outliers", -c("Sample_Name", "Total")) %>%
              dplyr::group_by(Metric) %>%
              dplyr::summarize(Nr_Outliers = sum(Outliers), .groups = "drop") %>%
              dplyr::mutate(Metric = stringr::str_replace(Metric, "_outlier$", ""))

            outliers_unique <- outliers %>%
              tidyr::gather(key = "Metric", value = "Outliers", dplyr::ends_with("unique")) %>%
              dplyr::group_by(Metric) %>%
              dplyr::summarize(Nr_Unique_Outliers = sum(Outliers), .groups = "drop") %>%
              dplyr::mutate(Metric = stringr::str_replace(Metric, "_outlier_unique$|_unique$", ""))

            outliers_total <- outliers_total %>%
              dplyr::left_join(outliers_unique[,c("Metric", "Nr_Unique_Outliers")], by = "Metric")

            ## Add a total number
            outliers_total <- dplyr::bind_rows(
              outliers_total,
              tibble::tibble(Metric = c("Total"), Nr_Outliers = nrow(getSampleOutliers(qcmetrics)), Nr_Unique_Outliers = NA)
            )
            outliers_total
          }
)


#' Convert QCmetrics class to QCmetricsBG class, dropping al information not needed for background representation.
#'
#' Convert QCmetrics class to QCmetricsBG class, dropping al information not needed for background representation.
#' @param qcmetrics QCmetrics object
#' @export

setGeneric("makeBackground",
           function(qcmetrics)
           standardGeneric("makeBackground")
)


#' @rdname makeBackground
#'
setMethod("makeBackground", "QCmetrics",
          function(qcmetrics) {
            new("QCmetricsBG",
                Sample_Metrics = qcmetrics@Sample_Metrics,
                metadata = qcmetrics@metadata)
          }
)


#' S4 class to represent QCmetrics background
#'
#' Returned by the \code{\link{makeBackground}} method.
#' @slot Sample_Metrics Contains all sample QC metrics
#' @slot metadata metadata
#' @export
setClass("QCmetricsBG",
         representation(Sample_Metrics = "data.frame",
                        metadata = "list"
         ))

#' @describeIn QCmetricsBG show summary info
#' @import methods
#' @param object a QCmetricsBG object
#' @export
setMethod("show", "QCmetricsBG",
          function(object) cat("class: QCmetricsBG",
                               sprintf("array: %s", object@metadata$array),
                               sprintf("N samples: %s",BiocGenerics::ncol(object)),
                               sep = "\n"))

#' @describeIn QCmetricsBG number of samples
#' @param x a \code{QCmetricsBG} objects
#' @export
setMethod("ncol", "QCmetricsBG",
          function(x) nrow(x@Sample_Metrics))


#' Combine QCmetricsBG objects
#'
#' Combine QCmetricsBG objects
#' @describeIn QCmetricsBG number of samples
#' @param x \code{QCmetricsBG} object
#' @param y \code{QCmetricsBG} object
#' @param force force combining when the metric-files do not have the same package version, probe_qc_method or detection P-value/beadnr threshold.
#' @param name name
#' @param ... Additional arguments.
#' @export

#'
setMethod("combineQCmetrics", signature = c("QCmetricsBG", "QCmetricsBG"),
          function(x, y, ..., force = FALSE, name = NULL) {
            combine_metrics_bg(x, y, ..., force = force, name = name)
          }
)

