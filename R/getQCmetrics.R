#' Create a QC-metrics file
#'
#' Extract QC metrics and outlier information from IDAT files.
#' Returns summarized QC metrics, sample and probe outliers, relatedness info and
#' DNA-methylation-based predictors (PMSs) of various traits/exposures
#' in a \linkS4class{QCmetrics} object, which can be queried interactively using the \code{\link{QCapp}} method.
#'
#'
#' Although many parameters can be tweaked,
#' only the `samplesheet` parameter is mandatory.
#' Note that the samplesheet must contain the following columns: Sample_Name','Sex' (M/F) and 'Basename'.
#' The 'Basename' column should contain paths to the IDAT files.
#' Optionally, 'Age', 'Smoking', 'BMI' and 'Alcohol' columns can be included.
#' If these are provided provided, these can be compared to DNA-methylation-based predictors (PMSs)
#' of the respective variables in the \code{\link{QCapp}}.
#' All other included columns are kept and can be used within the \code{\link{QCapp}} to color the variables by those columns
#' @param samplesheet Samplesheet. Should contain the following columns: Sample_Name','Sex' (M/F) and 'Basename'.
#' @param name Optionally, provide a name for the experiment (will be included in metadata).
#' @param experimenter Optionally, name of the experiment (will be included in metadata).
#' @param chunk_size Size of chunks, default = 40. Larger chunks will be a little bit faster
#'   but will require more RAM.
#' @param geno An optional matrix with genotypes (rows = snps, columns = samples)
#' to check for concordance between DNAm and genotype data. See the tutoriols for more information.
#' @param snplinker Two-column data frame that links DNAm probes to genotype probes (first column = DNAm probes, second column = genotype IDs).
#' If not provided the 65 designated SNP-probes will be used.
#' @param logfile Filepath to write logfile to, by default no logfile is written (`logfile = NULL`).
#' @param verbose If `TRUE` (the default), information is printed in the console.
#' @param probe_qc_method Should probe QC be performed after sample QC ('post_sampleqc') or prior to sample QC ('pre_sampleqc')?.
#' If 'post_sampleqc', sample outliers will be removed before performing site QC.
#' @param heatmap_variables Character vector specifying additional variables to test for association with PCs (default = `NULL`).
#' @param heatmap_variables_fct Relevant if the `heatmap_variables` argument is not `NULL`.
#' Which of the variables specified in `heatmap_variables` should be treated as factors?
#' @param sample_qc_thresholds List of sample QC thresholds, if `NULL` (the default),
#' default thresholds will be used.
#' Note that custom thresholds can be specified for a subset of metrics, default
#' thresholds will then be used for metrics that are not specified.
#' @param probe_qc_thresholds List of probe QC thresholds, if `NULL` (the default),
#' default thresholds will be used.
#' @param detp_threshold detection P-value threshold (default = 1e-16).
#' @param beadnr_threshold beadnr threshold (default = 3).
#' @param relatedness_thresholds Initial relatedness thresholds, can be modified in the \code{\link{QCapp}}.
#' @param custom_outliers Optional additional column name(s) in the samplesheet that
#' that indicate that samples should be treated as outliers.
#' The indicated columns should be of class `logical`, (`TRUE` = outlier).
#' @param epidish_method Which method to use in the EpiDish \link[EpiDISH]{epidish} cellcount prediction (RPC, CBS, or CP).
#' Default = 'RPC'.
#' @param tissue Tissue-type for cellcount-estimation.
#' Currently 'blood', 'breast' and 'epithelium' are implemented, following the \link[EpiDISH]{epidish} function.
#' @param calculate_PCs Should PCA be performed? (TRUE/FALSE).
#' Default = `FALSE`. For large sample sizes RAM could become a bottleneck.
#' To limit RAM usage (besides setting this parameter to `FALSE` to skip PCA), the
#' `n_probes_PCA` argument can be used to limit PCA to a subset of most variable probes.
#' @param n_probes_PCA Number of most variable probes to maintain for PCA,
#' defaults to 'NULL' (all probes are used).
#' @param n_estimate_PCA Relevant if `n_probes_PCA` is not `NULL`.
#'   This argument specifies the number of samples used to estimate the most variable probes (defaults to 'chunk_size').
#' @param nrPCs nr of PCs to retain. Default = 10
#' @param save_all_relationships Save all relatedness data (IBS)?
#' Relatedness output contains N^2 rows, and thus it will rapidly increase in size with large sample sizes.
#' If `FALSE`, relations defined as unrelated (based on the relatedness thresholds) will not be saved.
#' If `NULL` (the default), the parameter will be set to `FALSE` for sample sizes >= 1000, and to `TRUE` otherwise.
#' @param path_failed_measurements Filepath (.rds extension) for saving a matrix that indicates for each measurement whether it failed
#' based on the detP and beadNr cutoffs. Defaults to `NULL`, in which case this matrix is not saved.
#' This option is useful for the \code{\link{mask_measurements}} function,
#' which sets failed measurements in a beta-matrix to missing.
#' @param cores Number of cores to use, defaults to 1.
#' The `BPPARAM` argument controls parallelization (see: \href{https://bioconductor.org/packages/devel/bioc/vignettes/BiocParallel/inst/doc/Introduction_To_BiocParallel.pdf}{BiocParallel}).
#' @param BPPARAM Relevant if the `cores` argument > 1.
#' Specifies the parallel backend to use
#' (see \href{https://bioconductor.org/packages/devel/bioc/vignettes/BiocParallel/inst/doc/Introduction_To_BiocParallel.pdf}{BiocParallel}),
#' defaults to \code{\link[BiocParallel]{MulticoreParam}}
#' @param force See the \link[minfi]{minfi} documentation: Should reading different size IDAT files be forced?
#' Defaults to FALSE.
#' @import IlluminaHumanMethylation450kanno.ilmn12.hg19
#' @import IlluminaHumanMethylationEPICanno.ilm10b2.hg19
#' @importFrom graphics text
#' @importFrom stats kruskal.test lm p.adjust prcomp sd var wilcox.test
#' @importFrom utils capture.output sessionInfo
#' @return A \linkS4class{QCmetrics} object.
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
#' # Add some custom thresholds:
#' # Note that you can specify only the thresholds that you want to change,
#' # for the others the default values will be used
#'  thresholds <- list(MU = 2000, OP = 13, RG_ratio = 0.6)
#'  metrics <- getQCmetrics(samplesheet,
#'   name = "example minfiData",
#'   calculate_PCs = TRUE,
#'   n_probes_PCA = 20000,
#'   sample_qc_thresholds = thresholds,
#'   heatmap_variables = c("Array", "Slide"),
#'   heatmap_variables_fct = c("Array", "Slide"))
#' @export

getQCmetrics <- function(samplesheet,
                          name = NULL,
                          experimenter = NULL,
                          chunk_size = 40,
                          geno = NULL,
                          snplinker = NULL,
                          logfile = NULL,
                          verbose = TRUE,
                          probe_qc_method = c("post_sampleqc","pre_sampleqc"),
                          heatmap_variables = NULL,
                          heatmap_variables_fct = NULL,
                          sample_qc_thresholds = NULL,
                          probe_qc_thresholds = NULL,
                          detp_threshold = 1e-16,
                          beadnr_threshold = 3,
                          relatedness_thresholds = list(IBS_mean = 1.50, IBS_var = 0.35),
                          custom_outliers = NULL,
                          epidish_method = c("RPC", "CBS", "CP"),
                          tissue = c("blood", "breast", "epithelium"),
                          calculate_PCs = FALSE,
                          n_probes_PCA = NULL,
                          n_estimate_PCA = chunk_size,
                          nrPCs = 10,
                          save_all_relationships = NULL,
                          path_failed_measurements = NULL,
                          cores = 1,
                          BPPARAM = NULL,
                          force = FALSE
                         ) {

  ## Logfile
  if(!is.null(logfile)) file.create(logfile)
  if(!is.null(logfile)) cat(sprintf("Creating a QCmetrics object\nAnalysis started at: %s\n", Sys.time()),
                            file = logfile, append = TRUE)
  if(verbose) message(sprintf("Creating a QCmetrics object\nAnalysis started at: %s", Sys.time()))

  # If 'Sample_Name' and 'Sex' are not present in samplesheet, stop.
  check <- all(c("Sample_Name", "Sex") %in% colnames(samplesheet))
  if(!check) {
    stop("'Sample_Name' and/or 'Sex' column are not present in samplesheet!")
  }

  if(is.null(save_all_relationships) && nrow(samplesheet) >= 1000) {
    save_all_relationships <- FALSE
  } else if(is.null(save_all_relationships) && nrow(samplesheet) < 1000)  {
    save_all_relationships <- TRUE
  }

  # Chunksize
  stopifnot(chunk_size > 1)

  samplesheet$Sample_Name <- as.character(samplesheet$Sample_Name)

  check_sex <- unique(samplesheet$Sex)
  check_sex <- check_sex[!is.na(check_sex)]
  if(length(check_sex) > 2 || !all(check_sex %in% c("F", "M"))) {
    stop("'Sex' column should be coded as 'F'/'M', with missings coded as NA")
  }

  # Some input checks
  probe_qc_method <- match.arg(probe_qc_method)
  epidish_method <- match.arg(epidish_method)
  tissue <- match.arg(tissue)

  if(!is.null(heatmap_variables)) {
    # Check if variables are present in samplesheet
    check <- mean(heatmap_variables %in% colnames(samplesheet))
    if(check < 1) {
      stop("Heatmap variables that are not present in samplesheet are specified!")
    }
  }

  # Load one sample to check which array is used
  rgset <- minfi::read.metharray.exp(targets = data.frame(samplesheet[1,]),
                                     verbose = FALSE,
                                     extended = FALSE)
  array <- .guessArrayTypes(nrow(rgset))
  anno <- data.frame(minfi::getAnnotation(rgset)) %>%
    dplyr::select(chr, pos, strand, Name, AddressA, AddressB,Type,Color)

  # SNPlinker file
  if(!is.null(snplinker)) {
    if(ncol(snplinker) > 2) {
      # Assume that Probe and SNPID are in the first two columns
      snplinker <- snplinker[,1:2]
      colnames(snplinker) <- c("Probe", "snpID")
    } else if (ncol(snplinker) == 2) {
      colnames(snplinker) <- c("Probe", "snpID")
    } else {
      colnames(snplinker) <- "Probe"
    }
  } else {
    # Default snplinker table: designated SNP probes
    snpprobes <- minfi::getSnpBeta(rgset)
    snplinker <- tibble::tibble(Probe = rownames(snpprobes), snpID = rownames(snpprobes))
  }

  if(!is.null(logfile)) cat(sprintf("Array: %s\n-----------------------------\n", array), file = logfile, append = TRUE)
  if(verbose) message(sprintf("Array: %s\n-----------------------------", array))

  # Load default thresholds
  if(array == "450k") {
    thresholds_default <- list(MU = 2000, RG_ratio = 0.5, GR_ratio = 0.5, OP = 11.75, HC = 12, bscon = 80,
                               detP = 0.05, beadNr = 0.05, IBS_mean = 1.9, IBS_var = 0.1, XY_diff = -2, custom_outliers = NA)
    thresholds_probe_default <- list(detP = 0.05, beadNr = 0.05)
  } else if (array == "EPIC") {
    thresholds_default <- list(MU = 2000, RG_ratio = 0.4, GR_ratio = 0.4, OP = 12, HC = 12.75, bscon = 80,
                               detP = 0.05, beadNr = 0.05, IBS_mean = 1.9, IBS_var = 0.1, XY_diff = -2, custom_outliers = NA)
    thresholds_probe_default <- list(detP = 0.05, beadNr = 0.05)
  } else {
    stop("Array type not known!")
  }

  ## Thresholds, if thresholds are provided these will be used.
  # If user-provided thresholds are missing for some metrics, these will be filled in with default thresholds
  if(!is.null(sample_qc_thresholds)) {
    check <- names(thresholds_default) %in% names(sample_qc_thresholds)
    thresholds_not_provided <- names(thresholds_default)[!check]
    if(length(thresholds_not_provided) > 0) {
      if(!is.null(logfile)) {cat("\nNot all sample QC thresholds provided, defaults will be used for the following thresholds:\n", file = logfile,
          append = TRUE)}
      if(verbose) message("Not all sample QC thresholds provided, defaults will be used for the following thresholds:")
      for(i in thresholds_not_provided) {
        if(!is.null(logfile)) cat(sprintf("%s\n", i), file = logfile, append = TRUE)
        if(verbose)message(sprintf("%s", i))
      }
      thresholds_ <- c(sample_qc_thresholds, thresholds_default[thresholds_not_provided])
    } else {
      thresholds_ <- sample_qc_thresholds
    }
  } else {
    thresholds_ <- thresholds_default
  }

  if(!is.null(probe_qc_thresholds)) {
    check <- names(thresholds_probe_default) %in% names(probe_qc_thresholds)
    thresholds_not_provided <- names(thresholds_probe_default)[!check]
    if(length(thresholds_not_provided) > 0) {
      if(!is.null(logfile)) {cat("\nNot all probe QC thresholds provided, defaults will be used for the following thresholds:\n", file = logfile,
                                 append = TRUE)}
      if(verbose) message("Not all probe QC thresholds provided, defaults will be used for the following thresholds:")
      for(i in thresholds_not_provided) {
        if(!is.null(logfile)) cat(sprintf("%s\n", i), file = logfile, append = TRUE)
        if(verbose)message(sprintf("%s", i))
      }
      probe_thresholds_ <- c(probe_qc_thresholds, thresholds_probe_default[thresholds_not_provided])
    } else {
      probe_thresholds_ <- probe_qc_thresholds
    }
  } else {
    probe_thresholds_ <- thresholds_probe_default
  }

  ## Custom outliers
  if(!is.null(custom_outliers)) {
    if(!all(custom_outliers %in% colnames(samplesheet))) stop("Custom outliers should be present in samplesheet!")
    if(!all(purrr::map_lgl(samplesheet[,custom_outliers, drop = FALSE], .f = is.logical))) stop("Custom outliers should be logical (TRUE/FALSE)!")
    if(sum(is.na(samplesheet[,custom_outliers, drop = FALSE])) > 0) stop("Custom outliers should not contain missings!")
    thresholds_[["custom_outliers"]] <- custom_outliers
  }

  experiment_name <- name

  ## Extract metrics in chunks -------------------------------------------------------------

  chunks <- split(samplesheet$Sample_Name, ceiling(seq_along(samplesheet$Sample_Name)/chunk_size))
  # Check if any of the chunk sizes = 1, if so, append to previous chunk
  if(any(lengths(chunks) == 1)) {
    length_1 <- which(lengths(chunks) == 1)
    chunks[[length_1 - 1]] <- c(chunks[[length_1 - 1]], chunks[[length_1]])
    chunks[[length_1]] <- NULL
  }

  # n_probes_PCA is specified, estimate variances in specified number of samples
  if(calculate_PCs && !is.null(n_probes_PCA) && n_estimate_PCA < nrow(samplesheet)) {
    # Randomly sample from samplesheet
    set.seed(10)
    samples_estimate <- sample(samplesheet$Sample_Name, size = n_estimate_PCA)
    rgset <- minfi::read.metharray.exp(targets = data.frame(samplesheet %>% dplyr::filter(Sample_Name %in% samples_estimate)),
                                       verbose = FALSE,
                                       extended = FALSE,
                                       force = force
                                       )
    beta <- minfi::getBeta(rgset)
    beta <- beta[!rownames(beta) %in% c(xy.probes.450k, xy.probes.EPIC),,drop=FALSE]
    var <- matrixStats::rowVars(beta, na.rm = TRUE)
    names(var) <- rownames(beta)
    if(n_probes_PCA <= nrow(beta)) {
      probes_PCA_keep <- names(sort(var, decreasing=TRUE)[1:n_probes_PCA])
    } else {
      probes_PCA_keep <- NULL
    }
  } else {
    probes_PCA_keep <- NULL
  }

  if(!is.null(logfile)) cat(sprintf("Extracting QCmetrics in chunks of %s samples..\n", chunk_size),
                            file = logfile, append = TRUE)
  if(verbose) message(sprintf("Extracting QCmetrics in chunks of %s samples..", chunk_size))

  ## Run
  if(cores > 1) {
    if(is.null(BPPARAM)) BPPARAM <- BiocParallel::MulticoreParam(workers = cores)
    outputs <- BiocParallel::bplapply(chunks,
                          FUN = get_qcmetrics_chunk,
                          samplesheet = samplesheet,
                          all_chunks = chunks,
                          verbose = verbose,
                          logfile = logfile,
                          anno = anno,
                          geno = geno,
                          snplinker = snplinker,
                          thresholds = thresholds_,
                          probe_qc_method = probe_qc_method,
                          sample_qc_thresholds = thresholds_,
                          detp_threshold = detp_threshold,
                          beadnr_threshold = beadnr_threshold,
                          epidish_method = epidish_method,
                          tissue = tissue,
                          calculate_PCs = calculate_PCs,
                          probes_PCA = probes_PCA_keep,
                          BPPARAM = BPPARAM,
                          force = force
                          )
  } else {
    ## Previously used lapply/map here
    ## However, the bulk of the working memory is used by probe (detp/beadnr) matrix
    ## These are stored in a sparseMatrix, separate sparseMatrices use way more memory
    ## then combining them. It is therefore more memory-efficient to combine them along
    ## the way in a for loop
    outputs <- list()
    # Good 'ol for loop
    for(i in 1:length(chunks)) {
      output.tmp <- get_qcmetrics_chunk(chunks[[i]],
                                        all_chunks = chunks,
                                        samplesheet = samplesheet,
                                        verbose = verbose,
                                        logfile = logfile,
                                        anno = anno,
                                        geno = geno,
                                        snplinker = snplinker,
                                        thresholds = thresholds_,
                                        probe_qc_method = probe_qc_method,
                                        sample_qc_thresholds = thresholds_,
                                        detp_threshold = detp_threshold,
                                        beadnr_threshold = beadnr_threshold,
                                        epidish_method = epidish_method,
                                        tissue = tissue,
                                        calculate_PCs = calculate_PCs,
                                        probes_PCA = probes_PCA_keep,
                                        force = force
                                        )
      if(i == 1) {
        probe_qc_metrics <- output.tmp$probe_metrics
      } else {
        if (!force) {
          probe_qc_metrics <- cbind(probe_qc_metrics, output.tmp$probe_metrics)
        } else {
          intersec <- intersect(rownames(probe_qc_metrics), rownames(output.tmp$probe_metrics))
          probe_qc_metrics <- cbind(probe_qc_metrics[intersec,,drop=FALSE], output.tmp$probe_metrics[intersec,,drop=FALSE])
        }

      }
      # Set to NULL
      output.tmp$probe_metrics <- NULL
      outputs[[i]] <- output.tmp
    }
  }

  ## Extract metrics in chunks -------------------------------------------------------------

  if(!is.null(logfile)) cat("Done!\nMerging the chunks..\n", file = logfile, append = TRUE)
  if(verbose) message("Done!\nMerging the chunks..\n")

   # extract info
   arrays <- purrr::map_chr(outputs, "array")
   if(dplyr::n_distinct(arrays) > 1) {
     stop("Multiple array types detected, we recommend processing them separately.")
   } else {
     array <- unique(arrays)
   }
   methylaid_metrics <- dplyr::bind_rows(purrr::map(outputs, "methylaid_metrics"))
   bscon <- dplyr::bind_rows(purrr::map(outputs, "bscon"))
   sex <- dplyr::bind_rows(purrr::map(outputs, "sex"))
   rg <- dplyr::bind_rows(purrr::map(outputs, "rg"))
   cellcounts <- dplyr::bind_rows(purrr::map(outputs, "cellcounts"))
   pms <- dplyr::bind_rows(purrr::map(outputs, "pms"))
   controlMatrix <- do.call(rbind,(purrr::map(outputs, "controlMatrix")))
   betasnps <- do.call(cbind, (purrr::map(outputs, "betasnps")))

   if(calculate_PCs) {
     if (force) {
       intersec <- purrr::map(outputs, function(x) rownames(x[["beta"]]))
       intersec <- Reduce(intersect, intersec)
       beta <- do.call(cbind, purrr::map(outputs, function(x) x[["beta"]][intersec,,drop=FALSE]))
     } else {
       beta <- do.call(cbind, purrr::map(outputs, "beta"))
     }

     # calculate variances if  n_estimate_PCA >= nrow(samplesheet)
     if(!is.null(n_probes_PCA) && n_estimate_PCA >= nrow(samplesheet)) {
       beta <- beta[!rownames(beta) %in% c(xy.probes.450k, xy.probes.EPIC),,drop=FALSE]
       var <- matrixStats::rowVars(beta)
       names(var) <- rownames(beta)
       if(n_probes_PCA <= nrow(beta)) {
         probes_PCA_keep <- names(sort(var, decreasing=TRUE)[1:n_probes_PCA])
         beta <- beta[as.character(probes_PCA_keep),,drop=FALSE]
       } else {
         probes_PCA_keep <- NULL
       }
     }
   }

   ## IBS ----------------------------------------------------
   if(!is.null(logfile)) cat("Performing IBS..", file = logfile, append = TRUE)
   if(verbose) message("Performing IBS..")

   # Check overlap between DNAm data and genotype matrix
   if(!is.null(geno)) {
     overlap_geno_beta <- intersect(colnames(betasnps), colnames(geno))
     if(length(overlap_geno_beta) == 0) {
       if(!is.null(logfile)) cat("There is no overlap between the genotype matrix and the DNAm data!
           Genotype concordance will not be calculated.\n",
                                 file = logfile, append = TRUE)
       if(verbose) message("There is no overlap between the genotype matrix and the DNAm data!
               Genotype concordance will not be calculated.")

       geno <- NULL
     }
   }
   ibs <- .get_ibs(beta = betasnps, geno = geno, snplinker = snplinker, verbose = verbose)
   ibs$save_all_relationships <- save_all_relationships

   if(!save_all_relationships) {
    ibs <- .summarize_relationships(
     ibs = ibs,
     relatedness_thresholds = relatedness_thresholds)
   } else {
     ibs$hexagons_relatedness_dnam <- NULL
     ibs$hexagons_relatedness_geno <- NULL
   }

   ## Merge probe metrics -------------------------------------------------------------
   if(cores > 1) {
     probe_qc_metrics <- do.call(cbind, purrr::map(outputs, "probe_metrics"))
   }
   # If path_failed_measurements = specified, save failed measurements:
   if(!is.null(path_failed_measurements)) {
     if(!is.null(logfile)) cat(sprintf("\nSaving a matrix specifying failed measurements to %s", path_failed_measurements), file = logfile, append = TRUE)
     if(verbose) message(sprintf("\nSaving a matrix specifying failed measurements to %s", path_failed_measurements))
     saveRDS(as(probe_qc_metrics, "lgCMatrix"),path_failed_measurements)
     }
   probe_qc_metrics <- .get_probe_metrics(probe_qc_metrics, probe_qc_method = probe_qc_method,
                                          anno = anno, sex = sex, bscon = bscon, ibs = ibs,
                                          methylaid_metrics = methylaid_metrics, rg = rg,
                                          samplesheet = samplesheet, thresholds = thresholds_)

   ## PCA
   if(calculate_PCs) {
     if(!is.null(logfile)) cat("\nPerforming array-wide PCA..", file = logfile, append = TRUE)
     if(verbose) message("Performing array-wide PCA..")
     PCs <- .get_PCs(beta, nrPCs = nrPCs)
   }

   ## Perform control probe PCA ----------------------------------------------------
   if(!is.null(logfile)) cat("\nPerforming control probe PCA..", file = logfile, append = TRUE)
   if(verbose) message("Performing control probe PCA..")
   control_PCs <- .control_PCA(controlMatrix)

   ## Associations PCs and technical/biological variables ---------------------------

   # Array-wide PCs
   if(calculate_PCs) {
     if(!is.null(logfile)) cat("\nPerforming association tests between array-wide PCs and technical/biological variables..", file = logfile, append = TRUE)
     if(verbose) message("Performing association tests between array-wide PCs and technical/biological variables..")
     heatmap_pvals <- .get_heatmap_pvals(samplesheet = samplesheet,
                                         PCs = PCs,
                                         cellcounts = cellcounts,
                                         pms = pms,
                                         variables = heatmap_variables,
                                         factor = heatmap_variables_fct,
                                         logfile = NULL)
   }

   if(!is.null(logfile)) cat("\nPerforming association tests between control probe PCs and technical/biological variables..\n", file = logfile, append = TRUE)
   if(verbose) message("Performing association tests between control probe PCs and technical/biological variables..")
   heatmap_pvals_contolPCA <- .get_heatmap_pvals(samplesheet = samplesheet,
                                                 PCs = control_PCs,
                                                 cellcounts = cellcounts,
                                                 pms = pms,
                                                 variables = heatmap_variables,
                                                 factor = heatmap_variables_fct,
                                                 logfile = NULL)

  ## Merge QC measures ---------------------------
  if(!is.null(logfile)) cat("\nMerging the metrics..\n", file = logfile, append = TRUE)
  if(verbose) message("\nMerging the metrics..")
  sample_qc_metrics <- .join_sample_metrics(methylaid_metrics = methylaid_metrics,
                                            bscon = bscon,
                                            sex = sex,
                                            rg = rg,
                                            detp = probe_qc_metrics$detp_samples,
                                            beadnr = probe_qc_metrics$beadnr_samples,
                                            ibs = ibs)

  probe_qc_metrics <- probe_qc_metrics$beadnr_probes %>% dplyr::left_join(probe_qc_metrics$detp_probes, by = "Probe")
  ## Check number of outliers using default thresholds:
  sample_outliers <- .get_sample_outliers(sample_qc_metrics, thresholds_, samplesheet = samplesheet)
  outliers_percentage <- (nrow(sample_outliers) / nrow(sample_qc_metrics)) * 100

  # Write number of outliers to logfile
  if(!is.null(logfile)) cat(sprintf("\nNumber of outliers using %s thresholds:\n", if(!is.null(sample_qc_thresholds)) "user-defined" else "default"), file = logfile, append = TRUE)
  if(verbose) message(sprintf("\nNumber of outliers using %s thresholds:\n", if(!is.null(sample_qc_thresholds)) "user-defined" else "default"))
  for(i in colnames(sample_outliers)[colnames(sample_outliers) != "Sample_Name"]) {
    name <- stringr::str_replace(i, pattern = "_outlier", replacement = "")
    nr.outliers <- if(nrow(sample_outliers) > 0) sum(sample_outliers[[i]]) else 0
    percentage.outliers <- round((nr.outliers*100)/nrow(sample_qc_metrics),2)
    if(percentage.outliers > 10) {
      text <- sprintf("%s: %s outliers (%s%%) [WARNING] High percentage of outliers!",name, nr.outliers, percentage.outliers )
    } else {
      text <- sprintf("%s: %s outliers (%s%%)", name, nr.outliers, percentage.outliers)
    }
    if(verbose) message(text)
    if(!is.null(logfile)) cat(paste0(text, "\n"), file = logfile, append = TRUE)
  }

  # Probe outliers
  probe_outliers <- .get_probe_outliers(probe_qc_metrics,
                                        thresholds = probe_thresholds_)

  metadata <- list(experiment_name = experiment_name,
                   experimenter = experimenter,
                   array = array,
                   package_version = packageVersion("QCpG"),
                   analysis_date = Sys.Date()
                   )

  metrics <- new("QCmetrics",
                  Sample_Metrics = sample_qc_metrics,
                  Probe_Metrics = probe_qc_metrics,
                  samplesheet = samplesheet,
                  ibs = ibs,
                  PCs = if(calculate_PCs) PCs else NULL,
                  controlPCs = control_PCs,
                  heatmap_pvals = if(!calculate_PCs) NULL else heatmap_pvals,
                  heatmap_pvals_controlPCA = heatmap_pvals_contolPCA,
                  cellcounts = cellcounts,
                  controlmatrix = controlMatrix,
                  pms = pms,
                  probe_qc_method = probe_qc_method,
                  sample_qc_outliers = list(sample_outliers = sample_outliers, thresholds = thresholds_),
                  probe_qc_outliers = list(probe_outliers = probe_outliers, thresholds = probe_thresholds_),
                  detp_threshold = detp_threshold,
                  beadnr_threshold = beadnr_threshold,
                  nprobes_PCA = if(calculate_PCs && !is.null(probes_PCA_keep)) {as.character(nrow(beta))} else if(calculate_PCs && is.null(probes_PCA_keep)) {"all"} else {"none"},
                  relatedness_thresholds = relatedness_thresholds,
                  selected_samples = NA_character_,
                  metadata = metadata)

  ## Final info ---------------------------
  if(verbose) message("Done!\n")
  if(verbose) message(sprintf("Analysis finished at: %s\n\n", Sys.time()))

  if(!is.null(logfile)) {
    cat("Done!\n\n", file = logfile, append = TRUE)
    cat(sprintf("Analysis finished at: %s\n\n", Sys.time()), file = logfile, append = TRUE)
    cat("SessionInfo:\n\n",file = logfile, append = TRUE)
    readr::write_lines(capture.output(sessionInfo()), path = logfile, append = TRUE)
  }

  # Return metrics
  metrics
}

get_qcmetrics_chunk <- function(samples,
                                all_chunks,
                                samplesheet,
                                anno = NULL,
                                verbose = TRUE,
                                logfile = NULL,
                                geno = NULL,
                                snplinker = NULL,
                                thresholds = NULL,
                                probe_qc_method,
                                sample_qc_thresholds,
                                detp_threshold = 1e-16, beadnr_threshold = 3,
                                epidish_method = c("RPC", "CBS", "CP"),
                                tissue = c("blood", "breast", "epithelium"),
                                calculate_PCs = TRUE,
                                probes_PCA = NULL,
                                force = FALSE
                                ) {

  # Load iDATs
  samplesheet <- samplesheet %>% dplyr::filter(Sample_Name %in% samples)
  rgset <- minfi::read.metharray.exp(targets = data.frame(samplesheet),
                              verbose = FALSE,
                              extended = TRUE,
                              force = force
                              )

  # Guess array type
  array <- .guessArrayTypes(nrow(rgset))

  # Betas
  beta <- minfi::getBeta(rgset)

  # Snpprobes
  snpprobes <- minfi::getSnpBeta(rgset)

  ## Bscon metric ---------------------------
  bscon <- .get_bscon(rgset, array = array)

  ## Cellcounts ---------------------------
  cellcounts <- .get_cellcounts(beta, epidish_method = epidish_method, tissue = tissue)

  pms <- .get_pms(beta, verbose = if(any(samples %in% all_chunks[[1]])) verbose else FALSE, logfile = if(any(samples %in% all_chunks[[1]])) logfile else NULL)

  ## Median Red/Green type I signals  --------------------------
  rg <- .get_medianRG(rgset = rgset, anno = anno)

  ## Sex check ---------------------------
  sex <- suppressWarnings(.get_sex(rgset, cutoff = sample_qc_thresholds[["XY_diff"]]))

  ## MethylAid --------------------------
  ## The following code is adapted from the methylaid github (https://github.com/bbmri-nl/MethylAid/tree/master/R)
  methylaid_metrics <- .get_MethylAid(rgset)

  ## Extract genotype info ---------------------------
  if(!is.null(snplinker)) {
    betasnps <- rbind(beta, snpprobes)
    snplinker_ <- snplinker %>% dplyr::filter(Probe %in% rownames(betasnps))
    betasnps <- betasnps[as.character(snplinker_$Probe),, drop = FALSE]
  } else {
    # Default snplinker table: designated SNP probes
    betasnps <- snpprobes
  }

  ## Extract controlProbes  ---------------------------
  controlMatrix <- .get_control_probes(rgset)

  ## Extract detp and beadnr ---------------------------
  # probe_metrics <- .get_detp_beadnr(rgset, anno = anno, detp_threshold = detp_threshold, beadnr_threshold = beadnr_threshold)
  # probe_metrics <- .get_probe_metrics(rgset = rgset, probe_qc_method = probe_qc_method, anno = anno,
  #                                     sex = sex, detp_threshold = detp_threshold, beadnr_threshold = beadnr_threshold,
  #                                     bscon = bscon, ibs = NULL, methylaid_metrics = methylaid_metrics, rg = rg,
  #                                     samplesheet = samplesheet, thresholds = thresholds)
  probe_metrics <- .get_probe_matrix(rgset = rgset, anno = anno,
                                       detp_threshold = detp_threshold,
                                       beadnr_threshold = beadnr_threshold)

  output <- list(methylaid_metrics = methylaid_metrics, bscon = bscon, sex = sex,
                 rg = rg, probe_metrics = probe_metrics, cellcounts = cellcounts,
                 pms = pms,
                 controlMatrix = controlMatrix, betasnps = betasnps,
                 array = array,
                 beta = if(calculate_PCs && !is.null(probes_PCA)) {beta[as.character(probes_PCA),,drop=FALSE]} else if (calculate_PCs) {beta} else {NULL})
  # Return output
  output
}

.summarize_relationships <- function(ibs, relatedness_thresholds) {

  hexagons_relatedness_dnam <- hexbin::hexbin(ibs$ibs_dnam %>% dplyr::filter(IBS_mean < relatedness_thresholds$IBS_mean | IBS_var > relatedness_thresholds$IBS_var) %$% IBS_mean,
                                              ibs$ibs_dnam %>% dplyr::filter(IBS_mean < relatedness_thresholds$IBS_mean | IBS_var > relatedness_thresholds$IBS_var) %$% IBS_var,
                                              xbins = 30,
                                              xbnds=c(min(ibs$ibs_dnam$IBS_mean),2), ybnds=c(0,max(ibs$ibs_dnam$IBS_var)),
                                              IDs = TRUE)

  hexagons_relatedness_dnam <- data.frame(hexbin::hcell2xy(hexagons_relatedness_dnam),
                                          cell = hexagons_relatedness_dnam@cell,
                                          count = hexagons_relatedness_dnam@count)
  ibs$hexagons_relatedness_dnam <- hexagons_relatedness_dnam
  ibs$ibs_dnam <- ibs$ibs_dnam %>%
    dplyr::filter(IBS_mean > relatedness_thresholds$IBS_mean & IBS_var < relatedness_thresholds$IBS_var)

  if(!is.null(ibs$ibs_geno_relatedness)) {
    hexagons_relatedness_geno <- hexbin::hexbin(ibs$ibs_geno_relatedness %>% dplyr::filter(IBS_mean < relatedness_thresholds$IBS_mean | IBS_var > relatedness_thresholds$IBS_var) %$% IBS_mean,
                                                ibs$ibs_geno_relatedness %>% dplyr::filter(IBS_mean < relatedness_thresholds$IBS_mean | IBS_var > relatedness_thresholds$IBS_var) %$% IBS_var,
                                                xbins = 30,
                                                xbnds=c(min(ibs$ibs_geno_relatedness$IBS_mean),2), ybnds=c(0,max(ibs$ibs_geno_relatedness$IBS_var)),
                                                IDs = TRUE)

    hexagons_relatedness_geno <- data.frame(hexbin::hcell2xy(hexagons_relatedness_geno),
                                            cell = hexagons_relatedness_geno@cell,
                                            count = hexagons_relatedness_geno@count)
    ibs$hexagons_relatedness_geno <- hexagons_relatedness_geno
    ibs$ibs_geno_relatedness <- ibs$ibs_geno_relatedness %>%
      dplyr::filter(IBS_mean > relatedness_thresholds$IBS_mean & IBS_var < relatedness_thresholds$IBS_var)
  } else {
    ibs$hexagons_relatedness_geno <- NULL
  }
  return(ibs)
}

