## Update QCmetrics file with normalized beta-values,
# this will recalculate predicted cell type fractions, age, smokingscore, alcohol and BMI.
# + optionally calculated new PCs based on these beta-values

update_function <- function(qcmetrics,
                            beta,
                            samplesheet = NULL,
                            beta_format = NULL,
                            chunk_size = 500,
                            read_chunks = TRUE,
                            delimiter = "\t",
                            calculate_PCs = FALSE,
                            n_probes_PCA = NULL,
                            n_estimate_PCA = chunk_size,
                            nrPCs = 10,
                            epidish_method = c("RPC", "CBS", "CP"),
                            heatmap_variables = NULL,
                            heatmap_variables_fct = NULL,
                            verbose = TRUE
                            ) {

  epidish_method <- match.arg(epidish_method)

  if(is.null(samplesheet)) samplesheet <- getSamplesheet(qcmetrics)
  check <- all(c("Sample_Name", "Sex") %in% colnames(samplesheet))
  if(!check) {
    stop("'Sample_Name' and/or 'Sex' column are not present in samplesheet!")
  }

  check_sex <- unique(samplesheet$Sex)
  check_sex <- check_sex[!is.na(check_sex)]
  if(length(check_sex) > 2 || !all(check_sex %in% c("F", "M"))) {
    stop("'Sex' column should be coded as 'F'/'M', with missings coded as NA")
  }

  # Check type of beta argument
  if(is.character(beta) && is.null(beta_format))  {
    extension <- tools::file_ext(beta)
    if(!extension %in% c("txt", "rds")) {
      stop("A character string is passed to the `beta` argument, but the extension is not
           'rds' or 'txt'. If the file is an .rds or .txt file, but does the filepath does not
           end with '.rds' or '.txt' respectively, set the `beta_format` argument to 'rds' or 'txt'.
           ")
    } else {
      beta_format <- extension
    }
  }

  # Check tissue
  if("CD4T" %in% colnames(qcmetrics@cellcounts)) {
    tissue <- "blood"
  } else if ("Fat" %in% colnames(qcmetrics@cellcounts)) {
    tissue <- "breast"
  } else {tissue <- "epithelium"}

  if(is.character(beta) && !is.null(beta_format) && !beta_format %in% c("txt", "rds")) {
    stop("`beta_format` argument must be either 'txt' or 'rds'")
  }

  if(is.character(beta) && beta_format == "rds") {
    message(sprintf("Loading %s..\n", beta))
    beta <- readRDS(beta)
    samples <- colnames(beta)
    probes <- rownames(beta)
  } else if (is.character(beta)) {
    samples <- read.table(file=beta,
                          sep = delimiter,
                          nrows=1, header=TRUE,
                          check.names=FALSE)
    samples <- colnames(samples)[2:length(colnames(samples))]
  } else {
    samples <- colnames(beta)
    probes <- rownames(beta)
  }

  ## Check for overlap between QCmetrics object and beta-matrix
  overlapping_samples <- intersect(qcmetrics@samplesheet$Sample_Name, samples)
  if(length(overlapping_samples) == 0) {
    stop("There are no overlapping samples between the QCmetrics object and the beta-matrix")
  } else {
    message(sprintf("%s samples overlap between the QCmetrics object and the beta-matrix",
                    length(overlapping_samples)))
    message(sprintf("%s samples will be removed from the QC metrics object",
                    ncol(qcmetrics) - length(overlapping_samples)))
    samplesheet <- samplesheet %>% dplyr::filter(Sample_Name %in% overlapping_samples)
  }


  if(is.character(beta) && beta_format == "txt" && read_chunks && length(samples) > chunk_size) {
    message(sprintf("Reading %s in chunks of %s samples..\n", beta, chunk_size))
    used_chunks <- TRUE
    chunks <- split(samples, ceiling(seq_along(samples)/chunk_size))

    # Check if any of the chunk sizes = 1
    if(any(lengths(chunks) == 1)) {
      length_1 <- which(lengths(chunks) == 1)
      chunks[[length_1 - 1]] <- c(chunks[[length_1 - 1]], chunks[[length_1]])
      chunks[[length_1]] <- NULL
    }

    ## PCA
    if(calculate_PCs && !is.null(n_probes_PCA) && n_estimate_PCA < nrow(samplesheet)) {
      # Randomly sample from samplesheet
      set.seed(10)
      samples_estimate <- sample(samples, size = n_estimate_PCA)
      colclasses <- c("character", ifelse(samples %in% samples_estimate, "double", "NULL"))
      beta_ <- read.table(file=beta,
                         sep = delimiter,
                         header=TRUE,
                         check.names=FALSE,
                         colClasses=colclasses)
      probes <- beta_[[1]]
      beta_[[1]] <- NULL
      beta_ <- as.matrix(beta_)
      rownames(beta_) <- probes
      beta_ <- beta_[!rownames(beta_) %in% c(xy.probes.450k, xy.probes.EPIC),,drop=FALSE]
      var <- matrixStats::rowVars(beta_, na.rm = TRUE)
      names(var) <- rownames(beta_)
      if(n_probes_PCA <= nrow(beta_)) {
        probes_PCA_keep <- names(sort(var, decreasing=TRUE)[1:n_probes_PCA])
      } else {
        probes_PCA_keep <- NULL
      }
    } else {
      probes_PCA_keep <- NULL
    }

    results_chunks <- purrr::map(
      chunks,
      .f = .update_qcmetrics_chunk,
      path = beta,
      delimiter = delimiter,
      all_samples = samples,
      verbose = verbose,
      logfile = NULL, #might implement later
      epidish_method = epidish_method,
      calculate_PCs = calculate_PCs,
      probes_PCA = probes_PCA_keep,
      tissue = tissue
    )
    probes <- results_chunks[[1]]$probes
    cellcounts <- purrr::map_df(results_chunks, "cellcounts")
    pms <- purrr::map_df(results_chunks, "pms")

  } else {
    used_chunks <- FALSE
    if(is.character(beta) && beta_format == "txt") {
      # beta <- readr::read_delim(file=beta,
      #            delim = delimiter,
      #            col_names = TRUE)
      colclasses <- c("character", rep("double", length(samples)))
      beta <- read.table(file = beta,
                         sep = delimiter,
                         header = TRUE,
                         check.names=FALSE,
                         colClasses = colclasses
                         )
      probes <- beta[[1]]
      beta[[1]] <- NULL
      beta <- as.matrix(beta)
      rownames(beta) <- probes
    }

    ## Cellcounts ---------------------------

    cellcounts <- .get_cellcounts(beta,
                                  epidish_method = epidish_method,
                                  tissue = tissue)

    ## PMSs
    pms <- .get_pms(beta, verbose = verbose, logfile = NULL)
  }

  if(calculate_PCs) {
    if(used_chunks) {
      beta <- do.call(cbind, purrr::map(results_chunks, "beta"))
    }

    if(!used_chunks && !is.null(n_probes_PCA)) {
      beta <- beta[!rownames(beta) %in% c(xy.probes.450k, xy.probes.EPIC),,drop=FALSE]
      var <- matrixStats::rowVars(beta)
      names(var) <- rownames(beta)
      if(n_probes_PCA <= nrow(beta)) {
        probes_PCA_keep <- names(sort(var, decreasing=TRUE)[1:n_probes_PCA])
        beta <- beta[as.character(probes_PCA_keep),,drop=FALSE]
      } else {
        probes_PCA_keep <- NULL
      }
    } else {
      probes_PCA_keep <- NULL
    }

    PCs <- .get_PCs(beta, nrPCs = nrPCs)
    heatmap_pvals <- .get_heatmap_pvals(samplesheet = samplesheet,
                                        PCs = PCs,
                                        cellcounts = cellcounts,
                                        pms = pms,
                                        variables = heatmap_variables,
                                        factor = heatmap_variables_fct,
                                        logfile = NULL)
  }

  ## Recalculate control probe matrix
  controlmatrix <- qcmetrics@controlmatrix[samples,, drop = FALSE]
  controlmatrix_scaled <- scale(controlmatrix)

  # Perform PCA
  controlPCs <- prcomp(controlmatrix_scaled)$x

  # If rgset contains fewer samples than npcs, select all PCs
  if(nrPCs > nrow(controlmatrix)) {
    controlPCs <- controlPCs
  } else {
    controlPCs <- controlPCs[,1:nrPCs,drop=FALSE]
  }
  samples <- rownames(controlPCs)
  controlPCs <- dplyr::as_tibble(controlPCs)
  colnames(controlPCs) <- paste0(colnames(controlPCs), "_control")
  controlPCs <- controlPCs %>%
    dplyr::mutate(Sample_Name = samples) %>%
    dplyr::select(Sample_Name, dplyr::everything())

  heatmap_pvals_contolPCA <- .get_heatmap_pvals(samplesheet = samplesheet,
                                                PCs = controlPCs,
                                                cellcounts = cellcounts,
                                                pms = pms,
                                                variables = heatmap_variables,
                                                factor = heatmap_variables_fct,
                                                logfile = NULL)

  # Subset ibs
  ibs <- qcmetrics@ibs
  if(!is.null(ibs$ibs_geno_identical)) {
    ibs$ibs_geno_identical <- ibs$ibs_geno_identical %>% dplyr::filter(Sample_Name %in% overlapping_samples)
  }
  if(!is.null(ibs$ibs_geno_relatedness)) {
    ibs$ibs_geno_relatedness <- ibs$ibs_geno_relatedness %>%
      dplyr::filter(Sample_Name %in% overlapping_samples, Sample_Name_geno %in% overlapping_samples)
  }
  if(!is.null(ibs$ibs_dnam)) {
    ibs$ibs_dnam <- ibs$ibs_dnam %>% dplyr::filter(Sample_Name.x %in% overlapping_samples,
                                                   Sample_Name.y %in% overlapping_samples)
  }

  qcmetrics@metadata$analysis_date <- Sys.Date()
  qcmetrics@sample_qc_outliers$sample_outliers <- qcmetrics@sample_qc_outliers$sample_outliers %>%
    dplyr::filter(Sample_Name %in% overlapping_samples)
  qcmetrics@probe_qc_outliers$probe_outliers <- qcmetrics@probe_qc_outliers$probe_outliers %>%
    dplyr::filter(Probe %in% probes)

  ## New QCmetrics object
  metrics <- new("QCmetrics",
                 Sample_Metrics = qcmetrics@Sample_Metrics %>% dplyr::filter(Sample_Name %in% overlapping_samples),
                 Probe_Metrics = qcmetrics@Probe_Metrics %>% dplyr::filter(Probe %in% probes),
                 samplesheet = samplesheet %>% dplyr::filter(Sample_Name %in% overlapping_samples),
                 ibs = ibs,
                 PCs = if(calculate_PCs) PCs else NULL,
                 controlPCs = controlPCs,
                 heatmap_pvals = if(calculate_PCs) heatmap_pvals else NULL,
                 heatmap_pvals_controlPCA = heatmap_pvals_contolPCA,
                 cellcounts = cellcounts,
                 controlmatrix = controlmatrix,
                 pms = pms,
                 probe_qc_method = qcmetrics@probe_qc_method,
                 sample_qc_outliers = qcmetrics@sample_qc_outliers,
                 probe_qc_outliers = qcmetrics@probe_qc_outliers,
                 detp_threshold = qcmetrics@detp_threshold,
                 beadnr_threshold = qcmetrics@beadnr_threshold,
                 nprobes_PCA = if(calculate_PCs && !is.null(probes_PCA_keep)) {as.character(nrow(beta))} else if(calculate_PCs && is.null(probes_PCA_keep)) {"all"} else {"none"},
                 relatedness_thresholds = qcmetrics@relatedness_thresholds,
                 selected_samples = NA_character_,
                 metadata = qcmetrics@metadata)
}


.update_qcmetrics_chunk <- function(samples, path, delimiter, all_samples, verbose, logfile, epidish_method, calculate_PCs, probes_PCA, tissue) {

  ## Read chunk
  colclasses <- c("character", ifelse(all_samples %in% samples, "double", "NULL"))
  beta <- read.table(file=path,
                     sep = delimiter,
                     header=TRUE,
                     check.names=FALSE,
                     colClasses=colclasses)
  probes <- beta[[1]]
  beta[[1]] <- NULL
  beta <- as.matrix(beta)
  rownames(beta) <- probes

  ## Cellcounts
  cellcounts <- .get_cellcounts(beta, epidish_method, tissue = tissue)

  ## PMSs
  pms <- .get_pms(beta, verbose = if(all_samples[1] %in% samples) verbose else FALSE, logfile = if(all_samples[1] %in% samples) logfile else NULL)

  list(cellcounts = cellcounts,
       pms = pms,
       beta = if(calculate_PCs && !is.null(probes_PCA)) {beta[as.character(probes_PCA),,drop=FALSE]} else if (calculate_PCs) {beta} else {NULL},
       probes = if(all_samples[1] %in% samples) probes else NULL)
}
