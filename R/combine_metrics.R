combine_metrics <- function(..., QCmetricslist = NULL, force = FALSE, include_probe_metrics = TRUE, npcs = 30, name = NULL,
                            PCs = NULL, heatmap_variables = NULL, heatmap_variables_fct = NULL) {

  if(!is.null(QCmetricslist)) {
    objects <- QCmetricslist
  } else {
    objects <- list(...)
  }
  # # If length is 1 this probably means that the user provided the metric-files
  # # in a list -> and thus, unpack it
  # if(length(objects) == 1) {
  #   objects <- unlist(objects)
  # }

  ## Checks
  # Array types
  array <- purrr::map_chr(objects, function(x) slot(x, "metadata")$array)
  if(length(unique(array)) > 1) {
    stop("Different array types cannot be merged")
  } else {
    message(sprintf("Array type: %s", unique(array)))
  }

  # Package versions
  package_version <- purrr::map_chr(objects, function(x) as.character(slot(x, "metadata")$package_version))
  if(length(unique(package_version)) > 1) {
    if(!force) {
      stop("Package versions differ between metric-files!")
    } else {
      message("Not all package versions are identical, will continue merging")
    }
  } else {
    message("[Check] Same package version was used for all objects")
  }

  message("Checking whether thresholds are identical across objects..")

  # Detection P-value
  detp_threshold <- purrr::map_dbl(objects, function(x) slot(x, "detp_threshold"))
  if(length(unique(detp_threshold)) > 1) {
    if(!force) {
      stop("Detection P-value thresholds differ between metric-files!")
    } else {
      message("Not all detection P-value thresholds are identical, will continue merging")
    }
  } else {
    message("[Check] Detection P-value thresholds are identical")
  }

  # Beadnr threshold
  beadnr_threshold <- purrr::map_dbl(objects, function(x) slot(x, "beadnr_threshold"))
  if(length(unique(beadnr_threshold)) > 1) {
    if(!force) {
      stop("Beadnr thresholds differ between metric-files!")
    } else {
      message("Not all detection Beadnr thresholds are identical, will continue merging")
    }
  } else {
    message("[Check] Beadnr thresholds are identical")
  }

  # Probe qc method
  probe_qc_method <-  purrr::map_chr(objects, function(x) slot(x, "probe_qc_method"))
  if(length(unique(probe_qc_method)) > 1) {
    if(!force) {
      stop("Probe qc methods differ between metric-files!")
    } else {
      message(sprintf("Not all probe qc methods are identical, will continue merging and use the first: %s", probe_qc_method[1]))
      probe_qc_method <- probe_qc_method[1]
    }
  } else {
    message("[Check] Probe QC methods are identical")
  }

  ## Sample QC thresholds
  sample_qc_thresholds <- purrr::map(objects, function(x) slot(x, "sample_qc_outliers")$thresholds)
  check <- purrr::map_lgl(sample_qc_thresholds, .f = function(x) identical(sample_qc_thresholds[[1]], x))
  if(!all(check)) {
    message("Sample QC thresholds are not identical across QCmetrics objects, default thresholds will be used")
    if(array[1] == "450k") {
      sample_qc_thresholds <- list(MU = 2000, RG_ratio = 0.5, GR_ratio = 0.5, OP = 11.75, HC = 12, bscon = 80,
                                   detP = 0.05, beadNr = 0.05, IBS_mean = 1.9, IBS_var = 0.1, XY_diff = -2, custom_outliers = NA)
    } else if (array[1] == "EPIC") {
      sample_qc_thresholds <- list(MU = 2000, RG_ratio = 0.4, GR_ratio = 0.4, OP = 12, HC = 12.75, bscon = 80,
                                   detP = 0.05, beadNr = 0.05, IBS_mean = 1.9, IBS_var = 0.1, XY_diff = -2, custom_outliers = NA)
    } else {
      stop("Array type not known!")
    }
  } else {
    sample_qc_thresholds <- sample_qc_thresholds[[1]]
    message("[Check] Sample QC thresholds are identical")
  }

  ## Probe QC thresholds
  probe_qc_thresholds <- purrr::map(objects, function(x) slot(x, "probe_qc_outliers")$thresholds)
  check <- purrr::map_lgl(probe_qc_thresholds, .f = function(x) identical(probe_qc_thresholds[[1]], x))
  if(!all(check)) {
    message("Probe QC thresholds are not identical across QCmetrics objects, default thresholds will be used")
    probe_qc_thresholds <- list(detP = 0.05, beadNr = 0.05)
  } else {
    probe_qc_thresholds <- probe_qc_thresholds[[1]]
    message("[Check] Probe QC thresholds are identical")
  }
  ## Relatedness thresholds
  relatedness_thresholds <- purrr::map(objects, function(x) slot(x, "relatedness_thresholds"))
  check <- purrr::map_lgl(relatedness_thresholds, .f = function(x) identical(relatedness_thresholds[[1]], x))
  if(!all(check)) {
    message("Relatedness thresholds are not identical across QCmetrics objects, default thresholds will be used")
    relatedness_thresholds <- list(IBS_mean = 1.5, IBS_var = 0.35)
  } else {
    relatedness_thresholds <- relatedness_thresholds[[1]]
  }

  ## Combining

  # Number of metric files, samples, males, females
  message("Combining objects..")
  nr_objects <- length(objects)

  ##
  if (length(unique(probe_qc_method)) > 1) {
    message(sprintf("`probe_qc_method` is not identical across QCmetrics, will use the first: %s", probe_qc_method[1]))
    probe_qc_method <- probe_qc_method
  }

  if(unique(probe_qc_method) == "pre_sampleqc") {
    n <- purrr::map_dbl(objects, .f = function(x) BiocGenerics::ncol(x))
    n_female <-  purrr::map_dbl(objects, .f = .get_nr_sex, sex = "F", probe_qc_method = "pre_sampleqc")
    n_male <- purrr::map_dbl(objects, .f = .get_nr_sex, sex = "M", probe_qc_method = "pre_sampleqc")
  } else {
    n <- purrr::map_dbl(objects, .f = function(x) BiocGenerics::ncol(x) - nrow(x@sample_qc_outliers$sample_outliers))
    n_female <- purrr::map_dbl(objects, .f = .get_nr_sex, sex = "F", probe_qc_method = "post_sampleqc")
    n_male <- purrr::map_dbl(objects, .f = .get_nr_sex, sex = "M", probe_qc_method = "post_sampleqc")
  }

  # Combine sample QC
  sample_qc <- purrr::map_df(objects, function(x) slot(x, "Sample_Metrics"))

  # Genotypes (if available), relatedness is dropped (for now)
  geno <- purrr::map_df(objects, .f = function(x) slot(x, "ibs")$ibs_geno_identical)

  # Set to NULL if no data is available
  if(nrow(geno) == 0) geno <- NULL

  # Cellcounts
  cellcounts <- purrr::map_df(objects, .f = function(x) slot(x, "cellcounts"))

  # PMSs
  pms <- purrr::map_df(objects, .f = function(x) slot(x, "pms"))

  # Samplesheets
  samplesheets <- purrr::map_df(objects, .f = function(x) slot(x, "samplesheet"))

  # Control probes
  controlmatrix <- purrr::map(objects, .f = function(x) slot(x, "controlmatrix"))
  cols <- purrr::map(controlmatrix, .f = colnames)
  test <- all(purrr::map_lgl(cols, .f = function(x) identical(x, cols[[1]])))

  message("Recalculating control probe PCA..")
  if(!test) {
    stop("Control probe matrices do not contain the same columns!")
  } else {
    controlmatrix <- do.call(rbind, controlmatrix)
    controlmatrix_scaled <- scale(controlmatrix)
    # Perform PCA
    controlPCs <- prcomp(controlmatrix_scaled)$x

    # If controlmatrix contains fewer samples than npcs, select all PCs
    if(npcs > nrow(controlmatrix)) {
      controlPCs <- controlPCs
    } else {
      controlPCs <- controlPCs[,1:npcs,drop=FALSE]
    }
    samples <- rownames(controlPCs)
    controlPCs <- dplyr::as_tibble(controlPCs)
    colnames(controlPCs) <- paste0(colnames(controlPCs), "_control")
    controlPCs <- controlPCs %>%
      dplyr::mutate(Sample_Name = samples) %>%
      dplyr::select(Sample_Name, dplyr::everything())

    ## Recalculate heatmap P-values (using default values)
    heatmap_pvals_contolPCA <- .get_heatmap_pvals(samplesheet = samplesheets,
                                                  PCs = controlPCs,
                                                  cellcounts = cellcounts,
                                                  pms = pms,
                                                  variables = heatmap_variables,
                                                  factor = heatmap_variables_fct,
                                                  logfile = NULL)
  }


  ## Sample qc outliers
  sample_qc_outliers <- .get_sample_outliers(sample_qc,
                                             thresholds = sample_qc_thresholds,
                                             samplesheet = samplesheets)
  ## Probe qc outliers
  if(include_probe_metrics) {
    probe_qc <- .combine_probeqc(objects, n = n, n_male = n_male, n_female = n_female, array = unique(array))
    if(all(n == 0)) {
      probe_qc$BeadNr_Percentage <-probe_qc$detP_Percentage <- 1
    }
    probe_qc_outliers <- .get_probe_outliers(probe_qc, thresholds = probe_qc_thresholds)
  }

  experimenters <- unique(unlist(purrr::map(objects, .f = function(x) x@metadata$experimenter)))
  if(length(experimenters) > 1) experimenters <- paste(experimenters, collapse = ";")
  metadata <- list(experiment_name = name,
                   experimenter = experimenters,
                   array = unique(array),
                   package_version = if(!force) unique(package_version) else package_version[[1]],
                   analysis_date = Sys.Date()
  )

  metrics <- new("QCmetrics",
                 Sample_Metrics = sample_qc,
                 Probe_Metrics = if(include_probe_metrics) probe_qc else NULL,
                 samplesheet = samplesheets,
                 ibs = list(ibs_geno_identical = geno,
                            ibs_geno_relatedness = NULL,
                            ibs_dnam = NULL,
                            called_probes = NULL,
                            called_probes_geno = NULL,
                            save_all_relationships = FALSE
                 ),
                 PCs = NULL,
                 controlPCs = controlPCs,
                 heatmap_pvals = NULL,
                 heatmap_pvals_controlPCA = heatmap_pvals_contolPCA,
                 cellcounts = cellcounts,
                 controlmatrix = controlmatrix,
                 pms = pms,
                 probe_qc_method = if(!force) unique(probe_qc_method) else probe_qc_method[[1]],
                 sample_qc_outliers = list(sample_outliers = sample_qc_outliers,
                                           thresholds = sample_qc_thresholds),
                 probe_qc_outliers = list(probe_outliers = probe_qc_outliers,
                                          thresholds = probe_qc_thresholds),
                 detp_threshold = if(!force) unique(detp_threshold) else detp_threshold[[1]],
                 beadnr_threshold = if(!force) unique(beadnr_threshold) else beadnr_threshold[[1]],
                 nprobes_PCA = "none",
                 relatedness_thresholds = relatedness_thresholds,
                 selected_samples = NA_character_,
                 metadata = metadata)

  # Adding PCs
  if(!is.null(PCs)) {
    message("Adding PCs")
    metrics <- addPCs(metrics, PCs = PCs, heatmap_variables = heatmap_variables, heatmap_variables_fct = heatmap_variables_fct)
  }

  # Return
  message("Done!")
  metrics
}


combine_metrics_bg <- function(..., force = FALSE, name = NULL) {

  objects <- list(...)
  # If length is 1 this probably means that the user provided the metric-files
  # in a list -> and thus, unpack it
  if(length(objects) == 1) {
    objects <- list(x, y, unlist(objects))
  }

  ## Checks
  # Array types
  array <- purrr::map_chr(objects, function(x) slot(x, "metadata")$array)
  if(length(unique(array)) > 1) {
    stop("Different array types cannot be merged")
  } else {
    message(sprintf("Array type: %s", unique(array)))
  }

  # Package versions
  package_version <- purrr::map_chr(objects, function(x) as.character(slot(x, "metadata")$package_version))
  if(length(unique(package_version)) > 1) {
    if(!force) {
      stop("Package versions differ between metric-files!")
    } else {
      message("Not all package versions are identical, will continue merging")
    }
  } else {
    message("[Check] Same package version was used for all objects")
  }

  ## Combining

  # Number of metric files, samples, males, females
  message("Combining objects..")
  nr_objects <- length(objects)

  # Combine sample QC
  sample_qc <- purrr::map_df(objects, function(x) slot(x, "Sample_Metrics"))

  experimenters <- unique(unlist(purrr::map(objects, .f = function(x) x@metadata$experimenter)))
  if(length(experimenters) > 1) experimenters <- paste(experimenters, collapse = ";")
  metadata <- list(experiment_name = name,
                   experimenter = experimenters,
                   array = unique(array),
                   package_version = if(!force) unique(package_version) else package_version[[1]],
                   analysis_date = Sys.Date()
  )

  metrics <- new("QCmetricsBG",
                 Sample_Metrics = sample_qc,
                 metadata = metadata)
  metrics

}

.get_nr_sex <- function(metrics, sex, probe_qc_method) {
  if(probe_qc_method == "pre_sampleqc") {
    predictedSex <- metrics@Sample_Metrics$predictedSex
    sum(predictedSex == sex)
  } else {
    sample_metrics <- metrics@Sample_Metrics %>%
      dplyr::filter(!Sample_Name %in% metrics@sample_qc_outliers$sample_outliers$Sample_Name)
    sum(sample_metrics$predictedSex == sex)
  }
}

.combine_probeqc <- function(objects, n, n_male, n_female, array) {
  # Extract probe metrics
  probe_qc <- purrr::map(objects, function(x) slot(x,"Probe_Metrics"))
  detp <- purrr::map(probe_qc, .f = function(x) x[,c("Probe", "detP_Percentage")])
  beadnr <- purrr::map(probe_qc, .f = function(x) x[,c("Probe", "BeadNr_Percentage")])
  # probes <- purrr::map(probe_qc, .f = function(x) x[,c("Probe"), drop = FALSE])
  # probes <- Reduce(function(...) dplyr::left_join(..., by='Probe'), probes)
  probe_names <- purrr::map(probe_qc, "Probe")
  # Check if probe names are identical, for EPIC, this is the case for different versions of the array
  test <- mean(purrr::map_dbl(probe_names, .f = function(x) mean(x %in% probe_names[[1]])) == 1)
  # if(test != 1) {
  #   stop("Probe names are not identical!")
  # }

  if(array == "450k") {
    y.probes <- y.probes.450k
  } else if (array == "EPIC") {
    y.probes <- y.probes.EPIC
  }

  # Combine
  detp <- base::Reduce(function(...) dplyr::left_join(..., by='Probe'), detp)
  beadnr <- base::Reduce(function(...) dplyr::left_join(..., by='Probe'), beadnr)
  detp.y <- detp %>% dplyr::filter(Probe %in% y.probes)
  beadnr.y <- beadnr %>% dplyr::filter(Probe %in% y.probes)
  detp <- detp %>% dplyr::filter(!(Probe %in% y.probes))
  beadnr <- beadnr %>% dplyr::filter(!(Probe %in% y.probes))
  detp.probes <- detp$Probe
  beadnr.probes <- beadnr$Probe
  detp.y.probes <- detp.y$Probe
  beadnr.y.probes <- beadnr.y$Probe

  detp <- as.matrix(detp %>% dplyr::select(-Probe))
  beadnr <- as.matrix(beadnr %>% dplyr::select(-Probe))
  detp.y <- as.matrix(detp.y %>% dplyr::select(-Probe))
  beadnr.y <- as.matrix(beadnr.y %>% dplyr::select(-Probe))
  # Set NAs to 1 (completely missing)
  detp[is.na(detp)] <- 1
  beadnr[is.na(beadnr)] <- 1
  detp.y[is.na(detp.y)] <- 1
  beadnr.y[is.na(beadnr.y)] <- 1


  # Fractions
  fraction <- n/sum(n)
  fraction_female <- n_female/sum(n_female)
  fraction_male <- n_male/sum(n_male)
  detp.updated <- detp %*% fraction
  beadnr.updated <- beadnr %*% fraction

  if(all(n_male == 0)) {
    detp.y.updated <- detp.y %*% fraction_female
    beadnr.y.updated <- beadnr.y %*% fraction_female
  } else {
    detp.y.updated <- detp.y %*% fraction_male
    beadnr.y.updated <- beadnr.y %*% fraction_male
  }

  detp.updated <- tibble::tibble(Probe = c(detp.probes, detp.y.probes), detP_Percentage = rbind(detp.updated, detp.y.updated)[,1])
  beadnr.updated <- tibble::tibble(Probe = c(beadnr.probes, beadnr.y.probes), BeadNr_Percentage = rbind(beadnr.updated, beadnr.y.updated)[,1])
  probe_qc <- detp.updated %>% dplyr::left_join(beadnr.updated, by = "Probe")

  ## Return
  probe_qc
}
