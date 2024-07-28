setGeneric("prepareData",
           function(qcmetrics, samplesheet = NULL, variables = "all")
             standardGeneric("prepareData")
)


setMethod("prepareData", "QCmetrics",
          function(qcmetrics, samplesheet = NULL, variables = "all") {
            .prepareData(qcmetrics = qcmetrics, samplesheet = samplesheet, variables = variables)
          }
)

.prepareData <- function(qcmetrics, samplesheet = NULL, variables = "all") {

  if(is.null(samplesheet)) {
    samplesheet <- qcmetrics@samplesheet
  } else {
    samplesheet <- samplesheet
  }

  check <- all(c("Sample_Name", "Sex") %in% colnames(samplesheet))
  if(!check) {
    stop("'Sample_Name' and/or 'Sex' column are not present in samplesheet!")
  }

  check_sex <- unique(samplesheet$Sex)
  check_sex <- check_sex[!is.na(check_sex)]
  if(length(check_sex) > 2 || !all(check_sex %in% c("F", "M"))) {
    stop("'Sex' column should be coded as 'F'/'M', with missings coded as NA")
  }

  sample_metrics <- qcmetrics@Sample_Metrics

  ## Check if samplesheet and metrics file contain the same samples
  ## If not, stop.
  if(!(all(samplesheet$Sample_Name %in% sample_metrics$Sample_Name))) {
    message("Not all samples in the samplesheet are present in the metrics-file, the intersecting samples are selected.")
    samplesheet <- samplesheet %>% dplyr::filter(Sample_Name %in% sample_metrics$Sample_Name)
  }
  if(!(all(sample_metrics$Sample_Name %in% samplesheet$Sample_Name))) {
    stop("Not all samples in the metrics-file are present in the samplesheet")
  }

  ## Add PCs to sample metrics file (if provided)
  if(!is.null(qcmetrics@PCs)) {
    sample_metrics <- sample_metrics %>% dplyr::left_join(qcmetrics@PCs, by = "Sample_Name")
  }

  if(!is.null(qcmetrics@controlPCs)) {
    sample_metrics <- sample_metrics %>% dplyr::left_join(qcmetrics@controlPCs, by = "Sample_Name")
  }

  ## Add PMS
  sample_metrics <- sample_metrics %>% dplyr::left_join(qcmetrics@pms, by = "Sample_Name")

  if(variables[1] != "all") {
    variables <- variables
  } else {
    variables <- colnames(samplesheet)[colnames(samplesheet) != "Sample_Name"]
    cellcount_variables <- setdiff(colnames(qcmetrics@cellcounts), "Sample_Name")
    variables <- unique(c(variables, cellcount_variables))
  }

  ## All variables (includes PCA outliers)
  variables <- c('None', variables)

  if(!is.null(qcmetrics@ibs$ibs_geno_identical)) {
    outlier_variables <- c("None","MU_outlier", "RG_ratio_outlier", "GR_ratio_outlier",
                           "OP_outlier",  "HC_outlier", "bscon_outlier", "detectionP_outlier",
                           "beadNr_outlier", "sex_outlier",
                           "IBS_outlier", "PCA_outlier", "controlPCA_outlier")
  } else {
    outlier_variables <- c("None","MU_outlier" ,"OP_outlier",
                           "RG_ratio_outlier", "GR_ratio_outlier",
                           "HC_outlier", "bscon_outlier", "detectionP_outlier",
                           "beadNr_outlier", "sex_outlier", "PCA_outlier", "controlPCA_outlier")
  }

  ## Joining will cause a problem when the samplesheet contain column names that are
  ## used in the metrics file or outlier_variables. '.y' will be added to these column names
  tst <- colnames(samplesheet)[colnames(samplesheet) %in% c(colnames(sample_metrics),outlier_variables)]
  tst <- tst[tst != "Sample_Name"]
  if(length(tst) > 0) {
    colnames(samplesheet)[colnames(samplesheet) %in% tst] <- paste0(tst, ".y")
    variables[variables %in% tst] <- paste0(tst, ".y")
    print("Some column names in the samplesheet conflict, '.y' will be added to these column names")
  }

  # Join sample metrics with samplesheet
  data <- sample_metrics %>% dplyr::left_join(samplesheet, by = "Sample_Name")

  # Check if cellcounts are not already in samplesheet..
  sect <- setdiff(colnames(qcmetrics@cellcounts),colnames(samplesheet))
  # Add cellcounts
  data <- data %>% dplyr::left_join(qcmetrics@cellcounts[,c("Sample_Name", sect)], by = "Sample_Name")

  # Return
  list(data = data, samplesheet = samplesheet, outlier_variables = outlier_variables, variables = variables)
}
