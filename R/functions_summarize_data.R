## Function adapted from minfi github (https://github.com/hansenlab/minfi/blob/master/R/read.meth.R)
## To infer whether data is EPIC/450k/27k
.guessArrayTypes <- function(nProbes) {
    if (nProbes >= 622000 && nProbes <= 623000) {
            array = "450k"
    } else if (nProbes >= 1050000 && nProbes <= 1053000) {
        # NOTE: "Current EPIC scan type"
            array = "EPIC"
    } else if (nProbes >= 1032000 && nProbes <= 1033000) {
        # NOTE: "Old EPIC scan type"
            array = "EPIC"
    } else if (nProbes >= 54000 && nProbes <= 56000) {
            array = "27k"
    } else {
        array = "Unknown"
        array
    }
}

## Load annotation for the array used
.get_annotation <- function(array) {
  if(array == "450k") {
      anno <- data.frame(minfi::getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19::IlluminaHumanMethylation450kanno.ilmn12.hg19))
      anno
  } else if (array == "EPIC"){
      anno <- data.frame(minfi::getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19::IlluminaHumanMethylationEPICanno.ilm10b2.hg19))
      anno
  } else {
      stop("Array type not known")
  }
}

## bscon metric from wateRmelon package
# Fix for EPIC data -> solved in the latest wateRmelon package (github)
# .get_bscon <- function(rgset, array) {
#   if(array == "EPIC") {
#     coldata <- SummarizedExperiment::colData(rgset)
#     barcodes <- coldata$Basename
#     path <- paste0("/",paste(stringr::str_split(barcodes, pattern = "/", simplify=TRUE)[1,2:7], collapse = "/"))
#     mset <- wateRmelon::readEPIC(path)
#     bscon <- wateRmelon::bscon(mset)
#     rm(mset); gc()
#     names(bscon) <- stringr::str_split(names(bscon), pattern = "/", simplify = TRUE)[,2]
#     bscon <- tibble::tibble(Sample_Name = names(bscon), bscon = bscon)
#     bscon
#     } else {
#         bscon <- wateRmelon::bscon(rgset)
#         bscon <- tibble::tibble(Sample_Name = names(bscon), bscon = bscon)
#         bscon
#     }
# }
#update: problem is now fixed in Bioconductor version of wateRmelon
.get_bscon <- function(rgset, array) {
  bscon <- wateRmelon::bscon(rgset)
  bscon <- tibble::tibble(Sample_Name = names(bscon), bscon = unname(bscon))
  bscon
}

# From minfi
.get_control_probes <- function(rgset) {
  extractedData <- .extractFromRGSet450k(rgset)
  controlMatrix <- .buildControlMatrix450k(extractedData)
  controlMatrix
}

.control_PCA <- function(controlMatrix, npcs = 30) {
  controlMatrix_scaled <- scale(controlMatrix)
  controlPCs <- prcomp(controlMatrix_scaled)$x
  # If rgset contains fewer samples than npcs, select all PCs
  if(npcs > nrow(controlMatrix)) {
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

  controlPCs
}

## Calculate median red/green signals in type I probes
.get_medianRG <- function(rgset, anno) {
  anno_typeI_red <- anno %>% dplyr::filter(Type == "I", Color == "Red")
  anno_typeI_green <- anno %>% dplyr::filter(Type == "I", Color == "Grn")

  # Red signals
  red <- minfi::getRed(rgset)
  red_addressA <- red[rownames(red) %in% anno_typeI_red$AddressA,]
  red_addressB <- red[rownames(red) %in% anno_typeI_red$AddressB,]
  red_total <- (red_addressA + red_addressB)/2
  medians_red <- matrixStats::colMedians(red_total)
  rm(red_addressA, red_addressB, red_total); gc()

  # Green signals
  green <- minfi::getGreen(rgset)
  green_addressA <- green[rownames(green) %in% anno_typeI_green$AddressA,]
  green_addressB <- green[rownames(green) %in% anno_typeI_green$AddressB,]
  green_total <- (green_addressA + green_addressB)/2
  medians_green <- matrixStats::colMedians(green_total)
  rm(green_addressA, green_addressB, green_total); gc()

  rg <- tibble::tibble(Sample_Name = colnames(rgset), Green = medians_green,
               Red = medians_red) %>% dplyr::mutate(RG_ratio = Red/Green, GR_ratio = Green/Red)
  rg
}

## Calculate ibs based on DNAm-inferred SNPs and (optionally) genotype data
.get_ibs <- function(beta, geno = NULL, snplinker = NULL, logfile = NULL, verbose = NULL) {

  # If genotype matrix is provided, concordance between DNAm-inferred SNPs and
  # genotype SNPs will be calculated. Otherwise, only relatedness will be calculated
  # using the DNAm-inferred SNPs
  if(!is.null(geno)) {

    # Check overlapping samples
    overlap_samples <- intersect(colnames(geno), colnames(beta))
    non_overlapping_samples <- setdiff(colnames(beta), colnames(geno))
    geno <- geno[,as.character(overlap_samples), drop = FALSE]

    if(!is.null(logfile)) {cat(sprintf("%s SNPs out of %s SNPs provided in the snplinker table are present in the genotype matrix\n", sum(rownames(geno) %in% snplinker$snpID), nrow(snplinker)),
                               file = logfile, append = TRUE)}
    if(verbose) message(sprintf("%s SNPs out of %s SNPs provided in the snplinker table are present in the genotype matrix", sum(rownames(geno) %in% snplinker$snpID), nrow(snplinker)))
    if(!is.null(logfile)) {cat(sprintf("%s out of %s probes provided in the snplinker table are present in the DNAm data\n", sum(rownames(beta) %in% snplinker$Probe), nrow(snplinker)),
                               file = logfile, append = TRUE)}
    if(verbose) message(sprintf("%s out of %s probes provided in the snplinker table are present in the DNAm data", sum(rownames(beta) %in% snplinker$Probe), nrow(snplinker)))

    # Make sure that probes in snplinker file are present in beta-matrix
    snplinker <- snplinker %>% dplyr::filter(Probe %in% rownames(beta))

    # Call genotypes from DNAm data
    betasnps <- beta[as.character(snplinker$Probe),, drop = FALSE]
    dnamCalls <- omicsPrint::beta2genotype(betasnps)
    called_probes <- rownames(dnamCalls)
    snplinker_called <- snplinker %>% dplyr::filter(Probe %in% called_probes)
    if(!is.null(logfile)) cat(sprintf("%s out of %s probes could be accurately called\n", nrow(dnamCalls), nrow(betasnps)), file = logfile, append = TRUE)
    if(verbose) message(sprintf("%s out of %s probes could be accurately called\n", nrow(dnamCalls), nrow(betasnps)))

    if(length(intersect(snplinker_called$snpID, rownames(geno))) == 0) {
      if(!is.null(logfile)) {cat("There is no overlap between the called SNPs and the genotype data. Might be different SNP IDs?\n",
                                 file = logfile, append = TRUE)}
      if(verbose) message("There is no overlap between the snplinker table and the genotype data. Might be different SNP IDs?")
      if(!is.null(logfile)) {cat("Since no overlapping genotype data is available, only relatedness will be calculated. Using the designated SNP-probes\n",
                                 file = logfile, append = TRUE)}
      if(verbose) message("Since no overlapping genotype data is available, only relatedness will be calculated")

      if(!is.null(logfile)) cat("Calculating relatedness using only the DNAm-inferred SNPs..\n",file = logfile, append = TRUE)
      if(verbose) message("Calculating relatedness using only the DNAm-inferred SNPs..")

      # Calculate IBS
      ibs_dnam <- omicsPrint::alleleSharing(dnamCalls,verbose = FALSE)
      ibs_dnam <- ibs_dnam %>% dplyr::filter(relation == "unrelated")

      ## Output of allelesharing contains factors, update to character
      ibs_dnam <- ibs_dnam %>%
        dplyr::mutate(
          colnames.x = as.character(colnames.x),
          colnames.y = as.character(colnames.y)
        ) %>%
        dplyr::rename(Sample_Name.x = colnames.x, Sample_Name.y = colnames.y,
                      IBS_mean = mean, IBS_var = var) %>%
        dplyr::select(-relation)

      list(ibs_geno_identical = NULL, ibs_geno_relatedness = NULL,
           ibs_dnam = ibs_dnam, called_probes = called_probes, called_probes_geno = NULL)

    } else {
      ## Calculate relatedness using DNAm data
      if(!is.null(logfile)) cat("Calculating relatedness using only the DNAm-inferred SNPs..\n",file = logfile, append = TRUE)
      if(verbose) message("Calculating relatedness using only the DNAm-inferred SNPs..")
      ibs_dnam <- omicsPrint::alleleSharing(dnamCalls, verbose = FALSE)
      ibs_dnam <- ibs_dnam %>% dplyr::filter(relation == "unrelated")

      # Call snps
      snplinker_called <- snplinker_called %>% dplyr::filter(snpID %in% rownames(geno))
      if(!is.null(logfile)) {cat(sprintf("%s out of %s called probes are present in genotype matrix\n",
                                         sum(snplinker_called$Probe %in% rownames(dnamCalls)), nrow(dnamCalls)),
                                 file = logfile, append = TRUE)}
      if(verbose) {message(sprintf("%s out of %s called probes are present in genotype matrix",
                                   sum(snplinker_called$Probe %in% rownames(dnamCalls)), nrow(dnamCalls)))}

      # Create final geno and beta matrix
      dnamCalls <- dnamCalls[as.character(snplinker_called$Probe),, drop = FALSE]
      geno <- geno[as.character(snplinker_called$snpID),, drop = FALSE]
      rownames(dnamCalls) <- rownames(geno)

      # Relabel SNP matrix (omicsPrint uses 1,2,3 instead of 0,1,2)
      geno_relabeled <- ifelse(geno == 2, 3, ifelse(
        geno == 1, 2, ifelse(geno == 0, 1, NA)))

      # Match
      if(!is.null(logfile)) cat("Calculating genotype concordance..\n",file = logfile, append = TRUE)
      if(verbose) message("Calculating genotype concordance..")
      ibs <- omicsPrint::alleleSharing(dnamCalls, geno_relabeled, verbose = FALSE, alignment = TRUE) # Alignment is important!!

      ibs_dnam <- ibs_dnam %>%
        dplyr::mutate(
          colnames.x = as.character(colnames.x),
          colnames.y = as.character(colnames.y)
        ) %>%
        dplyr::rename(Sample_Name.x = colnames.x, Sample_Name.y = colnames.y,
                      IBS_mean = mean, IBS_var = var) %>%
        dplyr::select(-relation)

      ## Output of allelesharing contains factors, update to character
      ibs <- ibs %>%
        dplyr::mutate(
          colnames.x = as.character(colnames.x),
          colnames.y = as.character(colnames.y)
        ) %>%
        dplyr::rename(Sample_Name = colnames.x,
                      IBS_mean = mean, IBS_var = var)

      list(ibs_geno_identical = ibs %>% dplyr::filter(relation == "identical") %>% dplyr::select(-c(colnames.y, relation)),
           ibs_geno_relatedness = ibs %>% dplyr::filter(relation == "unrelated") %>%
             dplyr::rename(Sample_Name_geno = colnames.y) %>% dplyr::select(-relation),
           ibs_dnam = ibs_dnam, called_probes = called_probes, called_probes_geno = snplinker_called$Probe)
    }

  } else {
    # Only calculate DNAm relatedness

    betasnps <- beta[as.character(snplinker$Probe),, drop = FALSE ]
    dnamCalls <- omicsPrint::beta2genotype(betasnps)
    called_probes <- rownames(dnamCalls)
    if(!is.null(logfile)) cat(sprintf("%s out of %s probes could be accurately called\n", nrow(dnamCalls), nrow(betasnps)), file = logfile, append = TRUE)
    if(verbose) message(sprintf("%s out of %s probes could be accurately called", nrow(dnamCalls), nrow(betasnps)))

    if(!is.null(logfile)) cat("Calculating relatedness using the DNAm-inferred SNPs\n",file = logfile, append = TRUE)
    if(verbose) message("Calculating relatedness using the DNAm-inferred SNPs")
    ibs_dnam  <- omicsPrint::alleleSharing(dnamCalls,verbose = FALSE)
    ibs_dnam <- ibs_dnam %>% dplyr::filter(relation == "unrelated")
    ## Output of allelesharing contains factors, update to character
    ibs_dnam <- ibs_dnam %>%
      dplyr::mutate(
        colnames.x = as.character(colnames.x),
        colnames.y = as.character(colnames.y)
      ) %>%
      dplyr::rename(Sample_Name.x = colnames.x, Sample_Name.y = colnames.y,
                    IBS_mean = mean, IBS_var = var) %>%
      dplyr::select(-relation)

    list(ibs_geno_identical = NULL, ibs_geno_relatedness = NULL,
         ibs_dnam = ibs_dnam, called_probes = called_probes, called_probes_geno = NULL)
  }

}

.merge_probe_metrics <- function(detp, beadnr, n, n_male, n_female, array) {
  if(array == "450k") {
    y.probes <- y.probes.450k
  } else if (array == "EPIC") {
    y.probes <- y.probes.EPIC
  }
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

  ## Update! -> what if fraction males is 0
  # Two options:
  # Set to 1 -> thus will be removed
  # Calculate based on female values (how it was done previous versions)

  # For now: calculate in females, probes that have values < 0.05 are probably cross-reactive probes
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

## Based on code from methyylaid github! https://github.com/bbmri-nl/MethylAid
.get_MethylAid <- function(rgset) {

  rgset <- .replaceZero(rgset)
  ## Get control probe information
  TypeControl <- minfi::getProbeInfo(rgset, type = "Control")
  TypeControl <- as.data.frame(TypeControl)
  colnames(TypeControl)

  R <- minfi::getRed(rgset)
  G <- minfi::getGreen(rgset)

  ##maybe notify when controls are not on the array
  id <- intersect(TypeControl$Address, rownames(R))
  R <- R[rownames(R) %in% id,]
  G <- G[rownames(G) %in% id,]
  TypeControl <- TypeControl[TypeControl$Address %in% id,]
  TypeControl <- TypeControl[order(TypeControl$Address), ]

  R <- log2(R)
  G <- log2(G)
  TypeControl <- TypeControl[!(TypeControl$Type %in%
  c("NORM_A", "NORM_G", "NORM_C", "NORM_T")), ]
  data <- data.frame(Address=rep(rownames(R), ncol(R)),
                     Samples=rep(colnames(R), each=nrow(R)),
                     IntRed=as.vector(R),
                     IntGrn=as.vector(G))
  sdata <- merge(TypeControl, data)

  ## Median MU intensities
  MU <- matrix(0.0, nrow=2, ncol=ncol(rgset))
  rgset <- minfi::preprocessRaw(rgset)
  M <- minfi::getMeth(rgset)
  MU[1,] <- matrixStats::colMedians(M, na.rm=TRUE)
  U <- minfi::getUnmeth(rgset)
  MU[2,] <- matrixStats::colMedians(U, na.rm=TRUE)
  colnames(MU) <- colnames(rgset)
  rownames(MU) <- c("Methylated", "Unmethylated")
  MU <- data.frame(t(MU))
  MU$Sample_Name <- rownames(MU)
  MU <- dplyr::as_tibble(MU)
  MU <- MU %>% dplyr::select(Sample_Name, dplyr::everything())

  qcProbes = list(
    BSI = "^BISULFITE CONVERSION I$",
    BSII = "^BISULFITE CONVERSION II$",
    SPI = "^SPECIFICITY I$",
    SPII = "^SPECIFICITY II$",
    NP = "^NON-POLYMORPHIC$",
    NC = "^NEGATIVE$",
    SC = "^STAINING$",
    EC = "^EXTENSION$",
    TR = "^TARGET REMOVAL$",
    HYB = "^HYBRIDIZATION$"
  ) ## we don't use the normalization controls NORM_A, NORM_G, NORM_C or NORM_T


  ## Non-polymorphic controls
  rotateData <- function(data, columns) {
    data[,columns] <- c(0.5*(data[,columns[1]] + data[,columns[2]]),
                        data[,columns[1]] - data[,columns[2]])
    data
  }

    data <- sdata

    d <- data[grepl(qcProbes["NP"], data$Type),]

    dGrn <- d[d$ExtendedType %in% c("NP (C)", "NP (G)"), c(1:5,7)]
    OP_x <- tapply(dGrn$IntGrn, dGrn$Samples, mean)

    dRed <- d[d$ExtendedType %in% c("NP (A)", "NP (T)"), c(1:6)]
    OP_y <- tapply(dRed$IntRed, dRed$Samples, mean)

    data <- data.frame(OP_x, OP_y)
    data <- rotateData(data, columns=c("OP_x", "OP_y"))
    OP <- data
    samples <- rownames(OP)
    OP <- dplyr::as_tibble(OP)
    OP$Sample_Name <- samples

    ## Bisulfite conversion controls
    # data <- sdata
    # d <- data[grepl(qcProbes["BSI"], data$Type),]
    #
    # dGrn <- d[grepl("C1|C2|C3", d$ExtendedType), c(1:5,7)]
    # BS_x <- tapply(dGrn$IntGrn, dGrn$Samples, mean)
    #
    # dRed <- d[grepl("C4|C5|C6", d$ExtendedType), c(1:6)]  ##EPIC is missing I-C6 Bisulfite control probe and corresponding I-U6
    # BS_y <- tapply(dRed$IntRed, dRed$Samples, mean)
    #
    # data <- data.frame(BS_x, BS_y)
    #
    # data <- rotateData(data, columns=c("BS_x", "BS_y"))
    # BS <- data
    # samples <- rownames(BS)
    # BS <- dplyr::as_tibble(BS)
    # BS$Sample_Name <- samples

    ## Hybridization controls
    data <- sdata
    d <- data[grepl(qcProbes["HYB"], data$Type),]
    d <- d[order(d$Samples),]
    HC_x <- 0.5*(d$IntGrn[grepl("High", d$ExtendedType)] + d$IntGrn[grepl("Low", d$ExtendedType)])
    HC_y <- d$IntGrn[grepl("High", d$ExtendedType)] - d$IntGrn[grepl("Low", d$ExtendedType)]

    data <- data.frame(HC_x, HC_y, row.names=d$Samples[grepl("High", d$ExtendedType)])
    HC <- data
    samples <- rownames(HC)
    HC <- dplyr::as_tibble(HC)
    HC$Sample_Name <- samples

    methylaid_metrics <- MU %>% dplyr::left_join(OP, by = "Sample_Name") %>%
                                dplyr::left_join(HC, by = "Sample_Name")

}

.get_PCs <- function(beta, nrPCs = nrPCs) {
  # Remove xy-probes
  xy_probes <- unique(c(xy.probes.450k, xy.probes.EPIC))
  beta <- beta[!(rownames(beta) %in% xy_probes),]
  # Impute missing values
  invisible(capture.output(beta <- impute::impute.knn(beta))) # Hide output by inpute.knn
  beta <- beta$data

  nrPCs <- if(nrPCs > ncol(beta)) ncol(beta) else nrPCs
  frac <- nrPCs/ncol(beta)

  if(frac < 0.5) {
      PCs <- irlba::prcomp_irlba(t(beta), n = nrPCs)
      PCs <- PCs$x
      PCs <- dplyr::as_tibble(PCs)
      PCs <- PCs %>% dplyr::mutate(Sample_Name = colnames(beta)) %>% dplyr::select(Sample_Name, dplyr::everything())
      PCs
  } else {
      PCs <- prcomp(t(beta))
      PCs <- PCs$x
      PCs <- PCs[,1:nrPCs]
      PCs <- dplyr::as_tibble(PCs)
      PCs <- PCs %>% dplyr::mutate(Sample_Name = colnames(beta)) %>% dplyr::select(Sample_Name, dplyr::everything())
      PCs
  }
}

.get_sex <- function(rgset, cutoff) {
  sex <- .getSex(object = minfi::mapToGenome(rgset), cutoff = cutoff)
  # samples <- rownames(sex)
  # sex <- dplyr::as_tibble(as.data.frame(sex))
  sex <- sex %>% dplyr::select(Sample_Name, dplyr::everything()) %>%
                 dplyr::mutate(XYdiff = yMed - xMed)
  sex
}

.join_sample_metrics <- function(methylaid_metrics, bscon, sex, rg, detp, beadnr, ibs) {

  sample_qc_metrics <- methylaid_metrics %>%
                  dplyr::left_join(rg, by = "Sample_Name") %>%
                  dplyr::left_join(bscon, by = "Sample_Name") %>%
                  dplyr::left_join(beadnr, by = "Sample_Name") %>%
                  dplyr::left_join(detp, by = "Sample_Name") %>%
                  dplyr::left_join(sex, by = "Sample_Name")
   if(!is.null(ibs$ibs_geno_identical)) {
     sample_qc_metrics <- sample_qc_metrics %>% dplyr::left_join(ibs$ibs_geno_identical[,c("Sample_Name", "IBS_mean", "IBS_var")],
                                                          by = c("Sample_Name"))
   }
   sample_qc_metrics
}

.get_sample_outliers <- function(sample_metrics, thresholds, samplesheet) {

  sample_outliers <- sample_metrics %>%
                        dplyr::left_join(samplesheet[,c("Sample_Name", "Sex")], by = "Sample_Name") %>%
                        dplyr::mutate(
                              MU_outlier = Methylated < thresholds[["MU"]] | Unmethylated < thresholds[["MU"]],
                              RG_ratio_outlier = RG_ratio < thresholds[["RG_ratio"]],
                              GR_ratio_outlier = GR_ratio < thresholds[["GR_ratio"]],
                              OP_outlier = ifelse(is.na(OP_x), FALSE, ifelse(OP_x < thresholds[["OP"]], TRUE, FALSE)),
                              HC_outlier = HC_x < thresholds[["HC"]],
                              bscon_outlier = bscon < thresholds[["bscon"]],
                              detectionP_outlier = detP_Percentage > thresholds[["detP"]],
                              beadNr_outlier = BeadNr_Percentage > thresholds[["beadNr"]],
                              sex_outlier = (XYdiff < thresholds[["XY_diff"]] & Sex == "M") | (XYdiff > thresholds[["XY_diff"]] & Sex == "F"),
                              sex_outlier = ifelse(is.na(sex_outlier), TRUE, sex_outlier)
                             )
   if(all(!is.na(thresholds[["custom_outliers"]]))) {
     sample_outliers <- sample_outliers %>%
       dplyr::left_join(samplesheet[,c("Sample_Name", thresholds[["custom_outliers"]])], by = "Sample_Name")
   }

   base_vars <- c("Sample_Name", "MU_outlier",
                  "RG_ratio_outlier", "GR_ratio_outlier",
                  "OP_outlier", "HC_outlier",
                  "bscon_outlier", "detectionP_outlier", "beadNr_outlier", "sex_outlier")

   if(all(c("IBS_mean", "IBS_var") %in% colnames(sample_metrics))) {
     sample_outliers <- sample_outliers %>% dplyr::mutate(IBS_outlier = IBS_mean < thresholds[["IBS_mean"]] | IBS_var > thresholds[["IBS_var"]])
     sample_outliers <- sample_outliers %>% dplyr::mutate(IBS_outlier = ifelse(is.na(IBS_outlier), FALSE, IBS_outlier))

     if(all(!is.na(thresholds[["custom_outliers"]]))) {
       sample_outliers <- sample_outliers %>%
         dplyr::select(c(base_vars, "IBS_outlier",
                         thresholds[["custom_outliers"]]))
     } else {
       sample_outliers <- sample_outliers %>%
         dplyr::select(c(base_vars, "IBS_outlier"))
     }

   } else {

     if(all(!is.na(thresholds[["custom_outliers"]]))) {
       sample_outliers <- sample_outliers %>%
         dplyr::select(c(base_vars, thresholds[["custom_outliers"]]))
     } else {
       sample_outliers <- sample_outliers %>%
         dplyr::select(base_vars)
     }
   }

   ## Create 'Total' column
   outlier_table_tmp <- sample_outliers %>% dplyr::select(-Sample_Name)
   total <- apply(outlier_table_tmp, 1, function(x) any(x))
   sample_outliers$Total <- total
   sample_outliers %>% dplyr::filter(Total)
}

.get_probe_outliers <- function(probe_metrics, thresholds) {

  probe_metrics <- probe_metrics %>%
    dplyr::mutate(BeadNr_Percentage = BeadNr_Percentage > thresholds$beadNr,
                  detP_Percentage = detP_Percentage > thresholds$detP,
                  Total = BeadNr_Percentage | detP_Percentage
                  ) %>%
    dplyr::filter(Total)
  probe_metrics
}

.replaceZero <- function(rgset) {
    ##replace zero intensity values of all probes in the Red and Green
    ##channels with NA's
    Red <- minfi::getRed(rgset)
    Green <- minfi::getGreen(rgset)

    ##replace the zeros
    Red[which(Red == 0)] <- NA
    Green[which(Green == 0)] <- NA

    minfi::RGChannelSet(Green = Green,
                 Red = Red,
                 colData = SummarizedExperiment::colData(rgset),
                 annotation = minfi::annotation(rgset))
}

.get_probe_matrix <- function(rgset, anno, detp_threshold,
                              beadnr_threshold) {

  ## detp
  #detp <- minfi::detectionP(rgset)
  detp <- .detectionP_fix(rgset)

  ## Calculate beadnr
  beadnr <- minfi::getNBeads(rgset)
  anno_typeI <- anno %>% dplyr::filter(Type == "I")
  anno_typeI <- anno_typeI %>% dplyr::filter(AddressA %in% rownames(beadnr) & AddressB %in% rownames(beadnr))
  anno_typeII <- anno %>% dplyr::filter(Type == "II")
  anno_typeII <- anno_typeII %>% dplyr::filter(AddressA %in% rownames(beadnr))
  beadnr_M_typeI <- beadnr[anno_typeI$AddressA,]
  beadnr_U_typeI <- beadnr[anno_typeI$AddressB,]
  beadnr_typeII <- beadnr[anno_typeII$AddressA,]

  # Take minimum beadnr for type II probes
  beadnr_typeI_combined <- lapply(1:ncol(beadnr_M_typeI), function(n) pmin(beadnr_M_typeI[,n], beadnr_U_typeI[,n]))
  beadnr_typeI_combined <- do.call(cbind, beadnr_typeI_combined)
  colnames(beadnr_typeI_combined) <- colnames(beadnr_M_typeI)
  rownames(beadnr_typeII) <- anno_typeII$Name
  rownames(beadnr_typeI_combined) <- anno_typeI$Name
  beadnr <- rbind(beadnr_typeII, beadnr_typeI_combined)
  beadnr <- beadnr[match(rownames(detp), rownames(beadnr)),]

  # Create a sparseMatrix
  mtrx <- detp
  mtrx[detp > detp_threshold & beadnr < beadnr_threshold] <- 3 # both fail
  mtrx[detp > detp_threshold & beadnr >= beadnr_threshold] <- 2 # detp fail
  mtrx[detp < detp_threshold & beadnr < beadnr_threshold] <- 1 # beadnr fail
  mtrx[detp < detp_threshold & beadnr >= beadnr_threshold] <- 0 # no fail
  mtrx <- as(mtrx, "sparseMatrix")

  # Return
  mtrx
}

.get_probe_metrics <- function(mtrx, probe_qc_method, anno, sex, bscon = NULL,
                               ibs = NULL, methylaid_metrics = NULL, rg = NULL, samplesheet = NULL,
                               thresholds = NULL) {
  # y.probes, will be removed for females
  y.probes <- anno %>% dplyr::filter(chr == "chrY") %$% Name

  ## Check number of males/females
  sex_male <- sex %>% dplyr::filter(predictedSex == "M")
  sex_female <- sex %>% dplyr::filter(predictedSex == "F")

  ## Get beadnr and detp metrics
  # Note: for y-probes only males (based on predicted Sex) are used to calculate probe metrics
  if(nrow(sex_female) == 0) {
    # Detp sample metrics
    detp_percentage_samples_total <- Matrix::colMeans(mtrx == 2 | mtrx == 3)
    detp_percentage_samples_total <- tibble::tibble(Sample_Name = names(detp_percentage_samples_total), detP_Percentage = detp_percentage_samples_total)

    # Beadnr sample metrics
    beadnr_percentage_total <- Matrix::colMeans(mtrx == 1 | mtrx == 3)
    beadnr_percentage_samples_total <- tibble::tibble(Sample_Name = names(beadnr_percentage_total), BeadNr_Percentage = beadnr_percentage_total)

    if(probe_qc_method == "post_sampleqc") {
      sample_metrics <- .join_sample_metrics(methylaid_metrics = methylaid_metrics, bscon = bscon, sex = sex,
                                             rg = rg, detp = detp_percentage_samples_total,
                                             beadnr = beadnr_percentage_samples_total,  ibs = ibs)
      sample_outliers <- .get_sample_outliers(sample_metrics, thresholds, samplesheet = samplesheet)

      # Number of outliers
      outliers <- if(nrow(sample_outliers) == 0) c() else sample_outliers$Sample_Name
      outliers_percentage <- if(nrow(sample_outliers) == 0) 0 else (length(outliers)*100)/nrow(sample_outliers)

      ## Update matrix
      mtrx <- mtrx[,!(colnames(mtrx) %in% outliers), drop = FALSE]
      sex_male <- sex_male %>% dplyr::filter(!(Sample_Name %in% outliers))
      sex_female <- sex_female %>% dplyr::filter(!(Sample_Name %in% outliers))
      probe_qc_method_ <- "post_sampleqc"


    } else {
      probe_qc_method_ <- "pre_sampleqc"
    }
    # Probe metrics - beadnr
    beadnr_percentage_probes_total <- Matrix::rowMeans(mtrx == 1 | mtrx == 3)
    beadnr_percentage_probes_total <- tibble::tibble(Probe = names(beadnr_percentage_probes_total), BeadNr_Percentage = beadnr_percentage_probes_total)

    # Probe Metrics - detp
    detp_percentage_probes <- Matrix::rowMeans(mtrx == 2 | mtrx == 3)
    detp_percentage_probes_total <- tibble::tibble(Probe = names(detp_percentage_probes), detP_Percentage = detp_percentage_probes)

  } else if(nrow(sex_male) == 0) {
    # Sample Metrics - detp
    mtrx_no_y <- mtrx[!(rownames(mtrx) %in% y.probes), ]
    detp_percentage_samples <- Matrix::colMeans(mtrx_no_y == 2 | mtrx_no_y == 3)
    detp_percentage_samples_total <- tibble::tibble(Sample_Name = names(detp_percentage_samples),
                                                    detP_Percentage = detp_percentage_samples)

    # Sample metrics - beadnr
    beadnr_percentage_samples_total <- Matrix::colMeans(mtrx_no_y == 1 | mtrx_no_y == 3)
    beadnr_percentage_samples_total <- tibble::tibble(Sample_Name = names(beadnr_percentage_samples_total),
                                                      BeadNr_Percentage = beadnr_percentage_samples_total)

    if(probe_qc_method == "post_sampleqc") {
      sample_metrics <- .join_sample_metrics(methylaid_metrics = methylaid_metrics, bscon = bscon, sex = sex,
                                             rg = rg, detp = detp_percentage_samples_total,
                                             beadnr = beadnr_percentage_samples_total, ibs = ibs)
      sample_outliers <- .get_sample_outliers(sample_metrics, thresholds, samplesheet=samplesheet)

      # Number of outliers
      outliers <- if(nrow(sample_outliers) == 0) c() else sample_outliers$Sample_Name
      outliers_percentage <- if(nrow(sample_outliers) == 0) 0 else (length(outliers)*100)/nrow(sample_outliers)

      ## Remove outliers
      ## Update matrix
      mtrx <- mtrx[,!(colnames(mtrx) %in% outliers), drop = FALSE]
      sex_male <- sex_male %>% dplyr::filter(!(Sample_Name %in% outliers))
      sex_female <- sex_female %>% dplyr::filter(!(Sample_Name %in% outliers))
      probe_qc_method_ <- "post_sampleqc"

    } else {
      probe_qc_method_ <- "pre_sampleqc"
    }

    # Probe metrics - beadnr
    beadnr_percentage_probes_total <- Matrix::rowMeans(mtrx == 1 | mtrx == 3)
    beadnr_percentage_probes_total <- tibble::tibble(Probe = names(beadnr_percentage_probes_total),
                                                     BeadNr_Percentage = beadnr_percentage_probes_total)

    # Probe Metrics - detp
    detp_percentage_probes <- Matrix::rowMeans(mtrx == 2 | mtrx == 3)
    detp_percentage_probes_total <- tibble::tibble(Probe = names(detp_percentage_probes),
                                                   detP_Percentage = detp_percentage_probes)
  } else {
    # Sample metrics
    mtrx_no_y <- mtrx[!(rownames(mtrx) %in% y.probes),, drop = FALSE]
    mtrx_male <- mtrx[,as.character(sex_male$Sample_Name), drop = FALSE]

    # Detp
    detp_percentage_samples_male <- Matrix::colMeans(mtrx_male == 2 | mtrx_male == 3)
    mtrx_no_y_female <- mtrx_no_y[,as.character(sex_female$Sample_Name),drop=FALSE]
    detp_percentage_samples_female <- Matrix::colMeans(mtrx_no_y_female == 2 | mtrx_no_y_female == 3)
    detp_percentage_samples_total <- c(detp_percentage_samples_male, detp_percentage_samples_female)
    detp_percentage_samples_total <- tibble::tibble(Sample_Name = names(detp_percentage_samples_total),
                                                    detP_Percentage = detp_percentage_samples_total)

    # Sample metrics - beadnr
    beadnr_percentage_samples_male <- Matrix::colMeans(mtrx_male == 1 | mtrx_male == 3)
    beadnr_no_y_female <- mtrx_no_y[,as.character(sex_female$Sample_Name), drop = FALSE]
    beadnr_percentage_samples_female <- Matrix::colMeans(mtrx_no_y_female == 1 | mtrx_no_y_female == 3)
    beadnr_percentage_samples_total <- c(beadnr_percentage_samples_male, beadnr_percentage_samples_female)
    beadnr_percentage_samples_total <- tibble::tibble(Sample_Name = names(beadnr_percentage_samples_total),
                                                      BeadNr_Percentage = beadnr_percentage_samples_total)

    if(probe_qc_method == "post_sampleqc") {
      sample_metrics <- .join_sample_metrics(methylaid_metrics = methylaid_metrics, bscon = bscon, sex = sex,
                                             rg = rg, detp = detp_percentage_samples_total,
                                             beadnr = beadnr_percentage_samples_total, ibs = ibs)
      sample_outliers <- .get_sample_outliers(sample_metrics, thresholds, samplesheet = samplesheet)

      # Number of outliers
      outliers <- if(nrow(sample_outliers) == 0) c() else sample_outliers$Sample_Name
      outliers_percentage <- if(nrow(sample_outliers) == 0) 0 else (length(outliers)*100)/nrow(sample_outliers)

      ## Remove outliers
      ## Update detp and beadnr
      mtrx <- mtrx[,!(colnames(mtrx) %in% outliers), drop = FALSE]
      sex_male <- sex_male %>% dplyr::filter(!(Sample_Name %in% outliers))
      sex_female <- sex_female %>% dplyr::filter(!(Sample_Name %in% outliers))
      probe_qc_method_ <- "post_sampleqc"

    } else {
      probe_qc_method_ <- "pre_sampleqc"
    }

    # Note that removing sample qc failures can cause the dataset to become all-female or all-male!
    if(nrow(sex_male) == 0 || nrow(sex_female) == 0) {
      # Probe Metrics - detp
      detp_percentage_probes <- Matrix::rowMeans(mtrx == 2 | mtrx == 3)
      detp_percentage_probes_total <- tibble::tibble(Probe = names(detp_percentage_probes),
                                                     detP_Percentage = detp_percentage_probes)

      # Probe metrics - beadnr
      beadnr_percentage_probes_total <- Matrix::rowMeans(mtrx == 1 | mtrx == 3)
      beadnr_percentage_probes_total <- tibble::tibble(Probe = names(beadnr_percentage_probes_total),
                                                       BeadNr_Percentage = beadnr_percentage_probes_total)
    } else {
      # Probe metrics - detp
      mtrx_male_y <- mtrx[(rownames(mtrx) %in% y.probes),as.character(sex_male$Sample_Name), drop = FALSE]
      mtrx_no_y <- mtrx[!(rownames(mtrx) %in% y.probes), ]
      detp_percentage_y_male <- Matrix::rowMeans(mtrx_male_y == 2 | mtrx_male_y == 3)
      detp_percentage_probes <- Matrix::rowMeans(mtrx_no_y == 2 | mtrx_no_y == 3)
      detp_percentage_probes_total <- c(detp_percentage_y_male, detp_percentage_probes)
      detp_percentage_probes_total <- tibble::tibble(Probe = names(detp_percentage_probes_total),
                                                     detP_Percentage = detp_percentage_probes_total)

      # Probe metrics - beadnr
      beadnr_percentage_y_male <- Matrix::rowMeans(mtrx_male_y == 1 | mtrx_male_y == 3)
      beadnr_percentage_probes <- Matrix::rowMeans(mtrx_no_y == 1 | mtrx_no_y == 3)
      beadnr_percentage_probes_total <- c(beadnr_percentage_y_male, beadnr_percentage_probes)
      beadnr_percentage_probes_total <- tibble::tibble(Probe = names(beadnr_percentage_probes_total),
                                                       BeadNr_Percentage = beadnr_percentage_probes_total)
    }
  }
  if(probe_qc_method_ == "post_sampleqc" && nrow(sample_outliers) == nrow(samplesheet)) {
    beadnr_percentage_probes_total <- tibble::tibble(Probe = rownames(mtrx),
                                                     BeadNr_Percentage = 1)
    detp_percentage_probes_total <- tibble::tibble(Probe = rownames(mtrx),
                                                   detP_Percentage = 1)
  }
  detp_percentage_samples_total$detP_Percentage <- unname(detp_percentage_samples_total$detP_Percentage)
  beadnr_percentage_probes_total$BeadNr_Percentage <- unname(beadnr_percentage_probes_total$BeadNr_Percentage)
  list(detp_samples = detp_percentage_samples_total, beadnr_samples = beadnr_percentage_samples_total,
       detp_probes = detp_percentage_probes_total, beadnr_probes = beadnr_percentage_probes_total,
       probe_qc_method = probe_qc_method_)
}


.get_cellcounts <- function(beta, epidish_method, tissue) {
  if(tissue == "blood") {
    ref <- EpiDISH::centDHSbloodDMC.m[,c("B", "NK", "CD4T", "CD8T", "Mono", "Neutro", "Eosino")]
  } else if (tissue == "breast") {
    ref <- EpiDISH::centEpiFibFatIC.m
  } else if (tissue == "epithelium") {
    ref <- EpiDISH::centEpiFibIC.m
  }
  cellcounts <- EpiDISH::epidish(beta.m = beta, ref.m = ref, method = epidish_method)
  cellcounts <- cellcounts$estF
  samples <- rownames(cellcounts)
  cellcounts <- dplyr::as_tibble(cellcounts) %>%
    dplyr::mutate(Sample_Name = samples) %>%
    dplyr::select(Sample_Name, dplyr::everything())
  cellcounts
}

.get_predicted_age <- function(beta) {
  horvath_age <- .get_horvath_age(beta)
  zhang_age <- .get_zhang_age(beta)
  combined <- horvath_age %>%
    dplyr::left_join(zhang_age, by = "Sample_Name")
  combined
}

.get_horvath_age <- function(beta) {
  beta[is.na(beta)] <- 0
  intercept <- coefs[["age.horvath"]]["(Intercept)"]
  coef <- coefs[["age.horvath"]][names(coefs[["age.horvath"]]) != "(Intercept)"]
  coef <- coef[names(coef) %in% rownames(beta)]
  age <- intercept + as.vector(coef %*% beta[names(coef),,drop=FALSE])
  age <- ifelse(age < 0, (1 + 20)*exp(age) - 1, (1 + 20)*age + 20)
  age <- tibble::tibble(
    Sample_Name = colnames(beta),
    Predicted_Age =  age
  )
  age
}

# Adapted from:
# https://github.com/qzhang314/DNAm-based-age-predictor/blob/master/pred.R
.get_zhang_age <- function(beta) {

  # Scale
  beta_norm <- apply(beta, 2, scale)
  rownames(beta_norm) <- rownames(beta)
  # Get overlapping probes
  probes_intersect_en <- intersect(rownames(beta_norm), coefs[["age.en.coef"]]$Probe)
  probes_intersect_blup <- intersect(rownames(beta_norm), coefs[["age.blup.coef"]]$Probe)
  en_int <- coefs[["age.en.coef"]][1,2]
  blup_int <- coefs[["age.blup.coef"]][1,2]
  en_coef <- coefs[["age.en.coef"]] %>% dplyr::filter(Probe %in% probes_intersect_en)
  blup_coef <- coefs[["age.blup.coef"]] %>% dplyr::filter(Probe %in% probes_intersect_blup)

  # Set missings to zero
  beta_norm[is.na(beta_norm)] <- 0

  # Predict
  age <- tibble::tibble(
    Sample_Name = colnames(beta_norm),
    Predicted_Age_Zhang_en =  as.vector(en_coef[["coef"]] %*% beta_norm[en_coef[["Probe"]],,drop=FALSE]) + as.numeric(en_int),
    Predicted_Age_Zhang_blup = as.vector(blup_coef[["coef"]]  %*% beta_norm[blup_coef[["Probe"]],,drop=FALSE]) + as.numeric(blup_int)
  )
  age
}

# .get_epismoker <- function(beta, samplesheet) {
#   # Recode sex column to 2,1,0
#   # Convert samplesheet to data.frame with rownames (the epismoker functions expects te Sample_Name column to be rownames)
#   samplesheet <- samplesheet %>%
#     dplyr::mutate(Sex = dplyr::case_when(
#       Sex == "F" ~ 2,
#       Sex == "M" ~ 1,
#       TRUE ~ 0
#     ))  %>%
#     as.data.frame()
#   rownames(samplesheet) <- samplesheet$Sample_Name
#
#   # Replace missing beta-values with zero
#   beta[is.na(beta)] <- 0
#
#   # Run epismoker
#   epismoker <- suppressMessages(EpiSmokEr::epismoker(beta, samplesheet = samplesheet, method = "all",
#                                                      ref.Elliott = EpiSmokEr::Illig_data,
#                                                      ref.Zhang = EpiSmokEr::Zhangetal_cpgs,
#                                                      ref.CS = EpiSmokEr::CS_final_coefs,
#                                                      ref.FS = EpiSmokEr::FS_final_coefs,
#                                                      ref.NS = EpiSmokEr::NS_final_coefs))
#
#   # Reformat and change factors to character columns
#   epismoker <- epismoker %>%
#     dplyr::rename(Sample_Name = SampleName) %>%
#     purrr::map_if(is.factor, as.character) %>%
#     dplyr::as_tibble()
#
#   epismoker
# }

.get_pms <- function(beta, verbose, logfile) {
  alcohol <- .get_alcohol(beta, verbose=verbose, logfile=logfile)
  smoking <- .get_smoking(beta, verbose=verbose, logfile=logfile)
  predicted_age <- .get_predicted_age(beta)
  BMI <- .get_BMI(beta, verbose=verbose, logfile=logfile)
  pms <- predicted_age %>%
    dplyr::left_join(smoking, by = "Sample_Name") %>%
    dplyr::left_join(alcohol, by = "Sample_Name") %>%
    dplyr::left_join(BMI, by = "Sample_Name")
  pms
}

.get_alcohol <- function(beta, verbose, logfile) {
  coef_ <- coefs[["alcohol"]] %>%
    dplyr::filter(Probe %in% rownames(beta))

  if(nrow(coef_) > 0) {
    beta_ <- beta[coef_[["Probe"]],,drop=FALSE]
    beta_[is.na(beta_)] <- 0
    if(verbose) message(sprintf("Using %s/%s probes for alcohol prediction", nrow(coef_), length(coefs[["alcohol"]][["Probe"]])))
    if(!is.null(logfile)) {cat(sprintf("Using %s/%s probes for alcohol prediction\n", nrow(coef_), length(coefs[["alcohol"]][["Probe"]])),
                               file = logfile, append = TRUE)}

    alcohol <- tibble::tibble(Sample_Name = colnames(beta_),
                      Alcohol_PMS = (t(beta_[coef_[["Probe"]],]) %*% coef_[["coef"]])[,1]
    )
    alcohol
  } else {
    if(verbose) message("No probes left for alcohol prediction")
    if(!is.null(logfile)) {cat("No probes left for alcohol prediction\n",
                               file = logfile, append = TRUE)}
  }
}

.get_BMI <- function(beta, verbose, logfile) {

  intercept <- coefs[["BMI"]] %>% dplyr::filter(Probe == "(Intercept)") %$% coef
  coef_ <- coefs[["BMI"]] %>%
    dplyr::filter(Probe %in% rownames(beta))

  if(nrow(coef_) > 0) {
    beta_ <- beta[coef_[["Probe"]],,drop=FALSE]
    beta_[is.na(beta_)] <- 0
    if(verbose) message(sprintf("Using %s/%s probes for BMI prediction", nrow(coef_), length(coefs[["BMI"]][["Probe"]])-1))
    if(!is.null(logfile)) {cat(sprintf("Using %s/%s probes for BMI prediction\n", nrow(coef_), length(coefs[["BMI"]][["Probe"]])-1),
                               file = logfile, append = TRUE)}

    BMI <- tibble::tibble(Sample_Name = colnames(beta_),
                  BMI_PMS = intercept + (t(beta_[coef_[["Probe"]],]) %*% coef_[["coef"]])[,1]
    )
    BMI

  } else {
    if(verbose) message("No probes left for BMI prediction")
    if(!is.null(logfile)) {cat("No probes left for BMI prediction\n",
                               file = logfile, append = TRUE)}
  }
}

.get_smoking <- function(beta, verbose, logfile) {
  coef_ <- coefs[["smoking"]] %>%
    dplyr::filter(Probe %in% rownames(beta))
  refdata <- coefs[["smoking_ref"]] %>%
    dplyr::filter(cpgs %in% coef_$Probe)

  if(nrow(coef_) > 0) {
    beta_ <- beta[coef_[["Probe"]],,drop=FALSE]
    beta_[is.na(beta_)] <- 0
    if(verbose) message(sprintf("Using %s/%s probes for smoking prediction", nrow(coef_), length(coefs[["smoking"]][["Probe"]])))
    if(!is.null(logfile)) {cat(sprintf("Using %s/%s probes for smoking prediction\n", nrow(coef_), length(coefs[["smoking"]][["Probe"]])),
                               file = logfile, append = TRUE)}

    beta_ <- rbind(
      beta_[(refdata %>% dplyr::filter(all_effect >= 0) %$% cpgs),,drop=FALSE ] - (refdata %>% dplyr::filter(all_effect >= 0) %$% reference_never_median_beta_all),
      matrix(rep(refdata %>% dplyr::filter(all_effect < 0) %$% reference_never_median_beta_all, ncol(beta_)), nrow=nrow(refdata %>% dplyr::filter(all_effect < 0))) - beta_[(refdata %>% dplyr::filter(all_effect < 0) %$% cpgs),,drop=FALSE ]
    )
    smoking <- tibble::tibble(Sample_Name = colnames(beta_),
                      Smoking_PMS = (t(beta_[coef_[["Probe"]],]) %*% coef_[["coef"]])[,1]
    )
    smoking

  } else {
    if(verbose) message("No probes left for smoking prediction")
    if(!is.null(logfile)) {cat("No probes left for smoking prediction\n",
                               file = logfile, append = TRUE)}
  }
}

.get_heatmap_pvals <- function(samplesheet, PCs, cellcounts, pms, variables, factor, logfile) {

  ## Subset samplesheet
  samplesheet_ <- samplesheet %>% dplyr::filter(Sample_Name %in% PCs$Sample_Name)
  PCs_ <- PCs[match(samplesheet$Sample_Name, PCs$Sample_Name), ] # Match

  ## Select chosen number of PCs
  if(ncol(PCs) > 11) {
    PCs_ <- PCs_[,2:11]
  } else {
    PCs_ <- PCs_[,2:ncol(PCs)]
  }

  # Default = cellcounts + Sex
  default_variables <- c("Sex", "Predicted_Age", "Smoking_PMS", "Alcohol_PMS", "BMI_PMS", setdiff(colnames(cellcounts), "Sample_Name"))
  default_fct <- c("Sex")

  # Combine default and user-specified variables
  if(!is.null(variables)) {
    variables <- stringr::str_trim(unlist(stringr::str_split(variables, pattern = ",")))
    variables_total <- unique(c(default_variables, variables))
    if(!is.null(factor)) {
      fct_variables <- stringr::str_trim(unlist(stringr::str_split(factor, pattern = ",")))
      fct_total <- unique(c(default_fct, fct_variables))
    } else {
      fct_total <- default_fct
    }
  } else {
    variables_total <- default_variables
    fct_total <- default_fct
  }

  # Join samplesheet with cellcounts
  samplesheet_ <- samplesheet_[,!colnames(samplesheet_) %in% setdiff(default_variables, "Sex"),drop=FALSE] %>%
    dplyr::left_join(cellcounts, by = "Sample_Name") %>%
    dplyr::left_join(pms, by = "Sample_Name")

  ## Make factors
  samplesheet_ <- samplesheet_[variables_total]
  samplesheet_[fct_total] <- lapply(samplesheet_[fct_total], base::factor)
  samplesheet_[!colnames(samplesheet_) %in% fct_total] <- lapply(samplesheet_[!colnames(samplesheet_) %in% fct_total], function(x) as.numeric(as.character(x)))

  # Check if any of the factors have only one level
  nr.levels <- purrr::map_dbl(samplesheet_[fct_total], nlevels)
  nr.nonmissing <- purrr::map_dbl(samplesheet_[fct_total], function(x) sum(!is.na(x)))

  # Remove variables with only level (will give an error later on otherwise)
  if(any(nr.levels == 1)) {
    nr.levels <- nr.levels[nr.levels == 1]
    if(!is.null(logfile)) cat(sprintf("\nThe following variables have only one level, and will therefore be excluded from the association tests:\n"), file = logfile, append = TRUE)
    message(sprintf("The following variables have only one level, and will therefore be excluded from the association tests:"))
    for(name in names(nr.levels)) {
        if(!is.null(logfile)) cat(sprintf("%s\n", name), file = logfile, append = TRUE)
        message(sprintf("\t%s", name))
    }
    samplesheet_ <- samplesheet_[,!(colnames(samplesheet_) %in% names(nr.levels))]
  }

  # Remove quantitative variables with no variance
  vars <- purrr::map_dbl(samplesheet_[!colnames(samplesheet_) %in% fct_total], var)
  if(any(vars == 0)) {
   vars <- vars[vars == 0]
   if(!is.null(logfile)) {cat(sprintf("\nThe following variables have zero variance, and will therefore be excluded from the association tests:\n"),
        file = logfile, append = TRUE)}
    message(sprintf("The following variables have zero variance, and will therefore be excluded from the association tests:"))
    for(name in names(vars)) {
      if(!is.null(logfile)) cat(sprintf("%s\n", name), file = logfile, append = TRUE)
      message(sprintf("\t%s", name))
    }
    samplesheet_ <- samplesheet_[,!(colnames(samplesheet_) %in% names(vars))]
  }


  # Remove variables with too many missing values
  if(any(nr.nonmissing < 2)) {
    nr.nonmissing <- nr.nonmissing[nr.nonmissing < 2]
    if(!is.null(logfile)) cat(sprintf("The following variables too many missing values, and will therefore be excluded from the association tests:\n"), file = logfile, append = TRUE)
    message(sprintf("The following variables too many missing values, and will therefore be excluded from the association tests:"))
    for(name in names(nr.nonmissing)) {
        if(!is.null(logfile)) cat(sprintf("%s\n", name), file = logfile, append = TRUE)
        message(sprintf("\t%s", name))
    }
    samplesheet_ <- samplesheet_[,!(colnames(samplesheet_) %in% names(nr.nonmissing))]
  }

  ## Calculate P-values and plot
  PCs_ <- dplyr::as_tibble(PCs_)
  pvals <- .do.tests(PCs_, samplesheet_, twolevel_method = "lm")
  pvals_adjusted <- .padjust(pvals)
  pvals_adjusted_log <- as.matrix(-log10(pvals_adjusted))
  pvals_adjusted_log <- reshape2::melt(pvals_adjusted_log, varnames = c("Variable", "PC"))

  # If any -log10p is 'Inf' set it to the maximum pvalue
  # This occurs when p-value is exactly 0
  mx <- max(pvals_adjusted_log$value[!is.infinite(pvals_adjusted_log$value)])
  pvals_adjusted_log <- pvals_adjusted_log %>% dplyr::mutate(value = ifelse(is.infinite(value), mx, value))
  pvals_adjusted_log
}

.padjust <- function(pvals) {
   nr.tests <- nrow(pvals) * ncol(pvals)
   adjusted_pvals <- apply(pvals, 2, FUN = function(pval) p.adjust(pval, n = nr.tests, method = 'bonferroni'))
   adjusted_pvals
 }

 ## Association between a PC and a variable
 ## Factorize  categorical variables manually beforehand!
.do.test <- function(PC, covar, twolevel_method = 'wilcox') {
   if (nlevels(covar) > 2) {
     pval <- kruskal.test(PC ~ covar)$p.value
   } else if (nlevels(covar) == 2 & twolevel_method == 'wilcox') {
     pval <- wilcox.test(PC ~ covar, alternative = "two.sided")$p.value
   } else if (nlevels(covar) == 2 & twolevel_method == 'lm') {
     pval <- summary(lm(PC ~ covar))$coef[2,4]
   }
   else {
     pval <- summary(lm(PC ~ covar))$coef[2,4]
   }
   pval
 }

 ## Perform above associations tests for all PCs and all variables
.do.tests <- function(PCs, covars, twolevel_method = 'wilcox') {
   tests <- purrr::map(covars, function(x) purrr::map_dbl(PCs, .do.test, covar = x, twolevel_method = twolevel_method))
   tests <- as.matrix(do.call(dplyr::bind_rows, tests))
   rownames(tests) <- colnames(covars)
   tests
}


## Minfi internal code (https://github.com/hansenlab/minfi/blob/master/R/utils.R)

## getSex from:
.getsex <- function(CN = NULL, xIndex = NULL, yIndex = NULL, cutoff = -2) {
  if (is.null(CN) || is.null(xIndex) || is.null(yIndex)) {
    stop("must provide CN, xIndex, and yIndex")
  }
  # TODO: This does not handle only females or only males. This ought to be
  #       handled by the 'centers' (see below) being too close together
  xMed <- matrixStats::colMedians(CN, rows = xIndex, na.rm = TRUE)
  yMed <- matrixStats::colMedians(CN, rows = yIndex, na.rm = TRUE)
  dd <- yMed - xMed
 # k <- kmeans(dd, centers = c(min(dd), max(dd)))

  sex0 <- ifelse(dd < cutoff, "F", "M")
 # sex0 <- .checkSex(sex0)
#  sex1 <- ifelse(k$cluster == which.min(k$centers), "F", "M")
 # sex1 <- .checkSex(sex1)

  # if (!identical(sex0, sex1)) {
  #   warning("An inconsistency was encountered while determining sex. One ",
  #           "possibility is that only one sex is present. We recommend ",
  #           "further checks, for example with the plotSex function.")
  # }
  df <-tibble::tibble(Sample_Name = colnames(CN), xMed = xMed, yMed = yMed, predictedSex = sex0)
  #rownames(df) <- colnames(CN)
  df
}

.getSex <- function(object = NULL, cutoff = -2){
  .isGenomicOrStop(object)
  if (is(object, "GenomicMethylSet")) CN <- minfi::getCN(object)
  if (is(object, "GenomicRatioSet")) CN <- minfi::getCN(object)
  # TODO: Add test for logarithmic scale or non-log scale
  xIndex <- which(as.character(GenomicRanges::seqnames(object)) == "chrX")
  yIndex <- which(as.character(GenomicRanges::seqnames(object)) == "chrY")
  .getsex(
    CN = CN,
    xIndex = xIndex,
    yIndex = yIndex,
    cutoff = cutoff)
}

## fix in detectionP (from: https://github.com/hansenlab/minfi/blob/master/R/detectionP.R)

.detectionP_fix <- function(rgSet, type = "m+u") {
  .isRGOrStop(rgSet)

  # Extract data to pass to low-level function that constructs `detP`
  locusNames <- minfi::getManifestInfo(rgSet, "locusNames")
  controlIdx <- minfi::getControlAddress(rgSet, controlType = "NEGATIVE")
  Red <- minfi::getRed(rgSet)
  Green <- minfi::getGreen(rgSet)
  TypeI.Red <- minfi::getProbeInfo(rgSet, type = "I-Red")
  TypeI.Green <- minfi::getProbeInfo(rgSet, type = "I-Green")
  TypeII <- minfi::getProbeInfo(rgSet, type = "II")

  # Construct `detP`
  detP <- .detectionP_fix_internal(
    Red = Red,
    Green = Green,
    locusNames = locusNames,
    controlIdx = controlIdx,
    TypeI.Red = TypeI.Red,
    TypeI.Green = TypeI.Green,
    TypeII = TypeII)

  detP
}

.detectionP_fix_internal <- function(Red, Green, locusNames, controlIdx, TypeI.Red, TypeI.Green, TypeII,
         ...) {
  # Set up output matrix with appropriate dimensions and type
  detP <- matrix(NA_real_,
                 nrow = length(locusNames),
                 ncol = ncol(Red),
                 dimnames = list(locusNames, colnames(Red)))

  # Compute summary statistics needed for calculations
  rBg <- Red[controlIdx, , drop = FALSE]
  rMu <- matrixStats::colMedians(rBg)
  rSd <- matrixStats::colMads(rBg)
  gBg <- Green[controlIdx, , drop = FALSE]
  gMu <- matrixStats::colMedians(gBg)
  gSd <- matrixStats::colMads(gBg)

  # Fill output matrix
  for (j in seq_len(ncol(detP))) {
    # Type I Red
    intensity <- Red[TypeI.Red$AddressA, j] + Red[TypeI.Red$AddressB, j]
    detP[TypeI.Red$Name, j] <- pnorm(
      q = intensity,
      mean = 2 * rMu[j],
      sd = 2 * rSd[j],
      lower.tail = FALSE)
    # Type I Green
    intensity <- Green[TypeI.Green$AddressA, j] +
      Green[TypeI.Green$AddressB, j]
    detP[TypeI.Green$Name, j] <- pnorm(
      q = intensity,
      mean = 2 * gMu[j],
      sd = 2 * gSd[j],
      lower.tail = FALSE)
    # Type II
    intensity <- Red[TypeII$AddressA, j] + Green[TypeII$AddressA, j]
    detP[TypeII$Name, j] <- pnorm(
      q = intensity,
      mean = rMu[j] + gMu[j],
      sd = rSd[j] + gSd[j],
      lower.tail = FALSE)
  }

  # Return output matrix
  detP
}


## fix in detectionP (from: https://github.com/hansenlab/minfi/blob/master/R/utils.R)
.isGenomicOrStop <- function(object) {
  if (!is(object, "GenomicMethylSet") && !is(object, "GenomicRatioSet")) {
    stop("object is of class '", class(object), "', but needs to be of ",
         "class 'GenomicMethylSet' or 'GenomicRatioSet'")
  }
}

.isRGOrStop <- function(object) {
  if (!is(object, "RGChannelSet")) {
    stop("object is of class '", class(object), "', but needs to be of ",
         "class 'RGChannelSet' or 'RGChannelSetExtended'")
  }
}

## Extraction of the Control matrix; from: https://github.com/hansenlab/minfi/blob/master/R/preprocessFunnorm.R
.buildControlMatrix450k <- function(extractedData) {
  getCtrlsAddr <- function(exType, index) {
    ctrls <- ctrlsList[[index]]
    addr <- ctrls$Address
    names(addr) <- ctrls$ExtendedType
    na.omit(addr[exType])
  }

  array <- extractedData$array
  greenControls <- extractedData$greenControls
  redControls <- extractedData$redControls
  controlNames <- names(greenControls)
  ctrlsList <- extractedData$ctrlsList

  ## Bisulfite conversion extraction for probe type II:
  index <- match("BISULFITE CONVERSION II", controlNames)
  redControls.current <- redControls[[ index ]]
  bisulfite2 <- matrixStats::colMeans2(redControls.current, na.rm = TRUE)

  ## Bisulfite conversion extraction for probe type I:
  index <- match("BISULFITE CONVERSION I", controlNames)
  if (array=="IlluminaHumanMethylation450k"){
    addr <- getCtrlsAddr(exType = sprintf("BS Conversion I%sC%s", c(" ", "-", "-"), 1:3), index = index)
  } else {
    addr <- getCtrlsAddr(exType = sprintf("BS Conversion I%sC%s", c("-", "-"), 1:2), index = index)
  }
  greenControls.current <- greenControls[[ index ]][addr,,drop=FALSE]
  if (array=="IlluminaHumanMethylation450k"){
    addr <- getCtrlsAddr(exType = sprintf("BS Conversion I-C%s", 4:6), index = index)
  } else {
    addr <- getCtrlsAddr(exType = sprintf("BS Conversion I-C%s", 3:5), index = index)
  }
  redControls.current <- redControls[[ index ]][addr,, drop=FALSE]
  if (nrow(redControls.current)==nrow(greenControls.current)){
    bisulfite1 <- matrixStats::colMeans2(redControls.current + greenControls.current, na.rm = TRUE)
  } else {
    bisulfite1 <- matrixStats::colMeans2(redControls.current, na.rm=TRUE) + matrixStats::colMeans2(greenControls.current, na.rm = TRUE)
  }


  ## Staining
  index <- match("STAINING", controlNames)
  addr <- getCtrlsAddr(exType = "Biotin (High)", index = index)
  stain.green <- t(greenControls[[ index ]][addr,,drop=FALSE])
  addr <- getCtrlsAddr(exType = "DNP (High)", index = index)
  stain.red <- t(redControls[[ index ]][addr,, drop=FALSE ])

  ## Extension
  index <-    match("EXTENSION", controlNames)
  addr <- getCtrlsAddr(exType = sprintf("Extension (%s)", c("A", "T")), index = index)
  extension.red <- t(redControls[[index]][addr,,drop=FALSE])
  colnames(extension.red) <- paste0("extRed", 1:ncol(extension.red))
  addr <- getCtrlsAddr(exType = sprintf("Extension (%s)", c("C", "G")), index = index)
  extension.green <- t(greenControls[[index]][addr,,drop=FALSE])
  colnames(extension.green) <- paste0("extGrn", 1:ncol(extension.green))

  ## Hybridization should be monitored only in the green channel
  index <- match("HYBRIDIZATION", controlNames)
  hybe <- t(greenControls[[index]])
  colnames(hybe) <- paste0("hybe", 1:ncol(hybe))

  ## Target removal should be low compared to hybridization probes
  index <- match("TARGET REMOVAL", controlNames)
  targetrem <- t(greenControls[[index]])
  colnames(targetrem) <- paste0("targetrem", 1:ncol(targetrem))

  ## Non-polymorphic probes
  index <- match("NON-POLYMORPHIC", controlNames)
  addr <- getCtrlsAddr(exType = sprintf("NP (%s)", c("A", "T")), index = index)
  nonpoly.red <- t(redControls[[index]][addr, ,drop=FALSE])
  colnames(nonpoly.red) <- paste0("nonpolyRed", 1:ncol(nonpoly.red))
  addr <- getCtrlsAddr(exType = sprintf("NP (%s)", c("C", "G")), index = index)
  nonpoly.green <- t(greenControls[[index]][addr, ,drop=FALSE])
  colnames(nonpoly.green) <- paste0("nonpolyGrn", 1:ncol(nonpoly.green))

  ## Specificity II
  index <- match("SPECIFICITY II", controlNames)
  greenControls.current <- greenControls[[index]]
  redControls.current <- redControls[[index]]
  spec2.green <- t(greenControls.current)
  colnames(spec2.green) <- paste0("spec2Grn", 1:ncol(spec2.green))
  spec2.red <- t(redControls.current)
  colnames(spec2.red) <- paste0("spec2Red", 1:ncol(spec2.red))
  spec2.ratio <- matrixStats::colMeans2(greenControls.current, na.rm = TRUE) /
    matrixStats::colMeans2(redControls.current, na.rm = TRUE)

  ## Specificity I
  index <- match("SPECIFICITY I", controlNames)
  addr <- getCtrlsAddr(exType = sprintf("GT Mismatch %s (PM)", 1:3), index = index)
  greenControls.current <- greenControls[[index]][addr,,drop=FALSE]
  redControls.current <- redControls[[index]][addr,,drop=FALSE]
  spec1.green <- t(greenControls.current)
  colnames(spec1.green) <- paste0("spec1Grn", 1:ncol(spec1.green))
  spec1.ratio1 <- matrixStats::colMeans2(redControls.current, na.rm = TRUE) /
    matrixStats::colMeans2(greenControls.current, na.rm = TRUE)

  index <- match("SPECIFICITY I", controlNames) # Added that line
  addr <- getCtrlsAddr(exType = sprintf("GT Mismatch %s (PM)", 4:6), index = index)
  greenControls.current <- greenControls[[index]][addr,,drop=FALSE]
  redControls.current <- redControls[[index]][addr,,drop=FALSE]
  spec1.red <- t(redControls.current)
  colnames(spec1.red) <- paste0("spec1Red", 1:ncol(spec1.red))
  spec1.ratio2 <- matrixStats::colMeans2(greenControls.current, na.rm = TRUE) /
    matrixStats::colMeans2(redControls.current, na.rm = TRUE)
  spec1.ratio <- (spec1.ratio1 + spec1.ratio2) / 2

  ## Normalization probes:
  index <- match(c("NORM_A"), controlNames)
  normA <- matrixStats::colMeans2(redControls[[index]], na.rm = TRUE)
  index <- match(c("NORM_T"), controlNames)
  normT <- matrixStats::colMeans2(redControls[[index]], na.rm = TRUE)
  index <- match(c("NORM_C"), controlNames)
  normC <- matrixStats::colMeans2(greenControls[[index]], na.rm = TRUE)
  index <- match(c("NORM_G"), controlNames)
  normG <- matrixStats::colMeans2(greenControls[[index]], na.rm = TRUE)

  dyebias <- (normC + normG)/(normA + normT)

  oobG <- extractedData$oob$greenOOB
  oobR <- extractedData$oob$redOOB
  oob.ratio <- oobG[2,]/oobR[2,]
  oobG <- t(oobG)
  colnames(oobG) <- paste0("oob", c(1,50,99))

  model.matrix <- cbind(
    bisulfite1, bisulfite2, extension.green, extension.red, hybe,
    stain.green, stain.red, nonpoly.green, nonpoly.red,
    targetrem, spec1.green, spec1.red, spec2.green, spec2.red, spec1.ratio1,
    spec1.ratio, spec2.ratio, spec1.ratio2, normA, normC, normT, normG, dyebias,
    oobG, oob.ratio)


  ## Imputation
  for (colindex in 1:ncol(model.matrix)) {
    if(any(is.na(model.matrix[,colindex]))) {
      column <- model.matrix[,colindex]
      column[is.na(column)]    <- mean(column, na.rm = TRUE)
      model.matrix[ , colindex] <- column
    }
  }

  return(model.matrix)
}


.extractFromRGSet450k <- function(rgSet) {
  rgSet <- updateObject(rgSet)
  controlType <- c("BISULFITE CONVERSION I",
                   "BISULFITE CONVERSION II",
                   "EXTENSION",
                   "HYBRIDIZATION",
                   "NEGATIVE",
                   "NON-POLYMORPHIC",
                   "NORM_A",
                   "NORM_C",
                   "NORM_G",
                   "NORM_T",
                   "SPECIFICITY I",
                   "SPECIFICITY II",
                   "TARGET REMOVAL",
                   "STAINING")

  array <- minfi::annotation(rgSet)[["array"]]
  ## controlAddr <- getControlAddress(rgSet, controlType = controlType, asList = TRUE)
  ctrls <- minfi::getProbeInfo(rgSet, type = "Control")
  if(!all(controlType %in% ctrls$Type))
    stop("The `rgSet` does not contain all necessary control probes")
  ctrlsList <- S4Vectors::split(ctrls, ctrls$Type)[controlType]

  redControls <- minfi::getRed(rgSet)[ctrls$Address,,drop=FALSE]
  redControls <- lapply(ctrlsList, function(ctl) redControls[ctl$Address,,drop=FALSE])
  greenControls <- minfi::getGreen(rgSet)[ctrls$Address,,drop=FALSE]
  greenControls <- lapply(ctrlsList, function(ctl) greenControls[ctl$Address,,drop=FALSE])

  ## Extraction of the undefined negative control probes
  oobRaw <- minfi::getOOB(rgSet)
  probs <- c(0.01, 0.50, 0.99)
  greenOOB <- t(matrixStats::colQuantiles(oobRaw$Grn, na.rm = TRUE, probs = probs))
  redOOB   <- t(matrixStats::colQuantiles(oobRaw$Red, na.rm=TRUE,  probs = probs))
  oob      <- list(greenOOB = greenOOB, redOOB = redOOB)

  return(list(
    greenControls = greenControls,
    redControls = redControls,
    oob = oob, ctrlsList = ctrlsList,
    array = array))
}
