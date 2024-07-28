#' SNPlinker file used in MinE QC
#'
#' Probes used for inferring genotypes in the MinE data.
#' Probes that contain a SNP (MAF > 0.05 in EUR and AMR population) within the 'G' or 'C' base
#' were selected. Furthermore, probes that are present on both the EPIC and the 450k array were selected.
#' Then the beta2genotype() function from the omicsPrint package was used within the
#' MinE 450k data (6556 samples) to select probes that could be accurately called from these probes (426 probes in total).
#'
#' @docType data
#'
#' @usage data(snplinker_mine)
#'
#' @format An object of class "tibble".
#'
#' @keywords datasets
#'
#' @examples
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
#' data(snplinker_mine)
#' metrics <- getQCmetrics(samplesheet, snplinker = snplinker_mine)
"snplinker_mine"
