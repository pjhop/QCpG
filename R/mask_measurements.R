#' Mask measurements
#'
#' Mask measurements in a beta-matrix, based on a logical matrix that indicates which measurements
#' should be masked.
#'
#' @param beta beta-matrix or a filepath pointing to a txt-file where the first column contains the probe identifiers, and the remaining columns are samples.
#' @param mask_matrix logical matrix (TRUE/FALSE), where rows represent probes, and columns represent samples.
#' `TRUE` values indicate that the measurement should be masked.
#' Probes and samples present in the beta-matrix should also be present in this matrix.
#' A masking matrix can be generating by specifying the `path_failed_measurements` parameter in \code{\link{getQCmetrics}}.
#' @param output Ignored if `beta` is a beta-matrix, but should be specified if `beta` is a filepath.
#' In that case the beta-matrix will be loaded in chunks (of size `chunk_size`) and written to `output`.
#' @param chunk_size Relevant if `beta` is a filepath; specifies the number of probes to load at a time.
#' Defaults to 10,000.
#' @param delim Relevant if `beta` is a filepath, defaults to `\\t`.
#' @param show_progress Display a progress bar? (Relevant if `beta` is a filepath).
#' @export
#'
mask_measurements <- function(
  beta,
  mask_matrix,
  output = NULL,
  chunk_size = 10000,
  delim = "\t",
  show_progress = TRUE
) {
  ## Check if beta is a filepath or a matrix object
  if(is.matrix(beta)) {
    # Check if all probes and samples are present in mask_matrix
    check_rows <- mean(rownames(beta) %in% rownames(mask_matrix))
    check_cols <- mean(colnames(beta) %in% colnames(mask_matrix))
    if(check_rows < 1 || check_cols < 1) {
      if(check_rows < 1 && check_cols < 1) {
        stop("All probes and samples in the beta-matrix should be present in the masking matrix")
      } else if (check_rows < 1) {
        stop("All probes in the beta-matrix should be present in the masking matrix")
      } else {
        stop("All samples in the beta-matrix should be present in the masking matrix")
      }
    }

    mask_matrix_ <- mask_matrix[rownames(beta), colnames(beta)]
    message(sprintf("Setting %s out of %s measurements to missing",
                    sum(mask_matrix_),
                    ncol(mask_matrix_)*nrow(mask_matrix_)))
    beta[as.matrix(mask_matrix_)] <- NA
    return(beta)
  } else if(is.character(beta)) {
    if(is.null(output)) {stop("The `output` parameter should be specified when `beta` is a filepath.")}
    if(file.exists(output)) {stop(sprintf("%s already exists!", output))}

    readr::read_delim_chunked(
      beta,
      #callback = readr::SideEffectChunkCallback$new(.setNA.chunk),
      callback = readr::SideEffectChunkCallback$new(.setNA.chunk(mask_matrix=mask_matrix,output=output,chunk_size=chunk_size,delim=delim)),
      delim = delim,
      chunk_size = chunk_size,
      progress = show_progress    # optional, show progress bar
    )
  }
}

.setNA.chunk <- function(mask_matrix, output, chunk_size, delim) {
  function(df, pos) {
    probes <- df[,1] %>% dplyr::pull()
    df <- as.matrix(df[,2:ncol(df)])
    rownames(df) <- probes

    # Checks
    check_rows <- mean(rownames(df) %in% rownames(mask_matrix))
    check_cols <- mean(colnames(df) %in% colnames(mask_matrix))
    if(check_rows < 1 || check_cols < 1) {
      if(check_rows < 1 && check_cols < 1) {
        stop("All probes and samples in the beta-matrix should be present in the masking matrix")
      } else if (check_rows < 1) {
        stop("All probes in the beta-matrix should be present in the masking matrix")
      } else {
        stop("All samples in the beta-matrix should be present in the masking matrix")
      }
    }

    mask_matrix_ <- mask_matrix[rownames(df), colnames(df)]
    df[as.matrix(mask_matrix_)] <- NA
    df <- dplyr::bind_cols(
      Probe = rownames(df),
      df
    )
    is.append = ifelse(pos > 1, TRUE, FALSE)
    df %>% readr::write_delim(
      file = output,
      append = is.append,
      delim = delim
    )
  }
}

# .setNA.chunk = function(df, pos) {
#   probes <- df[,1] %>% dplyr::pull()
#   df <- as.matrix(df[,2:ncol(df)])
#   rownames(df) <- probes
#
#   # Checks
#   check_rows <- mean(rownames(df) %in% rownames(mask_matrix))
#   check_cols <- mean(colnames(df) %in% colnames(mask_matrix))
#   if(check_rows < 1 || check_cols < 1) {
#     if(check_rows < 1 && check_cols < 1) {
#       stop("All probes and samples in the beta-matrix should be present in the masking matrix")
#     } else if (check_rows < 1) {
#       stop("All probes in the beta-matrix should be present in the masking matrix")
#     } else {
#       stop("All samples in the beta-matrix should be present in the masking matrix")
#     }
#   }
#
#   mask_matrix_ <- mask_matrix[rownames(df), colnames(df)]
#   df[as.matrix(mask_matrix_)] <- NA
#   df <- dplyr::bind_cols(
#     Probe = rownames(df),
#     df
#   )
#   is.append = ifelse(pos > 1, T, F)
#   df %>% readr::write_delim(
#     file = output,
#     append = is.append,
#     delim = delim
#   )
# }
#
