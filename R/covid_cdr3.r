#' Adaptive COVID TCRB CDR3 data
#'
#' Unique TCRB CDR3 sequences from the Nolan et al. 2020. CDR3s were extracted via IgBLAST. The license for this data is Creative Commons Attribution 4.0 International License. 
#'
#' @docType data
#'
#' @usage data(covid_cdr3)
#'
#' @format A character vector of length 133,034. 
#'
#' @keywords datasets
#'
#' @references Nolan, Sean, et al. "A large-scale database of T-cell receptor beta (TCRB) sequences and binding associations from natural and synthetic exposure to SARS-CoV-2." (2020). doi: 10.21203/rs.3.rs-51964/v1. 
#'
#' @examples
#' data(covid_cdr3)
#' # Average CDR3 length
#' mean(nchar(covid_cdr3)) # [1] 43.56821
#' 
"covid_cdr3"