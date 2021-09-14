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
#' @references Nolan, Sean, et al. "A large-scale database of T-cell receptor beta (TCR??) sequences and binding associations from natural and synthetic exposure to SARS-CoV-2." (2020).
#' (\href{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7418738/}{PubMed})
#'
#' @source \href{https://clients.adaptivebiotech.com/pub/covid-2020}{ImmuneACCESS Data}
#'
#' @examples
#' data(covid_cdr3)
#' # Average CDR3 length
#' mean(nchar(covid_cdr3)) # [1] 43.56821
#' 
"covid_cdr3"