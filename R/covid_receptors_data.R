#' Adaptive COVID TCRB nucleotide rearrangements
#'
#' Unique TCRB nucleotide rearrangement sequences (the full "TCR Nucleotide
#' Sequence", spanning the V segment, CDR3 junction, and J segment).
#' Every sequence is 87 nucleotides long. Sequences containing the ambiguous base
#' `N` were removed, and duplicates were collapsed. The data are licensed under
#' the Creative Commons Attribution 4.0 International License.
#'
#' @docType data
#'
#' @usage data(covid_receptors)
#'
#' @format A character vector of length 139,667, each a 87-nucleotide DNA string.
#'
#' @keywords datasets
#'
#' @references Nolan, Sean, et al. "A large-scale database of T-cell receptor beta (TCRB) sequences and binding associations from natural and synthetic exposure to SARS-CoV-2." (2020). doi: 10.21203/rs.3.rs-51964/v1.
#'
#' @examples
#' data(covid_receptors)
#' # All sequences are the same length
#' table(nchar(covid_receptors)) # 87: 139667
#'
"covid_receptors"
