# Names used in non-standard evaluation inside seqtrie_plot_graph() (ggplot2
# aes() and igraph layout columns); declared here to satisfy R CMD check.
utils::globalVariables(c("V1", "V2", "parent_x", "parent_y",
                         "child_x", "child_y", "fill", "node"))

.onLoad <- function(libname, pkgname) {
  S7::methods_register()
}

.onAttach <- function(libname, pkgname) {
  # nothing
}
