#' Khan: Feature Selection on Graphs with FDR Control
#'
#' Implements the methodology from Liang, Zhang, and Neykov (2024) for
#' false-discovery-rate controlled inference on subgraph features and
#' persistent homology generators in Gaussian graphical models and
#' Ising models. The package provides:
#'
#' * Debiased GLasso estimation and edge-level p-values
#'   ([debias_glasso()], [compute_ggm_pvalues()])
#' * BHq selection for arbitrary subgraph features ([bhq_test()])
#' * The Discrete Gram-Schmidt algorithm for homology generators
#'   ([dgs()])
#' * The KHAN algorithm for uFDR-controlled persistent homology
#'   ([khan()])
#' * Gaussian multiplier bootstrap baselines
#'   ([bootstrap_feature_test()], [bootstrap_simultaneous_test()])
#' * Subgraph enumeration utilities ([find_all_features()])
#'
#' @keywords internal
"_PACKAGE"
