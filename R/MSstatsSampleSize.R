#' MSstatsSampleSize: A package for optimal design of high-dimensional MS-based proteomics experiment
#'
#' A set of functions for sample size calculation. The packages estimates the variance in the input protein abundance data
#' and simulates data with pre-defined number of biological replicates based on the variance estimation.
#' It reports the mean predictive accuracy of the classifier and mean protein importance over multiple iterations of the simulation.
#'
#' @section functions :
#' \itemize{
#'   \item \code{\link{estimateVar}} : estimate the mean abundance and variance of each protein in each condition.
#'   \item \code{\link{meanSDplot}} : draw the plot for the mean protein abundance vs standard deviation in each condition.
#'   \item \code{\link{simulateDataset}} : simulate datasets with the pre-defined size based on the preliminary data.
#'   \item \code{\link{designSampleSizeClassification}} : estimate the mean predictive accuracy and protein importance over all the simulated datasets.
#'   \item \code{\link{designSampleSizePCAplot}} : make PCA plots with the first two components for each simulated dataset.
#'   \item \code{\link{designSampleSizeClassificationPlots}} : visualization for sample size calculation in classification.
#'   \item \code{\link{designSampleSizeHypothesisTestingPlot}} : Sample size calculation plot for hypothesis testing.
#' }
#'
#' @docType package
#' @name MSstatsSampleSize
NULL
