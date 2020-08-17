
#' The training set from a study for subjects with ovarian cancer
#'
#' It is a protein abundance data matrix,
#' where rows are proteins and columns are samples.
#' It includes log2 protein intensities for 67 proteins
#' among 173 biological subjects from control and cancer groups.
#' It is the input for estimateVar and simulateDataset function,
#' with annotation file. It should be prepared by users.
#'
#' @references Huttenhain R and Choi M et al. (2019).
#'  A targeted mass spectrometry strategy for developing proteomic biomarkers:
#'  a case study of epithelial ovarian cancer.
#'  Mol Cell Proteomics 18(9):1836-1850.
#'  doi:10.1074/mcp.RA118.001221.
#' @format A numeric matrix with 67 rows and 173 columns.
#' @examples
#' head(OV_SRM_train)
#'
"OV_SRM_train"

#' Annotation file for \code{\link{OV_SRM_train}},
#'
#' Annotation of example data, \code{\link{OV_SRM_train}}, in this package.
#' It should be prepared by users. The variables are as follows:
#'
#' \itemize{
#'  \item Run : MS run ID. It should be the same as the column names of \code{\link{OV_SRM_train}}.
#'  Multiple `Run' may come from same `BioReplicate'.
#'  \item BioReplicate : Unique ID for biological subject.
#'  \item Condition : Condition for BioReplicate (ex. Healthy, Cancer, Time0)
#' }
#'
#' @references Huttenhain R and Choi M et al. (2019).
#'  A targeted mass spectrometry strategy for developing proteomic biomarkers:
#'  a case study of epithelial ovarian cancer.
#'  Mol Cell Proteomics 18(9):1836-1850.
#'  doi:10.1074/mcp.RA118.001221.
#' @format A data frame with 173 rows and 3 variables.
#' @examples
#' head(OV_SRM_train_annotation)
#'
"OV_SRM_train_annotation"

#' Example of output from estimateVar function
#'
#' It is the output of \code{\link{estimateVar}} function with two inputs:
#' \code{\link{OV_SRM_train}} and \code{\link{OV_SRM_train_annotation}}.
#' The list should include the required elements as below.
#'
#' \itemize{
#'   \item data : the input data matrix with log2 protein abundance.
#'   \item model : the list of linear models trained for each protein.
#'   \item mu : the mean abundance matrix of each protein in each condition
#'   \item sigma : the standard deviation matrix of each protein
#'   in each condition
#'   \item promean: the mean abundance vector of each protein
#'   across all the samples
#'   \item promean: the standard deviation vector of each protein
#'   across all the samples.
#'   \item protein : proteins, correpsonding to the rows
#'   in \emph{mu} and \emph{sigma} or the element of \emph{promean}
#' }
#'
#' @format A list with seven elements
#' @examples
#' head(variance_estimation$mu)
#' head(variance_estimation$sigma)
#' head(variance_estimation$promean)
#'
"variance_estimation"

#' Example of output from simulateDataset function
#'
#' It is the output of \code{\link{simulateDataset}} function with two inputs:
#' \code{\link{OV_SRM_train}} and \code{\link{OV_SRM_train_annotation}}.
#' The list should include the required elements as below.
#'
#' \itemize{
#'   \item num_proteins : the number of simulated proteins
#'   \item num_samples : a vector with the number of simulated samples
#'   in each condition
#'   \item simulation_train_Xs : the list of simulated protein abundance matrices.
#'   Each element of the list represents one simulation
#'   \item simulation_train_Ys : the list of simulated condition vectors(simulation_train_Xs).
#'   Each element of the list represents one simulation
#'   \item input_X : the input protein abundance matrix `OV_SRM_train'.
#'   \item input_Y : is the condition vector for the input `OV_SRM_train'.
#'   \item valid_X: the validation protein abundance matrix, which is used for classification
#'   \item valid_Y : the condition vector of validation samples (valid_X)
#' }
#'
#' @format A list with eight elements
#' @examples
#' simulated_datasets$num_proteins
#' simulated_datasets$num_samples
#' head(simulated_datasets$simulation_train_Xs[[1]])
#' head(simulated_datasets$simulation_train_Ys[[1]])
#'
"simulated_datasets"

#' Example of output from designSampleSizeClassification function
#'
#' It is the output of \code{\link{designSampleSizeClassification}} function
#' with a list of \code{\link{simulated_datasets}}
#' generated under same protein number and sample size.
#' The list should include the required elements as below.
#'
#' \itemize{
#'   \item num_proteins : the number of simulated proteins
#'   \item num_samples : a vector with the number of simulated samples in each condition
#'   \item mean_predictive_accuracy : the mean predictive accuracy over all the simulated datasets.
#'   \item mean_feature_importance : the mean protein importance vector over all the simulated datasets,
#'   the length of which is `num_proteins'.
#'   \item predictive_accuracy : a vector of predictive accuracy on each simulated dataset.
#'   \item feature_importance : a matrix of feature importance,
#'   where rows are proteins and columns are simulated datasets.
#'   the length of which is `num_proteins'.
#' }
#'
#' @format A list with six elements
#' @examples
#' classification_results$num_proteins
#' classification_results$num_samples
#' classification_results$mean_predictive_accuracy
#' head(classification_results$mean_feature_importance)
#'
"classification_results"
