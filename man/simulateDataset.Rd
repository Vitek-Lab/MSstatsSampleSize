% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulateDataset.R
\name{simulateDataset}
\alias{simulateDataset}
\title{Simulate datasets with the given number of biological replicates and
proteins based on the input \emph{data}}
\usage{
simulateDataset(
  data,
  annotation,
  log2Trans = FALSE,
  num_simulations = 10,
  samples_per_group = 50,
  protein_select = list(Mean = "high", SD = "high"),
  protein_quantile_cutoff = list(Mean = 0, SD = 0),
  expected_FC = "data",
  list_diff_proteins = NULL,
  simulate_validation = FALSE,
  valid_samples_per_group = 50,
  ...
)
}
\arguments{
\item{data}{Protein abundance data matrix.
Rows are proteins and columns are biological replicates (samples).}

\item{annotation}{Group information for samples in data.
`Run' for MS run, `BioReplicate' for biological subject ID and `Condition' for group information are required.
`Run' information should be the same with the column of `data'. Multiple `Run' may come from same `BioReplicate'.}

\item{log2Trans}{Default is FALSE. If TRUE, the input `data' is log-transformed with base 2.}

\item{num_simulations}{Number of times to repeat simulation experiments
(Number of simulated datasets). Default is 10.}

\item{samples_per_group}{Number of samples per group to simulate. Default is 50.}

\item{protein_select}{select proteins with "low" or "high" mean abundance or
standard deviation (variance) or their combination to simulate.
If `protein_order = "mean"' or protein_order = "sd"', `protein_select' should be "low" or "high".
Default is "high", indicating high abundance or standard deviation proteins are selected to simulate.
If `protein_order = "combined"', `protein_select' has two elements.
The first element corrresponds to the mean abundance.
The second element corrresponds to the standard deviation (variance).
Default is c("high", "low") (select proteins with high abundance and low variance).}

\item{protein_quantile_cutoff}{Quantile cutoff(s) for selecting protiens to simulate.
For example, when `protein_rank="mean"', and
`protein_select="high"', `protein_quantile_cutoff=0.1'
Proteins are ranked based on their mean abundance across all the samples.
Then, the top 10% highest abundant proteins are selected to simulate.
Default is 0.0, which means that all the proteins are used.
If `protein_rank = "combined"', `protein_quantile_cutoff'` has two cutoffs.
The first element corrresponds to the cutoff for mean abundance.
The second element corrresponds to the cutoff for the standard deviation (variance).
Default is c(0.0, 1.0), which means that all the proteins will be used.}

\item{expected_FC}{Expected fold change of proteins.
The first option (Default) is "data",
indicating the fold changes are directly estimated from the input `data'.
The second option is a vector with predefined fold changes of listed proteins.
The vector names must match with the unique information of Condition in `annotation'.
One group must be selected as a baseline and has fold change 1 in the vector.
The user should provide list_diff_proteins, which users expect to have the fold changes greater than 1.
Other proteins that are not available in `list_diff_proteins' will be expected to have fold change = 1}

\item{list_diff_proteins}{Vector of proteins names
which are set to have fold changes greater than 1 between conditions.
If user selected `expected_FC= "data" ', this should be NULL.}

\item{simulate_validation}{Default is FALSE. If TRUE, simulate the validation set;
otherwise, the input `data' will be used as the validation set.}

\item{valid_samples_per_group}{Number of validation samples per group to simulate.
This option works only when user selects `simulate_validation=TRUE'. Default is 50.}

\item{protein_rank}{The standard to rank the proteins in the input `data'.
It can be 1) "mean" of protein abundances over all the samples or
2) "sd" (standard deviation) of protein abundances over all the samples or
3) the "combined" of mean abundance and standard deviation.
The proteins in the input `data' are ranked based on `protein_rank'
and the user can select a subset of proteins to simulate.}
}
\value{
\emph{num_proteins} is the number of simulated proteins.
 It should be set up by parameters, named \emph{protein_proportion} or \emph{protein_number}

\emph{num_samples} is a vector with the number of simulated samples in each condition.
It should be same as the parameter, \emph{samples_per_group}

\emph{input_X} is the input protein abundance matrix `data'.

\emph{input_Y} is the condition vector for the input `data.

\emph{simulation_train_Xs} is the list of simulated protein abundance matrices.
Each element of the list represents one simulation.

\emph{simulation_train_Ys} is the list of simulated condition vectors.
Each element of the list represents one simulation.

\emph{valid_X} is the validation protein abundance matrix, which is used for classification.

\emph{valid_Y} is the condition vector of validation samples.
}
\description{
Simulate datasets with the given number of biological replicates and
proteins based on the input \emph{data}
}
\details{
This function simulate datasets with
the given numbers of biological replicates and
proteins based on the input dataset (input for this function).
The function fits intensity-based linear model on the input \emph{data}
in order to get variance and mean abundance, using \code{\link{estimateVar}} function.
Then it uses variance components and mean abundance to simulate new training data
with the given sample size and protein number.
It outputs the number of simulated proteins,
a vector with the number of simulated samples in a condition,
the list of simulated training datasets,
the input preliminary dataset and
the (simulated) validation dataset.
}
\examples{
data(OV_SRM_train)
data(OV_SRM_train_annotation)

# num_simulations = 10: simulate 10 times
# protein_rank = "mean", protein_select = "high", and protein_quantile_cutoff = 0.0:
# select the proteins with high mean abundance based on the protein_quantile_cutoff
# expected_FC = "data": fold change estimated from OV_SRM_train
# simulate_validation = FALSE: use input OV_SRM_train as validation set
# valid_samples_per_group = 50: 50 samples per condition
simulated_datasets <- simulateDataset(data = OV_SRM_train,
                                      annotation = OV_SRM_train_annotation,
                                      log2Trans = FALSE,
                                      num_simulations = 10,
                                      samples_per_group = 50,
                                      protein_rank = "mean",
                                      protein_select = "high",
                                      protein_quantile_cutoff = 0.0,
                                      expected_FC = "data",
                                      list_diff_proteins =  NULL,
                                      simulate_validation = FALSE,
                                      valid_samples_per_group = 50)

# the number of simulated proteins
simulated_datasets$num_proteins

# a vector with the number of simulated samples in each condition
simulated_datasets$num_samples

# the list of simulated protein abundance matrices
# Each element of the list represents one simulation
head(simulated_datasets$simulation_train_Xs[[1]]) # first simulation

# the list of simulated condition vectors
# Each element of the list represents one simulation
head(simulated_datasets$simulation_train_Ys[[1]]) # first simulation

}
\author{
Ting Huang, Meena Choi, Sumedh Sankhe, Olga Vitek.
}
