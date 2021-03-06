% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Methods.R
\name{.feature_importance}
\alias{.feature_importance}
\title{Feature Importance}
\usage{
.feature_importance(model, classifier, top_K)
}
\arguments{
\item{model}{A model object derived from the caret classification models}

\item{classifier}{The names of the classifier used to train the model on}

\item{top_K}{an Integer specifying the number of features to consider for 
selection}
}
\value{
A string vector of names of the top_K features
}
\description{
The method .feature_importance() finds variables importances
for the given model and selects the top n features
}
\keyword{internal}
