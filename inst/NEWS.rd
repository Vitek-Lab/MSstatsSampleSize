\name{MSstatsSampleSizenews}
\title{MSstatsSampleSize News}
\encoding{UTF-8}

\section{Version 1.2.2, Bioconductor 3.11 Release (July 15th 2020)}{\itemize{
\item
update the log info system
\item
designSampleSizeClassification(): fix the bugs in classifiers nnet, logreg, naive_bayes and svmLinear
}}

\section{Version 1.2.1, Bioconductor 3.11 Release (July 2020)}{\itemize{
\item
estimateVar() and simulateDataset(): Add log2Trans option
\item
simulateDataset(): change the standard to select proteins to simulate
\item
meanSDplot(): change the plot for the mean protein abundance (X-axis) vs standard deviation (Y-axis) across all the samples.
}}

\section{Version 1.0.2, Bioconductor 3.10 Release (March 2020)}{\itemize{
\item
Add MS run column to the input annotation file in order to cover more designs
\item
Report an optimal sample size for classification
}}

\section{Version 1.0.1, Bioconductor 3.10 Release (February 2020)}{\itemize{
\item
Add new function designSampleSizeHypothesisTestingPlot, which estiamtes the optimal sample size for hypothesis testing.
\item
Add new multivariate t-distribution distribution ellipse to pca plot.
\item
Calculate the predictive accuracy on the selected set of biomarker candidates, instead of all the proteins
\item
Re-define protein importance as the number of iterations that a protein is selected as predictive
}}

\section{Version 1.0.0, Bioconductor 3.10 Release (October 2019)}{\itemize{
\item
New package MSstatsSampleSize, for optimal design of high-dimensional MS-based proteomics experiment.
}}
