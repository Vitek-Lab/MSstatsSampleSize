# MSstatsSampleSize

MSstatsSampleSize uses as input a list of proteins quantified in mass spectrometry runs, and their annotations in terms of biological replicates and their membership in a group (such as a disease). The package fits intensity-based linear model on the input preliminary data. 

It estimates the protein abundance variance from the fitted model and simulates data with certain number of biological replicates based on the variance estimation. It reports the mean predictive accuracy of the classifier and mean protein importance over multiple iterations of the simulation. While varying the number of biological replicates to simulate, the sample size which generates the largest predictive accuracy is estimated. And the proteins which can best separate different conditions are reported.

The package also uses the fitted models and the fold changes estimated from the models to calculate sample size for hypothesis testing. It outputs the minimal number of biological replicates per condition to acquire the expected FDR and power under different fold changes.
