#' Estimate the mean predictive accuracy and mean protein importance over all the simulated datasets
#' @details This function fits the classification model,
#' in order to classify the subjects in each simulated training dataset (the output of \code{\link{simulateDataset}}).
#' Then the fitted model is validated on the (simulated) validation set (the output of \code{\link{simulateDataset}}).

#' Two performance are reported :
#'
#' (1) the mean predictive accuracy : The function trains classifier on each simulated training dataset and
#' reports the predictive accuracy of the trained classifier on the validation data (output of \code{\link{simulateDataset}} function).
#' Then these predictive accuracies are averaged over all the simulation.
#'
#' (2) the mean protein importance : It represents the importance of a protein in separating different groups.
#' It is estimated on each simulated training dataset using function `varImp' from package caret. Please refer to the help file of `varImp' about
#' how each classifier calculates the protein importance. Then these importance values for each protein are averaged over all the simulation.
#'
#' The list of classification models trained on each simulated dataset, the predictive accuracy
#' on the validation set predicted by the corresponding classification model and
#' the importance value for all the proteins estimated by the corresponding classification model
#' are also reported.
#'
#' @param simulations A list of simulated datasets It should be the name of the output of \code{\link{simulateDataset}} function.
#' @param classifier A string specifying which classfier to use. This function uses function `train' from package caret.
#' The options are 1) rf (random forest calssifier, default option). 2) nnet (neural network),
#' 3) svmLinear (support vector machines with linear kernel), 4) logreg (logistic regression), and 5) naive_bayes (naive_bayes).
#' @param top_K the number of proteins selected as important features (biomarker candidates).
#' All the proteins are ranked in descending order based on its importance to separate different groups and
#' the `top_K` proteins are selected as important features.
#' @param parallel Default is FALSE. If TRUE, parallel computation is performed.
#'
#' @return \emph{num_proteins} is the number of simulated proteins.
#'  It should be the same as one of the output from \emph{simulateDataset}, called \emph{num_proteins}
#' @return \emph{num_samples} is a vector with the number of simulated samples in each condition.
#'  It should be the same as one of the output from \emph{simulateDataset}, called \emph{num_samples}
#' @return \emph{mean_predictive_accuracy} is the mean predictive accuracy over all the simulated datasets,
#' which have same `num_proteins' and `num_samples'.
#' @return \emph{mean_feature_importance} is the mean protein importance vector over all the simulated datasets,
#' the length of which is `num_proteins'.
#' @return \emph{predictive_accuracy} is a vector of predictive accuracy on each simulated dataset.
#' @return \emph{feature_importance} is a matrix of feature importance, where rows are proteins and columns are simulated datasets.
#' @return \emph{results} is the list of classification models trained on each simulated dataset and
#' the predictive accuracy on the validation set predicted by the corresponding classification model.
#' @author Ting Huang, Meena Choi, Olga Vitek
#'
#' @examples data(OV_SRM_train)
#' data(OV_SRM_train_annotation)
#'
#' # num_simulations = 10: simulate 10 times
#' # expected_FC = "data": fold change estimated from OV_SRM_train
#' # select_simulated_proteins = "proportion":
#' # select the simulated proteins based on the proportion of total proteins
#' # simulate_valid = FALSE: use input OV_SRM_train as validation set
#' # valid_samples_per_group = 50: 50 samples per condition
#' simulated_datasets <- simulateDataset(data = OV_SRM_train,
#'                                       annotation = OV_SRM_train_annotation,
#'                                       num_simulations = 10,
#'                                       expected_FC = "data",
#'                                       list_diff_proteins =  NULL,
#'                                       select_simulated_proteins = "proportion",
#'                                       protein_proportion = 1.0,
#'                                       protein_number = 1000,
#'                                       samples_per_group = 50,
#'                                       simulate_valid = FALSE,
#'                                       valid_samples_per_group = 50)
#'
#' # run classification on simulated datasets without parallel computation
#' classification_results <- designSampleSizeClassification(simulations = simulated_datasets,
#'                                                          parallel = FALSE)
#'
#  # the number of simulated proteins
#' classification_results$num_proteins
#'
#' # a vector with the number of simulated samples in each condition
#' classification_results$num_samples
#'
#' # the mean predictive accuracy over all the simulated datasets,
#' # which have same 'num_proteins' and 'num_samples'
#' classification_results$mean_predictive_accuracy
#'
#' # the mean protein importance vector over all the simulated datasets,
#' # the length of which is 'num_proteins'.
#' head(classification_results$mean_feature_importance)
#'
#' @importFrom caret train trainControl varImp predict.train
#' @importFrom BiocParallel bplapply
#' @importFrom stats rnorm predict
#' @importFrom utils sessionInfo read.table write.table txtProgressBar setTxtProgressBar
#' @export
designSampleSizeClassification <- function(simulations,
                                           classifier = "rf",
                                           top_K = 10,
                                           parallel = FALSE) {

    ###############################################################################
    ## log file
    ## save process output in each step
    loginfo <- .logGeneration()
    finalfile <- loginfo$finalfile
    processout <- loginfo$processout

    processout <- rbind(processout,
                        as.matrix(c(" ", " ", "MSstatsSampleSize - designSampleSizeClassification function", " "), ncol=1))

    ###############################################################################
    ## Input and option checking

    ## 1. input  should be the output of SimulateDataset function
    if ( !is.element('simulation_train_Xs', names(simulations)) | !is.element('simulation_train_Ys', names(simulations)) ) {

        processout <- rbind(processout,
                            "The required input - simulations : did not process from SimulateDataset function. - stop")
        write.table(processout, file=finalfile, row.names=FALSE)

        stop("Please use 'SimulateDataset' first. Then use output of SimulateDataset function as input in designSampleSizeClassification.")
    }

    ## 2. input for classifier option
    if ( !any(classifier == c('rf', 'nnet', 'svmLinear', 'logreg', 'naive_bayes')) ) {
        processout <- rbind(processout, c("ERROR: `classifier` should be one of 'rf', 'nnet', 'svmLinear', 'logreg', and 'naive_bayes'. Please check it."))
        write.table(processout, file=finalfile, row.names=FALSE)

        stop("`classifier` should be one of 'rf', 'nnet', 'svmLinear', 'logreg', and 'naive_bayes'. Please check it. \n")
    }
    processout <- rbind(processout, c(paste0("classifier : ", classifier)))
    write.table(processout, file=finalfile, row.names=FALSE)
    message(" classifier: ", classifier)

    ## 3. input for top_K option
    if (is.null(top_K) ) {
        processout <- rbind(processout, c("ERROR : top_K is required. Please provide the value for top_K."))
        write.table(processout, file=finalfile, row.names=FALSE)

        stop("top_K is required. Please provide the value for top_K. \n")

    } else if ( top_K < 0 | top_K > simulations$num_proteins ) {
        processout <- rbind(processout,
                            c(paste0("ERROR : top_K should be between 0 and the total number of protein(", simulations$num_proteins,
                                     "). Please check the value for top_K")))
        write.table(processout, file=finalfile, row.names=FALSE)

        stop(paste0("top_K should be between 0 and the total number of protein(", simulations$num_proteins,
                    "). Please check the value for top_K \n"))
    }
    processout <- rbind(processout, c(paste0("top_K = ", top_K)))
    write.table(processout, file=finalfile, row.names=FALSE)


    ###############################################################################
    ## start to train classifier

    processout <- rbind(processout, c(" Start to train classifier..."))
    write.table(processout, file=finalfile, row.names=FALSE)
    message(" Start to train classifier...")

    ## get the validation set for prediction
    iter <- length(simulations$simulation_train_Xs) # number of simulations
    num_proteins <- simulations$num_proteins
    num_samples <- simulations$num_samples
    valid_x <- simulations$valid_X
    valid_y <- simulations$valid_Y

    ## if parallel TRUE,
    if(parallel){

        # ## create cluster for paralleled workflow
        # message(paste0("Cluster Size: ", threads,"\n"))
        #
        # ## allocate resource for parallel computation
        # param <- SnowParam(workers = threads, type = "SOCK")

        # ## fit the classifier for each simulation dataset
        # results <- bplapply(1:iter, .classificationPerformance,
        #                     classifier=classifier,
        #                     train_x_list = simulations$simulation_train_Xs,
        #                     train_y_list = simulations$simulation_train_Ys,
        #                     valid_x = valid_x,
        #                     valid_y = valid_y,
        #                     BPPARAM=param,
        #                     top_K = top_K)

        ## fit the classifier for each simulation dataset
        results <- bplapply(seq_len(iter), .classificationPerformance,
                            classifier=classifier,
                            train_x_list = simulations$simulation_train_Xs,
                            train_y_list = simulations$simulation_train_Ys,
                            valid_x = valid_x,
                            valid_y = valid_y,
                            top_K = top_K)


    } else { ## if parallel FALSE,
        # ## show progress
        # pb <- txtProgressBar(max = iter, style = 3)
        # ## progress
        # setTxtProgressBar(pb, i)
        # close(pb)

        ## fit the classifier for each simulation dataset
        results <- lapply(seq_len(iter),
                          .classificationPerformance,
                          classifier=classifier,
                          train_x_list = simulations$simulation_train_Xs,
                          train_y_list = simulations$simulation_train_Ys,
                          valid_x = valid_x,
                          valid_y = valid_y,
                          top_K = top_K)

    }

    ## calculate the mean predictive accuracy over all the simulations
    PA <- NULL

    ## calculate the frequency a protein is selected as important (biomarker candidates)
    FI <- NULL
    features <- rownames(varImp(results[[1]]$model, scale = TRUE)$importance)

    for (i in seq_len(iter)) {
        # record the importance of each protein
        PA <- c(PA, results[[i]]$acc)

        # record the importance of each protein
        imp <- varImp(results[[i]]$model, scale = TRUE)$importance
        # select the top-k important proteins
        imp.prots <- rownames(imp)[order(imp, decreasing=TRUE)][seq_len(top_K)]
        # if important, set 1; otherwise,set 0
        imp$Important <- ifelse(rownames(imp) %in% imp.prots, 1, 0)

        FI <- cbind(FI, imp[features, "Important"])

    }

    ## report the training and validating done
    processout <- rbind(processout, c(" Finish to train classifier and to check the performance."))
    write.table(processout, file=finalfile, row.names=FALSE)
    message(" Finish to train classifier and to check the performance.")

    # assign simulation index
    simulation_index <- paste0("simulation", seq_len(iter))
    rownames(FI) <- features
    colnames(FI) <- simulation_index
    names(PA) <- simulation_index

    # calculate mean predictive accuracy
    meanPA <-  mean(PA)

    # calculate mean feature importance
    meanFI <-  rowSums(FI)
    names(meanFI) <- features
    # sort in descending order
    meanFI <- sort(meanFI, decreasing=TRUE)

    processout <- rbind(processout, c(" Report the mean predictive accuracy and mean feature importance."))
    write.table(processout, file=finalfile, row.names=FALSE)
    message(" Report the mean predictive accuracy and mean feature importance.")

    return(list(num_proteins = num_proteins, # number of proteins
                num_samples = num_samples, # number of samples per group
                mean_predictive_accuracy = meanPA, # mean predictive accuracy
                mean_feature_importance = meanFI, # vector of mean feature importance
                predictive_accuracy = PA, # vector of predictive accuracy
                feature_importance = FI,  # matrix of feature importance
                results = results # fitted models
                ))

}

