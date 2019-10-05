#' Estimate the mean predictive accuracy and protein importance over all the simulated datasets
#' @details This function fits the statistical model for classification,
#' in order to classify the subjects in the simulated training datasets (the output of \code{\link{simulateDataset}}).
#' Then the fitted model is validated by the (simulated) validation set (the output of \code{\link{simulateDataset}}).

#' Two performance are reported :
#'
#' (1) the mean predictive accuracy : The function trains classifier on each simulated training dataset and
#' reports the predictive accuracy of the trained classifier on the validation data (output of \code{\link{simulateDataset}} function).
#' Then these predictive accuracies are averaged over all the simulation.
#'
#' (2) the protein importance : It represents the importance of a protein in separating different groups.
#' It is estimated using function `varImp' from package caret. Please refer to the help file of `varImp' about
#' how each classifier calculates the protein importance.
#'
#' The sample size per condition, which generates the largest predictive accuracy, is calculated,
#' while varying the number of biological replicates to simulate.
#' Also, the proteins, which can classify the conditions best, are reported.
#' The reported sample size per condition can be used to design future experiments.
#' The list of classification models trained on each simulated dataset and the predictive accuracy
#' on the validation set predicted by the corresponding classification model is also reported.
#'
#' @param simulations A list of simulated datasets It should be the name of the output of \code{\link{simulateDataset}} function.
#' @param classifier A string specifying which classfier to use. This function uses function `train' from package caret.
#' The options are 1) rf (random forest calssifier, default option). 2) nnet (neural network),
#' 3) svmLinear (support vector machines with linear kernel), 4) logreg (logistic regression), and 5) naive_bayes (naive_bayes).
#' @param threads A user needs to specify the number of threads (clusters) for computation. Default is NULL. User can specify it. For example, `threads=4'.
#'
#' @return \emph{num_proteins} the number of simulated proteins.
#' @return \emph{num_samples} a vector with the number of simulated samples in each condition.
#' @return \emph{results} the list of classification models trained on each simulated dataset and
#' the predictive accuracy on the validation set predicted by the corresponding classification model.
#' @return \emph{mean_predictive_accuracy} the mean predictive accuracy over all the simulated datasets,
#' which have same `num_proteins' and `num_samples'.
#' @return \emph{mean_feature_importance} the mean protein importance vector over all the simulated datasets,
#' the length of which is `num_proteins'.
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
#' classification_results <- designSampleSizeClassification(simulations = simulated_datasets)
#'
# the number of simulated proteins
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
#' @import doParallel
#' @importFrom caret train trainControl varImp predict.train
#' @importFrom foreach foreach %dopar%
#' @importFrom stats rnorm predict
#' @importFrom utils sessionInfo read.table write.table txtProgressBar setTxtProgressBar
#' @export
designSampleSizeClassification <- function(simulations,
                                           classifier = "rf",
                                           threads = NULL) {

    ###############################################################################
    ## log file
    ## save process output in each step

    allfiles <- list.files()
    filenaming <- "MSstatsSampleSize-ProgressReport"

    if (length(grep(filenaming, allfiles)) == 0) {

        finalfile <- "MSstatsSampleSize-ProgressReport.log"

        session <- sessionInfo()
        sink("sessionInfo.txt")
        print(session)
        sink()

        processout <- as.matrix(read.table("sessionInfo.txt", header=TRUE, sep="\t"))

    } else {

        num <- 0
        finalfile <- "MSstatsSampleSize-ProgressReport.log"

        while (is.element(finalfile, allfiles)) {
            num <- num + 1
            lastfilename <- finalfile ## in order to rea
            finalfile <- paste0(paste(filenaming, num, sep="-"), ".log")
        }

        finalfile <- lastfilename
        processout <- as.matrix(read.table(finalfile, header=TRUE, sep="\t"))
    }

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
        processout <- rbind(processout, c("ERROR : `classifier` should be one of 'rf', 'nnet', 'svmLinear', 'logreg', and 'naive_bayes'. Please check it."))
        write.table(processout, file=finalfile, row.names=FALSE)

        stop("ERROR : `classifier` should be one of 'rf', 'nnet', 'svmLinear', 'logreg', and 'naive_bayes'. Please check it. \n")
    }

    processout <- rbind(processout, c(paste0("classifier : ", classifier)))
    write.table(processout, file=finalfile, row.names=FALSE)
    message(" classifier: ", classifier)



    ###############################################################################
    ## start to train classifier

    processout <- rbind(processout, c(" Start to train classifier..."))
    write.table(processout, file=finalfile, row.names=FALSE)
    message(" Start to train classifier...")

    ## get the validation set for prediction
    iter <- length(simulations$simulation_train_Xs) # number of simulations
    num_proteins <- ncol(simulations$simulation_train_Xs[[1]])
    num_samples <- table(simulations$simulation_train_Ys[[1]])
    valid_x <- simulations$valid_X
    valid_y <- simulations$valid_Y

    ## if cluster available,
    if(!is.null(threads)){

        ## create cluster for paralleled workflow
        message(paste0("Cluster Size: ", threads,"\n"))

        ## allocate resource for parallel computation
        cl <- parallel::makeCluster(threads, outfile="")
        doParallel::registerDoParallel(cl)

        ## fit the classifier for each simulation dataset
        results <- foreach::foreach(i=1:iter) %dopar% {

            x <- as.data.frame(simulations$simulation_train_Xs[[i]])
            y <- as.factor(simulations$simulation_train_Ys[[i]])

            start_time <- Sys.time()

            if(classifier == "logreg"){
                ## Train logistic regression on training data
                model <- caret::train(x=x, y=make.names(y),
                                      method = classifier,
                                      trControl = caret::trainControl(method = "none", classProbs = TRUE),
                                      # trControl = trainControl(method = "cv", number=10, classProbs = TRUE),
                                      method = "glm",
                                      family = "binomial",
                                      verbose=TRUE)

            } else {
                ## set the parameters of classifier
                switch(classifier,
                   "rf"= {tunegrid=data.frame(mtry=2)},
                   "nnet"={tunegrid=data.frame(size=5, decay=0.1)},
                   "svmLinear"={tunegrid=data.frame(C=1)},
                   "naive_bayes"={tunegrid=data.frame(laplace=0, usekernel=FALSE, adjust=1)}
                )

                ## Train classifier on training data
                model <- caret::train(x=x, y=make.names(y),
                                  method = classifier,
                                  trControl = caret::trainControl(method = "none", classProbs = TRUE),
                                  # trControl = trainControl(method = "cv", number=10, classProbs = TRUE),
                                  tuneGrid = tunegrid,
                                  verbose=TRUE)
            }

            end_time <- Sys.time()
            runtime <- end_time - start_time
            message("    Took ", runtime, " seconds")

            ## Predict validation data
            model.pred <- predict(model, valid_x)

            ## Calculate predictive accuracy on validation data
            acc <- sum(diag(table(model.pred, valid_y))) / length(model.pred)

            ## record the predictive accuracy and train model
            run <- list(acc, model)
        }
        parallel::stopCluster(cl)

        ## calculate the mean predictive accuracy over all the simulations
        PA <- NULL

        ## calculate the mean protein importance over all the simulations
        FI <- NULL
        features <- rownames(varImp(results[[1]][[2]], scale = TRUE)$importance)

        for (i in 1:iter) {

          PA <- c(PA, results[[i]][[1]])

          # calculate the importance of each protein
          FI <- cbind(FI, varImp(results[[i]][[2]], scale = TRUE)$importance[features, ])
        }

    } else {
        ## if threads==NULL, no clustering for computing

        ## to record the predictive accuracy over all the simulations
        PA <- NULL

        ## to record the protein importance over all the simulations
        FI <- NULL
        features <- NULL

        ## show progress
        pb <- txtProgressBar(max = iter, style = 3)

        results <- list()
        ## fit the classifier for each simulation dataset
        for(i in 1:iter) {

            x <- as.data.frame(simulations$simulation_train_Xs[[i]])
            y <- as.factor(simulations$simulation_train_Ys[[i]])


            if(classifier == "logreg"){

                ## Train logistic regression on training data
                model <- caret::train(x=x, y=make.names(y),
                                method = classifier,
                                trControl = caret::trainControl(method = "none", classProbs = TRUE),
                                # trControl = trainControl(method = "cv", number=10, classProbs = TRUE),
                                method = "glm",
                                family = "binomial",
                                verbose=TRUE)

            } else {
                ## set the parameters of classifier
                switch(classifier,
                 "rf"= {tunegrid=data.frame(mtry=2)},
                 "nnet"={tunegrid=data.frame(size=5, decay=0.1)},
                 "svmLinear"={tunegrid=data.frame(C=1)},
                 "naive_bayes"={tunegrid=data.frame(laplace=0, usekernel=FALSE, adjust=1)}
                )

                ## Train classifier on training data
                model <- caret::train(x=x, y=make.names(y),
                                method = classifier,
                                trControl = caret::trainControl(method = "none", classProbs = TRUE),
                                # trControl = trainControl(method = "cv", number=10, classProbs = TRUE),
                                tuneGrid = tunegrid,
                                verbose=TRUE)
            }

            ## Predict validation data
            model.pred <- predict(model, valid_x)

            ## Calculate predictive accuracy on validation data
            acc <- sum(diag(table(model.pred, valid_y))) / length(model.pred)

            ## record the predictive accuracy
            PA <- c(PA, acc)

            ## record the train model
            ## first, keep the protein list, in case that the protein list is not consistent
            if (i == 1) {
                features <- rownames(varImp(model, scale = TRUE)$importance)
            }

            FI <- cbind(FI, varImp(model, scale = TRUE)$importance[features, ])

            results[[i]] <- list(acc, model)

            ## progress
            setTxtProgressBar(pb, i)
        }
        close(pb)
    }

    ## report the training and validating done
    processout <- rbind(processout, c(" Finish to train classifier and to check the performance."))
    write.table(processout, file=finalfile, row.names=FALSE)
    message(" Finish to train classifier and to check the performance.")

    # calculate mean predictive accuracy
    meanPA <-  mean(PA)

    # calculate mean feature importance
    meanFI <-  rowMeans(FI)
    names(meanFI) <- features

    # sort in descending order
    meanFI <- sort(meanFI, decreasing=TRUE)

    processout <- rbind(processout, c(" Report the mean predictive accuracy and mean feature importance."))
    write.table(processout, file=finalfile, row.names=FALSE)
    message(" Report the mean predictive accuracy and mean feature importance.")

    return(list(num_proteins = num_proteins, # number of proteins
                num_samples = num_samples, # number of samples per group
                mean_predictive_accuracy = meanPA, # mean predictive accuracy
                mean_feature_importance = meanFI, # mean feature importance
                results = results # fitted models
    ))

}


