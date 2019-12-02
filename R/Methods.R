#' Simulate extended datasets for sample size estimation
#'
#' @param m number of samples per group to simulate
#' @param mu a matrix of mean abundance in each phenotype group and each protein
#' @param sigma a matrix of variance in each phenotype group and each protein
#' @return \emph{X} A simulated matrix with required sample size
#' @return \emph{Y} Group information corresponding with \emph{X}
#' @keywords internal
.sampleSimulation <- function(m, mu, sigma) {

    # record number of proteins and groups
    nproteins <- nrow(mu)
    ngroup <- ncol(mu)
    # make sure the column names are same bewteen mu and sigma
    sigma <- sigma[, colnames(mu)]
    ## Determine the size of each phenotype group
    samplesize <- rep(m, ngroup)

    ## Simulate the data matrix
    sim_matrix <- matrix(rep(0, nproteins * sum(samplesize)), ncol=sum(samplesize))
    for (i in seq_len(nproteins)) {
        index <- 1
        for (j in seq_len(ngroup)) {
            sim_matrix[i, index:(index+samplesize[j]-1)] <- rnorm(samplesize[j], mu[i, j], sigma[i, j])
            index <- index + samplesize[j]

        }
    }

    sim_matrix <- t(sim_matrix)
    colnames(sim_matrix) <- rownames(mu)
    #Simulate the phenotype information
    group <- rep(colnames(mu), times=samplesize)

    index <- sample(length(group), length(group))
    sim_matrix <- sim_matrix[index, ]
    group <- group[index]

    return(list(X=sim_matrix,
                Y=as.factor(group)))

}



#' For each protein, impute the missing values based on the observed values
#'
#' @param data protein abundance data for one protein.
#' @return Imputed protein abundance data
#' @keywords internal
.randomImputation <- function (data){

    missing <- is.na(data) # count missing values
    n.missing <- sum(missing)
    data.obs <- data[!missing] # keep the observed values
    imputed <- data
    # impute missing values by randomly selecting observed values
    imputed[missing] <- sample(data.obs, n.missing, replace=TRUE)
    return (imputed)

}


#' Create log file
#'
#' @return \emph{finalfile} The log file name
#' @return \emph{processout} The current log information
#' @keywords internal
.logGeneration <- function (){

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

    return(list(finalfile = finalfile,
                processout = processout))

}


#' For each simulated dataset, calculate predictive accuracy on validation set
#'
#' @param index protein abundance data for one protein.
#' @return A list with (1) predictive accuracy on validation set and
#' (2) the trained classification model
#' @keywords internal
.classificationPerformance <- function(index, classifier, train_x_list, train_y_list, valid_x, valid_y, top_K){

    # record the train x and y
    x <- as.data.frame(train_x_list[[index]])
    y <- as.factor(train_y_list[[index]])

    model <- .classificationModel(x = x,
                                  y = y,
                                  classifier = classifier)

    pred.features <- rownames(varImp(model, scale = TRUE)$importance)[seq(top_K)]

    pred.model <- .classificationModel(x = x[, pred.features],
                                       y = y,
                                       classifier = classifier)

    ## Predict validation data
    pred_y <- predict(pred.model, valid_x[, pred.features])

    ## Calculate predictive accuracy on validation data
    acc <- sum(diag(table(pred_y, valid_y))) / length(pred_y)
    acc

    # ## Predict validation data with all the proteins
    # pred_y <- predict(model, valid_x)
    #
    # ## Calculate predictive accuracy on validation data
    # acc <- sum(diag(table(pred_y, valid_y))) / length(pred_y)
    # acc

    ## record the predictive accuracy and train model
    run <- list(acc = acc, model = model)

    return(run)

}


#' Fit a classification model
#'
#' @param x protein abundance data for one protein
#' @param y group information
#' @return trained classification model
#' @keywords internal
.classificationModel <- function(x, y, classifier){

    if(classifier == "logreg"){
        ## Train logistic regression on training data
        model <- caret::train(x=x, y=make.names(y),
                              trControl = caret::trainControl(method = "none", classProbs = TRUE),
                              # trControl = trainControl(method = "cv", number=10, classProbs = TRUE),
                              method = "glm",
                              family = "binomial",)

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

    return(model)

}

