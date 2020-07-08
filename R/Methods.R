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
.randomImputation <- function(data){

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
.logGeneration <- function(...){
    dots <- list(...)
    if(is.null(dots$file)){
        LOG_DIR <- file.path(getwd(),'logs')
        dir.create(LOG_DIR, showWarnings = F)
        FILE <- sprintf("MSstatsSSE_Log_%s.Rout", 
                        format(Sys.time(),"%Y%m%d%H%M%S"))
        assign("LOG_FILE", file.path(LOG_DIR, FILE))
        assign("FILE_CONN", file(LOG_FILE, open='a'))
        writeLines(capture.output(sessionInfo()), FILE_CONN)
        writeLines("\n\n ############## LOG ############# \n\n", FILE_CONN)
    } else{
        FILE_CONN <- file(dots$file, open = 'a')
        LOG_FILE <- dots$file
    }
    return(list(con = FILE_CONN, file = LOG_FILE))
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



#' Fit a classification model
#'
#' @param detail Details of the log entry
#' @param value value ranging from 0 - 1 for progress bars used only with shiny
#' @param session A session object for shiny operations
#' @return Console logging
#' @keywords internal
.status <- function(detail, ...){
    dots <- list(...)
    dots$func <- ifelse(is.null(dots$func),"'__'", as.character(dots$func))
    mess <- sprintf("CALL_%s_%s", dots$func, detail)
    if(!is.null(dots$log)){
        .log_write(log = mess, log_type = "INFO", conn = dots$log)
    }
    if(!is.null(dots$session) && !is.null(dots$value))
        shiny::setProgress(value = dots$value, message = "Progress:",
                           detail = detail, session = dots$session)
    message(Sys.time()," : ",mess,"...")
}


.catch_faults <- function(..., conn, session=NULL){
    warn <- err <- NULL
    res <- withCallingHandlers(tryCatch(..., error = function(e) {
        assign("LOG_FILE", conn$file, envir = .GlobalEnv)
        err <<- conditionMessage(e)
        .log_write(log = err, log_type = "ERROR", conn = conn$con)
        if(!is.null(session)){
            shiny::showNotification(as.character(err), duration = 20, 
                                    type = "error",session = session, 
                                    id = "error") 
            shiny::validate(shiny::need(is.null(err), as.character(err)))
        }
        NULL
    }), warning = function(w) {
        warn <<- append(warn, conditionMessage(w))
        .log_write(log = warn, log_type = "WARNING", conn = conn$con)
        invokeRestart("muffleWarning")
    })
    assign("LOG_FILE", conn$file, envir = .GlobalEnv)
    close(conn$con)
    return(res)
}


.log_write <- function(log, log_type, conn){
    sink(conn, type="message")
    message(Sys.time()," : ",log_type,"_",log,"...")
    sink(type="message")
    if(log_type == 'ERROR'){
        close(conn)
    }
}

.check_missing_values <- function(x){
    return(apply(x, 2, function(x){
        any(is.na(x) | is.infinite(x) | is.nan(x) | is.null(x) == T)
    }))
}

.data_checks <- function(data, annotation){
    packageStartupMessage(Sys.time()," : Checking Data for consistency...", 
                          appendLF = F)
    func <- as.list(sys.call())[[1]]
    dups <- any(duplicated(colnames(data)) == T)
    if(dups){
        packageStartupMessage(" Failure")
        stop("CALL_",func,
             "__Please check the column names of 'data'. There are duplicated 'Run'.")
    }
    required.annotation <- c('Condition', 'BioReplicate', 'Run')
    consistent_cols <- !setequal(required.annotation, colnames(annotation))
    if(consistent_cols){
        packageStartupMessage(" Failure")
        nf <- required.annotation[!colnames(annotation) %in% required.annotation]
        stop("CALL_", func,"_",nf,
             " is not provided in Annotation, please check annotation file")
    }
    ic <- setequal(annotation$Run, colnames(data)) && nrow(annotation) == ncol(data)
    if(!ic){
        packageStartupMessage(" Failure")
        stop("CALL_",func,"_",
             "Please check the annotation file. 'Run' must match with the column names of 'data'.") 
    }
    missing_values <- .check_missing_values(annotation)
    if(any(missing_values) == T){
        packageStartupMessage(" Failure")
        stop("CALL_",func,
             "_'NA' not permitted in 'Run', 'BioReplicate' or 'Condition' of annotaion.")   
    }
    packageStartupMessage(" Success")
}
