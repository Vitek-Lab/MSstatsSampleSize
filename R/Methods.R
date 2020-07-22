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


#' For each simulated dataset, calculate predictive accuracy on validation set
#'
#' @param index protein abundance data for one protein.
#' @param classifier A character string with the name of the classifier to be used
#' @return A list with (1) predictive accuracy on validation set and
#' (2) the important features
#' @keywords internal
.classification_performance <- function(index, classifier, train_x_list,
                                       train_y_list, valid_x, valid_y, top_K,
                                       tunegrid){
    # record the train x and y
    x <- as.data.frame(train_x_list[[index]])
    y <- as.factor(train_y_list[[index]])

    model <- .classification_model(x = x, y = y, classifier = classifier,
                                  tunegrid = tunegrid)

    pred.features <- .feature_importance(model, classifier, top_K)

    pred.model <- .classification_model(x = x[, pred.features], y = y,
                                       classifier = classifier, tunegrid = tunegrid)

    ## Predict validation data
    pred_y <- predict(pred.model, valid_x[, pred.features])
    ## Calculate predictive accuracy on validation data
    acc <- sum(diag(table(pred_y, valid_y))) / length(pred_y)
    ## record the predictive accuracy and train model
    run <- list(acc = acc, f_imp = pred.features)
    return(run)
}


#' @title Feature Importance
#' @description The method .feature_importance() finds variables importances
#' for the given model and selects the top n features
#' @param model A model object derived from the caret classification models
#' @param classifier The names of the classifier used to train the model on
#' @param top_K an Integer specifying the number of features to consider for 
#' selection
#' @return A string vector of names of the top_K features
#' @keywords internal
#' @import data.table 
.feature_importance <- function(model, classifier, top_K){
    f_imp <- caret::varImp(model, scale = T)
    i_ff <- as.data.table(f_imp$importance, keep.rownames = T)
    
    if(classifier %in% c("svmLinear", "naive_bayes")){
        i_ff[, Overall := rowMeans(i_ff[, -1], na.rm = T)]
    }
    setorder(i_ff, -Overall)
    sel_imp <- i_ff[seq_len(top_K),]$rn
    return(sel_imp)
}


#' Fit a classification model
#'
#' @param x protein abundance data for one protein
#' @param y group information
#' @param classifier A string of classifier 
#' @param tunegrid A data.frame with tuning parameters
#' @return trained classification model
#' @keywords internal
.classification_model <- function(x, y, classifier, tunegrid, ...){
    dots <- list(...)
    MaxNwts <- ifelse(is.null(dots$MaxNwts), 80000, dots$MaxNwts)
    maxit <- ifelse(is.null(dots$maxit), 100, dots$maxit)
    
    if(classifier == "logreg"){
        method <- ifelse(length(unique(y))>2, "multinom", "glm")
    }
    
    if(classifier == "logreg" && method =="glm"){
        ## Train logistic regression on training data
        model <- caret::train(x=x, y=make.names(y), method = method, 
                              trControl = caret::trainControl(method = "none", 
                                                              classProbs = TRUE),
                              maxit = maxit)
    } else {
        if(classifier == 'logreg'){
            classifier <- method
        }
        model <- caret::train(x=x, y=make.names(y),
                              method = classifier,
                              trControl = caret::trainControl(method = "none",
                                                              classProbs = TRUE),
                              tuneGrid = tunegrid, MaxNwts = MaxNwts,
                              maxit = maxit)
    }

    return(model)

}


.tuning_params <- function(classifier, mtry = 2, size = 5, decay = 0.1, C = 1,
                           laplace = 0, usekernel = F, adjust = 1, ...){
    req_classifier <- c('rf','nnet','svmLinear', 'naive_bayes')
    func <- as.list(sys.call())[[1]]
    dots <- list(...)
    for (name in names(dots) ) {
        assign(name, dots[[name]])
    }
    
    if(!classifier %in% req_classifier){
        stop("CALL_",func,"_Incorrect classifier, should be one of (",
             paste(req_classifier, collapse = ", "),")")
    }
    
    switch(classifier, 
           "rf" = {tunegrid=data.frame(mtry=mtry)},
           "nnet"={tunegrid=data.frame(size=size, decay=decay)},
           "svmLinear"={tunegrid=data.frame(C=C)},
           "naive_bayes"={tunegrid=data.frame(laplace=laplace, 
                                              usekernel=usekernel, 
                                              adjust=adjust)})
    return(tunegrid)
}


#' Create log file
#'
#' @return \emph{con} A file connection to write the logs into
#' @return \emph{file} Name of the file to which connection is open
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


#' @title Write Log
#' @param log A string log that has to be written to the log file
#' @param log_type A type of log to be written (purely for informational purpose)
#' @param conn A connection object to the logfile
#' @keywords internal
.log_write <- function(log, log_type, conn){
    sink(conn, type="message")
    message(Sys.time()," : ",log_type,"_",log,"...")
    sink(type="message")
    if(log_type == 'ERROR'){
        message(Sys.time()," : ",log_type,"_",log,"...")
    }
    
}


#' @title Status updates
#' @description Logs all messages comming through to the log file and to the console
#' also provided status updates to the shiny interface with the right options
#' @param detail Message that need to be logged/displayed
#' @keywords internal
.status <- function(detail, ...){
    dots <- list(...)
    dots$func <- ifelse(is.null(dots$func), as.character(as.list(
        sys.call(-1))[[1]]),
        as.character(dots$func))
    mess <- sprintf("CALL_%s_%s", dots$func, detail)
    if(!is.null(dots$log)){
        .log_write(log = mess, log_type = "INFO", conn = dots$log)
    }
    if(!is.null(dots$session) && !is.null(dots$value))
        shiny::setProgress(value = dots$value, message = "Progress:",
                           detail = detail, session = dots$session)
    message(Sys.time()," : ",mess,"...")
}


#' @title Catch Faults
#' @description A wrapper function that catches exceptions and logs them appropriately
#' @param conn A connection object generated by .logGeneration
#' @param session A session object for using with the shiny server
#' @keywords internal
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
    if(isOpen(conn$con) && is.null(session)){
        close(conn$con)
    }
    return(res)
}


.find_log <- function(...){
    dots <- list(...)
    session <- dots$session
    if(is.null(dots$log_conn)){
        conn = mget("LOG_FILE", envir = .GlobalEnv,
                    ifnotfound = NA)
        if(is.na(conn)){
            rm(conn)
            conn <- .logGeneration()
        } else{
            conn <- .logGeneration(file = conn$LOG_FILE)   
        }
    }else{
        conn <- dots$log_conn
    }
    return(conn)
}


.check_missing_values <- function(x){
    return(apply(x, 2, function(x){
        any(is.na(x) | is.infinite(x) | is.nan(x) | is.null(x) == T)
    }))
}


#' @title Data Checking
#' @description This function check for data consistency in the input dataset
#' and provided either a failure or success message when the data checks are
#' passed successfully
#' @param data A abundance file input 
#' @param annotation Annotation file for the abundance
#' @param conn A connection 
#' @keywords internal
.data_checks <- function(data, annotation, conn){
    packageStartupMessage(Sys.time()," : Checking Data for consistency...", 
                          appendLF = F)
    func <- as.list(sys.call())[[1]]
    
    #Check for any duplicated column names in the data provided
    dups <- any(duplicated(colnames(data)) == T)
    if(dups){
        packageStartupMessage(" Failure")
        stop("CALL_",func,
             "__Please check the column names of 'data'. There are duplicated 'Run'.")
    }
    
    #Check if input data is consistent with the required format
    required.annotation <- c('Condition', 'BioReplicate', 'Run')
    consistent_cols <- !setequal(required.annotation, colnames(annotation))
    if(consistent_cols){
        packageStartupMessage(" Failure")
        nf <- required.annotation[!colnames(annotation) %in% required.annotation]
        stop("CALL_", func,"_",nf,
             " is not provided in Annotation, please check annotation file")
    }

    ic <- setequal(annotation$BioReplicate, colnames(data)) && 
        nrow(annotation) == ncol(data)
    if(!ic){
        packageStartupMessage(" Failure")
        stop("CALL_",func,"_",
             "Please check the annotation file. 'Bioreplicate' must match with",
             "the column names of 'data'.")
    }

    
    if(length(unique(annotation$Condition)) < 2){
        stop("Need at least two conditions to do simulations")
    }
    
    missing_values <- .check_missing_values(annotation)
    if(any(missing_values) == T){
        packageStartupMessage(" Failure")
        stop("CALL_",func,
             "_'NA|NAN|NULL' not permitted in 'Run', 'BioReplicate' or",
             "'Condition' of annotaion.")
    }
    
    data <- data[, annotation$Run]
    group <- as.factor(as.character(annotation$Condition))
    
    if (nrow(data) == 0) {
        stop("CALL_",func,
             "_Please check the column names of data and Run in annotation.")
    }
    
    temp <- unique(annotation[annotation$Run %in% colnames(data),
                              c("BioReplicate", "Condition")])
    temp$Condition <- factor(temp$Condition)
    temp$BioReplicate <- factor(temp$BioReplicate)
    if(any(table(temp$Condition) < 3)){
        stop("CALL_",func,
             "_Each condition must have at least three biological replicates.")
    }
    packageStartupMessage(" Success")
    
    return(list("data" = data, "group" = group))
}


.log2_trans <- function(trans, data, conn, ...){
    func <- as.list(sys.call())[[1]]
    if(is.logical(trans)){
        if(trans){
            data <- log2(data)
        }
        .status(detail = "Negative values if any where replaced with NA",
                log = conn$con, func = func, ...)
        data[!is.na(data) & data < 0] <- NA
        return(data)
    }else{
        stop("CALL_", func,
             "_log2Trans should be logical. Please provide either TRUE or FALSE")
    }
}


#' @title Plot PCA outputs
#' @description A utility wrapper for plotting pca data as returned by the `do_prcomp()`
#' @param data A data frame containing the Principal components to be plotted
#' @param exp_var A vector containing numeric values of the expected variance
.pca_plot <- function(data, exp_var, ...){
    dots <- list(...)
    title <- ifelse(is.null(dots$title),"Input dataset", dots$title)
                    
    p <- ggplot(data = data,
                aes(x = PC1, y = PC2, color = group)) +
        geom_point(size =  ifelse(is.null(dots$dot_size), 3, dots$dot_size)) + 
        stat_ellipse()+
        labs(title = title,
             x = sprintf("PC1 (%s%% explained var.)", exp_var[1]),
             y = sprintf("PC2 (%s%% explained var.)", exp_var[2])) +
        theme_MSstats(...)
    return(p)
}


#' @title Do Principal Component Analysis
#' @description A wrapper function to the base `prcomp()` formats the results 
#' in a required output format
#' @param sim_x A dataframe of the feature variables
#' @param sim_y A dataframe of the predictor vairable
#' @return A named list of the outputs 
.do_prcomp <- function(sim_x, sim_y){
    result.pca <- prcomp(sim_x, scale. = TRUE)
    summary.pca <- summary(result.pca)
    important.pc <- result.pca$x[, 1:2]
    pc.result <- data.frame(important.pc, group = sim_y)
    exp.var <- summary.pca$importance
    exp.var <- format(exp.var[2, 1:2] * 100, digits = 3)
    return(list("pc.result" = pc.result, "exp.var" = exp.var))
}


#' @title MSstats Theme
#' @description A utility function that standardized all the ggplot themes
#' @param x.axis.size A numeric value for size for the elements on the x-axis
#' @param y.axis.size A numeric value for size for the elements on the y-axis
#' @param legend.size A numeric value fot the size of the legend
theme_MSstats <- function(x.axis.size = 10, y.axis.size = 10, 
                          legend.size = 11, margin = 0.5, leg.dir="horizontal",
                          download = F,...){
    
    dots <- list(...)
    x.axis.size <- ifelse(is.null(dots$x.axis.size), x.axis.size,
                          dots$x.axis.size)
    
    y.axis.size <- ifelse(is.null(dots$y.axis.size), y.axis.size,
                          dots$y.axis.size)
    
    legend.size <- ifelse(is.null(dots$legend.size), legend.size,
                          dots$legend.size)
    leg.dir <- ifelse(is.null(dots$leg.dir), leg.dir,
                      dots$leg.dir)
    download <- ifelse(is.null(dots$download), download,
                       dots$download)
    leg.pos <- ifelse(is.null(dots$leg.pos), "top", dots$leg.pos)
    
    th <- ggplot2::theme(panel.background = element_rect(fill = "white", 
                                                         colour = "black"),
                         panel.grid.major = element_line(colour = "gray95"), 
                         panel.grid.minor = element_blank(), 
                         strip.background = element_rect(fill = "gray95"), 
                         strip.text.x = element_text(colour = c("#00B0F6"), 
                                                     size = 14),
                         axis.text.x = element_text(size = x.axis.size, 
                                                    colour = "black"), 
                         axis.text.y = element_text(size = y.axis.size, 
                                                    colour = "black"),
                         axis.ticks = element_line(colour = "black"), 
                         axis.title.x = element_text(size = x.axis.size + 5,
                                                     vjust = -0.4), 
                         axis.title.y = element_text(size = y.axis.size + 5,
                                                     vjust = 0.3),
                         title = element_text(size = x.axis.size + 4,
                                              vjust = 1.5),
                         legend.key = element_rect(fill = "white", 
                                                   colour = "white"),
                         legend.direction = leg.dir,
                         legend.position = leg.pos, 
                         legend.text = element_text(size = legend.size), 
                         legend.title = element_blank(),
                         plot.margin = unit(rep(margin,4), "cm"))
    
    if(!download)
        th <- th + ggplot2::theme(legend.position=c(1, 1.05),
                                  legend.justification="right")
    
    return(th)
}


#' @title Identify optimal sample size
#' @param df A dataframe containing the sample size, accuracies
#' @param cutoff A numeric values to determine the threshold for identifying the
#' sample size
#' @param use_h2o A logical input to use logic for h2o data
#' @return \emph{df} A data.frame
#' @return \emph{opt} A integer with the optimal sample size
#' @return \emph{y_lim} A vector of limits for a ggplot2 object
#' @keywords internal
#' @import data.table
.identify_optimal <- function(df, cutoff, use_h2o = F){
    
    if(!is.logical(use_h2o)){
        stop("'use_h2o' needs to be logical")
    }
    
    df <- as.data.table(df)
    df[, mean_acc := mean(acc), sample]
    opt_df <- unique(df[, .(sample, mean_acc)])
    setorder(opt_df, -sample)
    if(use_h2o){
        setorder(opt_df, sample)
    }
    
    opt_df[, sample := as.numeric(as.character(sample))]
    opt_df[, lg := (mean_acc - shift(mean_acc))/(sample - shift(sample))]
    opt_df[, optimal := ifelse(lg >= cutoff, T, F)]
    if(nrow(opt_df[, .N, optimal][optimal == T]) != 0){
        opt_val <- opt_df[optimal == T][which.min(lg), sample]
    } else {
        opt_val <- opt_df[which.min(sample), sample]
    }
    
    y_lim <- c(df[,min(mean_acc, na.rm = T)]-0.1, 1)
    df[sample == opt_val, fill_col := "red"]
    
    return(list('df' = df, 'opt' = opt_val, 'y_lim' = y_lim))
}


#' Plot Accuracy 
#' @keywords internal
#' @import ggplot2
#' @importFrom scales pretty_breaks
.plot_acc <- function(df, y_lim, optimal_ss, ...){
    g <- ggplot(data = df, aes(x = sample))+
        geom_boxplot(aes(y = acc, group = sample, fill = fill_col), alpha = 0.5)+
        scale_fill_identity()+
        geom_point(aes(y = mean_acc))+
        geom_line(aes(y = mean_acc, group = 1), size = 0.75, color = "blue")+
        labs(x = "Simulated Sample Size Per Group", y = "Predictive Accuracy",
             #title = sprintf("Classifier %s", alg),
             subtitle = sprintf("Optimum accuracy achieved when sample size per group is : %s",
                                optimal_ss))+
        scale_y_continuous(breaks = scales::pretty_breaks(), limit = y_lim)+
        scale_x_continuous(breaks = unique(df$sample))+
        theme_MSstats(...)+
        theme(plot.subtitle = element_text(face = "italic", color = "red"))
    
    return(g)
}


#' Plot Variable Importance
#' @keywords internal
#' @import ggplot2
.plot_imp <- function(df, samp = NA, ylim, facet = F,...){
    if(!is.na(samp)){
        df <- subset(df, sample == samp)
    }
    g <- ggplot(data = df, aes(x = reorder(protein, frequency), 
                               y = frequency))+
        geom_col()+
        labs(x = "Protein", y = "Frequency")+
        scale_y_continuous(breaks = scales::pretty_breaks(),
                           limit = ylim)+
        theme_MSstats(...)+
        coord_flip()
    
    if(facet){
        titles <- sprintf("%s Samples/Group", unique(df$sample))
        names(titles) <- unique(df$sample)
        titles_lb <- as_labeller(titles)
        g <- g+
            facet_wrap(sample~., ncol = 2, scales = 'free', labeller = titles_lb)
    }else{
        g <- g+
            labs(title = sprintf("%s Samples/Group", samp))
    }
    return(g)
}


#' @title Plot QC boxplots
#' @description Plots an interactive plotly graphic with boxplots for all the 
#' proteins and their respective conditions with abundance indicated in the log
#' scale on the y axis
#' @param data A formatted data.frame with abundance, bioreplicate and conditions 
#' information
#' @return A plotly object
#' @keywords internal
#' @import data.table
qc_boxplot <- function(data = NULL ,annot = NULL){
    #create the interactive boxplot for all the different proteins found in the data
    data <- as.data.table(data, keep.rownames = T)
    setnames(data, 'rn', 'proteins')
    data <- melt(data, id.vars = 'proteins', variable.name = "BioReplicate",
                 value.name = "Abundance")
    annot <- as.data.table(annot)
    data <- merge(data, annot, by = "BioReplicate") 
    
    box_plot <- plotly::plot_ly(data = data[!is.na(Abundance)],
                                y = ~Abundance, x = ~BioReplicate, color = ~Condition,
                                type = "box") %>%
        plotly::layout(xaxis = list(title="Biological Replicate",showticklabels = TRUE,
                                    tickangle = -45 ), 
                       yaxis = list(title="Protein abundance"),
                       legend = list(orientation = "h", #position and of the legend
                                     xanchor = "center",
                                     x = 0.5, y = 1.1)) %>%
        plotly::config(displayModeBar = F) #hide controls of the plotly chart
    
    box_plot
}


.format_df <- function(dat, sample, top_n = NA){
    df <- stack(dat)
    df$sample <- sample
    if(is.na(top_n)){
        return(df)
    }else{
        return(df[seq_len(top_n),])   
    }
}


#' @title Get summary table for the annotation data
#' @description Get the summary for the unique bioreplicates and number of
#' MS runs for the data that is provided
#' @param data A data.frame with the annotation data containing the Bioreplicates
#' Runs and condition information
#' @return A data.frame with the counts of bioreplicates for each conditions
#' and the number of runs as well
.format_summary_table <- function(data = NULL){
    #create crosstable for the conditions vs bioreplicates
    data <- as.data.table(data)
    biorep <- unique(data[,.(BioReplicate, Condition)])
    biorep <- xtabs(~Condition, data = biorep)
    
    #create crosstable for the conditions vs runs if runs data exists
    if(any(c("run","Run") %in% names(data))){
        msruns <- unique(data[,.(Run, Condition)])
        msruns <- xtabs(~Condition, data = msruns)
    }else{
        msruns <- rep(0, length(names(biorep)))
        names(msruns) <- names(biorep) #make runs data 0 if not found
    }
    #format it correctly
    #summary <- rbind(biorep, msruns)
    summary <-matrix(biorep, nrow = 1)
    colnames(summary) <- names(biorep)
    summary <- summary[,which(colSums(summary, na.rm = T) > 0)]
    sum_table <- matrix(summary, nrow=1)
    #rownames(summary) <- c("# of Biological Replicates", "# of MS runs")
    dimnames(sum_table) <- list("# of Biological Replicates", names(summary))
    return(sum_table)
}


#' @title 
#' @description A wrapper function to classify multiple simulated datasets
#' @param n_samp A character vector consisting the of the different samples simulated
#' @param sim_data A list object, output out the `simulate_grid` function
#' @param classifier A string specifying the classifier to use
#' @param k A integer value for the number of features to select
#' @param session A session object for the shiny app
#' @keywords internal
ss_classify_caret <- function(n_samp, sim_data, classifier, k = 10,
                              session = NULL, ...){
    
    dots <- list(...)
    session <- dots$session
    conn <- .find_log(...)
    samp <- unlist(strsplit(n_samp,","))
    res <-  list()
    max_val <- length(sim_data)
    iter <- 0
    
    for(i in seq_along(samp)){
        if(!is.null(session)){
            #Provides progress bars to the shiny UI
            iter = iter + 1/max_val
        }
        .status(detail = sprintf("Classifying Sample %s of %s", i, length(samp)),
                log = conn$con, session = session, value = iter)
        res[[i]] <- designSampleSizeClassification(sim_data[[i]],
                                                   classifier = classifier, 
                                                   top_K = k, session = session,
                                                   ...)
    }
    return(res)
}


#' @title Simulate datasets to be tested out
#' @description A wrapper function for the `simulateDataset` function from the 
#' MSstatsSampleSize package which enables simulating datasets for running experiments
#' @param data 
#' @param annot
#' @param num_simulations
#' @param exp_fc
#' @param list_diff_proteins
#' @param sel_simulated_proteins
#' @param prot_proportion
#' @param prot_number
#' @param samples_per_group
#' @param sim_valid
#' @param valid_samples_per_grp
#' @param seed
#' @param session
#' @return 
simulate_grid <- function(data = NULL, annot = NULL, num_simulation, exp_fc,
                          list_diff_proteins, samples_per_group, sim_valid,
                          valid_samples_per_grp, seed, est_var, conn, 
                          session = NULL){
    #check if seed value required
    .status(detail = "Setting Up Data Simulation Runs", value = 0.1,
            session = session, log = conn$con)
    if(seed != -1)
        set.seed(seed)
    
    if(exp_fc != "data"){
        .status(detail = "Extracting Fold Change Informations", value = 0.15, 
               session = session, log = conn$con)
        .status(detail = sprintf("List of differential proteins selected: (%s)", 
                                list_diff_proteins), session = session, 
                log = conn$con)
        
        diff_prots <- unlist(strsplit(list_diff_proteins, ","))
        fc <- exp_fc$`Fold Change Value`
        names(fc) <- exp_fc$orig_group
    } else{
        diff_prots <- NULL
        fc <- exp_fc
    }
    .status(detail = "Extracting Number of Samples Information", value = 0.2,
            session = session, log = conn$con)
    .status(detail = sprintf("Number of samples per group: (%s)",
                             samples_per_group),
            session = session , log = conn$con)
    
    samp <- as.numeric(unlist(strsplit(samples_per_group, ",")))
    shiny::validate(shiny::need(all(!is.na(samp)),
                                sprintf("Samples Per Group need to be numeric values, Found : %s",
                                        samples_per_group)),
                    shiny::need(all(samp >= 1), "All samples Need to be >= 1"))
    
    if(sim_valid){
        .status(detail = "Validation Simulation requested", value = 0.2,
                session = session, log = conn$con)
    }
    
    .status(detail = "Starting Simulation", value = 0.3, session = session,
            log = conn$con)
    
    sim <- list()
    for(i in samp){
        .status(detail = sprintf("Running Simulation for sample %s of %s",
                                 which(i == samp), length(samp)),
               value = which(i==samp)/length(samp), 
               session = session, log = conn$con)
        
        sim[[paste(i)]] <- simulateDataset(data = data,
                                           annotation = annot,
                                           num_simulations = num_simulation,
                                           expected_FC = fc,
                                           list_diff_proteins =  diff_prots,
                                           samples_per_group = i,
                                           simulate_valid = as.logical(sim_valid),
                                           valid_samples_per_group = valid_samples_per_grp,
                                           session = session)
    }
    .status(detail = "Simulation Complete", value = 0.9, session = session,
            log = conn$con)
    return(sim)
}

