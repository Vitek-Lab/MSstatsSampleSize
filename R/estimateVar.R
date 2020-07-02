#' Estimate the mean abundance and variance of each protein in each condition.
#' @details The function fits intensity-based linear model on the input data `data'.
#' This function outputs variance components and mean abundance for each protein.
#'
#' @param data Data matrix with protein abundance.
#' Rows are proteins and columns are Biological replicates or samples.
#' @param annotation Group information for samples in data.
#' `Run' for MS run, `BioReplicate' for biological subject ID and `Condition' for group information are required.
#' `Run' information should be the same with the column of `data'. Multiple `Run' may come from same `BioReplicate'.
#' @return \emph{model} is the list of linear models trained for each protein.
#' @return \emph{mu} is the mean abundance matrix of each protein in each phenotype group.
#' @return \emph{sigma} is the sd matrix of each protein in each phenotype group.
#' @return \emph{promean} is the mean abundance vector of each protein across all the samples.
#' @return \emph{protein} is proteins, correpsonding to the rows in \emph{mu} and \emph{sigma} or the element of \emph{promean}.
#'
#' @author Ting Huang, Meena Choi, Olga Vitek
#'
#' @examples
#' data(OV_SRM_train)
#' data(OV_SRM_train_annotation)
#'
#' # estimate the mean protein abunadnce and variance in each condition
#' variance_estimation <- estimateVar(data = OV_SRM_train,
#'                                    annotation = OV_SRM_train_annotation)
#'
#' # the mean protein abundance in each condition
#' head(variance_estimation$mu)
#'
#' # the standard deviation in each condition
#' head(variance_estimation$sigma)
#'
#' # the mean protein abundance across all the conditions
#' head(variance_estimation$promean)
#'
#' @export
#' @importFrom utils sessionInfo read.table write.table
#' @importFrom stats lm coef anova
#'
estimateVar <- function(data,
                        annotation, ...) {

    ###############################################################################
    ## log file
    ## save process output in each step
    dots <- list(...)
    if(is.null(dots$log_conn)){
        conn <- .logGeneration()  
    }else{
        conn <- dots$log_conn
    }
    func <- as.list(sys.call())[[1]]
    .status("Estimating Variance", log = conn$con, func = func, ...)
    
    ###############################################################################
    ## input checking
    if(anyDuplicated(colnames(data)) != 0){
        .status("ERROR: Please check the column names of 'data'. 
                There are duplicated 'Run'", log = conn$con, func = func, ...)
        stop("Please check the column names of 'data'.
             There are duplicated 'Run'.")
    }

    ## check whether annotation has requried columns
    required.annotation <- c('Condition', 'BioReplicate', 'Run')

    if (!all(required.annotation %in% colnames(annotation))) {

        missedAnnotation <- which(!(required.annotation %in% colnames(annotation)))
        .status(sprintf("%s is not provided in Annotation. Please check the annotation file.",
                        toString(required.annotation[missedAnnotation])),
                log = conn$con, func = func, ...)

        stop(paste(toString(required.annotation[missedAnnotation]),
                   "is not provided in Annotation. Please check the annotation file."))
    }

    if (!all(annotation$Run %in% colnames(data)) &
        nrow(annotation) == ncol(data)) {
        .status("ERROR: Please check the annotation file.
                'Run' must match with the column names of 'data'.",
                log = conn$con,func = func, ...)
        stop("Please check the annotation file.
             'Run' must match with the column names of 'data'.")
    }

    if (length(unique(annotation$Condition)) < 2){
        .status("ERROR: Need at least two conditions to do simulation.",
                log = conn$con, func = func, ...)
        stop("Need at least two conditions to do simulation.")
    }

    if (any(is.na(annotation$Run)) |
        any(is.na(annotation$BioReplicate))|
        any(is.na(annotation$Condition))){
        .status("ERROR: NA not permitted in 'Run', 'BioReplicate' or 'Condition' of annotaion.",
                log = conn$con, func = func, ...)
        stop("NA not permitted in 'Run', 'BioReplicate' or 'Condition' of annotaion.")
    }

    ## match between data and annotation
    data <- data[, annotation$Run]
    group <- as.factor(as.character(annotation$Condition))

    if (nrow(data) == 0) {
        .status("ERROR: Please check the column names of data and Run in annotation.",
                log = conn$con, func = func, ...)
        stop("Please check the column names of data and Run in annotation. \n")
    }

    .status(sprintf("Summary : number of proteins in the input data = %s",
                    nrow(data)), log = conn$con, func = func)
    .status(sprintf("Summary : number of samples in the input data = %s",
                    ncol(data)), log = conn$con, func = func)
    
    ###############################################################################
    .status("Preparing variance analysis", log = conn$con, 
            func = func)
    ## unique groups
    groups <- as.character(unique(group))
    ngroup <- length(groups)
    nproteins <- nrow(data)

    GroupVar <- matrix(rep(NA, nproteins * ngroup), ncol = ngroup)
    GroupMean <- matrix(rep(NA, nproteins * ngroup), ncol = ngroup)
    SampleMean <- NULL # mean across all the samples
    Proteins <- NULL # record the proteins to simulate
    Models <- list()
    count = 0
    for (i in seq_len(nrow(data))) {
        sub<- data.frame(ABUNDANCE = unname(unlist(data[i, ])), GROUP = factor(group), row.names = NULL)
        sub <- sub[!is.na(sub$ABUNDANCE), ]

        ## train the one-way anova model
        df.full <- suppressMessages(try(lm(ABUNDANCE ~ GROUP , data = sub), TRUE))
        if(!inherits(df.full, "try-error")){
            abun <- coef(df.full)
            if(length(abun) == ngroup){

                count <- count + 1
                # total variance in protein abundance
                var <- anova(df.full)['Residuals', 'Mean Sq']
                # estimate mean abundance of each group
                abun[-1] <- abun[1] + abun[-1]
                names(abun) <- gsub("GROUP", "", names(abun))
                names(abun)[1] <- setdiff(as.character(groups), names(abun))
                abun <- abun[groups]

                # save group mean, group variance
                Models[[rownames(data)[i]]] <- df.full
                GroupVar[count, ] <- rep(sqrt(var), times=length(abun))
                GroupMean[count, ] <- abun
                SampleMean <- c(SampleMean, mean(sub$ABUNDANCE, na.rm = TRUE))
                Proteins <- c(Proteins, rownames(data)[i])
            }
        }
    }

    # only keep the rows with results
    GroupMean <- GroupMean[1:count, ]
    GroupMean <- GroupMean[1:count, ]
    # assign the row and column names
    rownames(GroupVar) <- Proteins
    colnames(GroupVar) <- groups
    rownames(GroupMean) <- Proteins
    colnames(GroupMean) <- groups
    names(SampleMean) <- Proteins

    .status("Variance analysis completed.", log = conn$con, func = func)
    close(conn$con)
    
    return(list(model = Models,
                protein = Proteins,
                promean = SampleMean,
                mu = GroupMean,
                sigma = GroupVar))
}

