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
                        annotation) {

    ###############################################################################
    ## log file
    ## save process output in each step
    allfiles <- list.files()

    num <- 0
    filenaming <- "MSstatsSampleSize-ProgressReport"
    finalfile <- "MSstatsSampleSize-ProgressReport.log"

    while(is.element(finalfile, allfiles)) {
        num <- num + 1
        finalfile <- paste0(paste(filenaming, num, sep="-"), ".log")
    }

    session <- sessionInfo()
    sink("sessionInfo.txt")
    print(session)
    sink()

    processout <- as.matrix(read.table("sessionInfo.txt", header=TRUE, sep="\t"))
    write.table(processout, file=finalfile, row.names=FALSE)

    processout <- rbind(processout, as.matrix(c(" "," ","MSstatsSampleSize - estimateVar function"," "), ncol=1))

    ###############################################################################
    ## input checking
    if(anyDuplicated(colnames(data)) != 0){
        processout <- rbind(processout, c("ERROR: Please check the column names of 'data'.
                                          There are duplicated 'BioReplicate'."))
        write.table(processout, file=finalfile, row.names=FALSE)

        stop("Please check the column names of 'data'.
             There are duplicated 'BioReplicate'.  \n")

    }

    if (!all(annotation$BioReplicate %in% colnames(data)) &
        nrow(annotation) == ncol(data)) {
        processout <- rbind(processout, c("ERROR: Please check the annotation file.
                                          'BioReplicate' must match with the column names of 'data'."))
        write.table(processout, file=finalfile, row.names=FALSE)

        stop("Please check the annotation file.
             'BioReplicate' must match with the column names of 'data'.  \n")

    }

    if (length(unique(annotation$Condition)) < 2){
        processout <- rbind(processout, c("ERROR: Need at least two conditions to do simulation."))
        write.table(processout, file=finalfile, row.names=FALSE)

        stop("Need at least two conditions to do simulation.  \n")

    }

    if (any(is.na(annotation$BioReplicate)) |
        any(is.na(annotation$BioReplicate))){
        processout <- rbind(processout, c("ERROR: NA not permitted in 'BioReplicate' or 'Condition' of annotaion."))
        write.table(processout, file=finalfile, row.names=FALSE)

        stop("NA not permitted in 'BioReplicate' or 'Condition' of annotaion.  \n")
    }

    ## match between data and annotation
    data <- data[, annotation$BioReplicate]
    group <- as.factor(as.character(annotation$Condition))

    if (nrow(data) == 0) {
        processout <- rbind(processout, c("ERROR: Please check the column names of data and BioReplicate in annotation."))
        write.table(processout, file=finalfile, row.names=FALSE)

        stop("Please check the column names of data and BioReplicate in annotation. \n")
    }

    processout <- rbind(processout, c(paste0("Summary : number of proteins in the input data = ", nrow(data) )))
    processout <- rbind(processout, c(paste0("Summary : number of samples in the input data = ", ncol(data) )))
    write.table(processout, file=finalfile, row.names=FALSE)

    ###############################################################################
    processout <- rbind(processout, c("Preparing variance analysis..."))
    write.table(processout, file=finalfile, row.names=FALSE)
    message(" Preparing variance analysis...")

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

    processout <- rbind(processout, c(" Variance analysis completed."))
    write.table(processout, file=finalfile, row.names=FALSE)
    message(" Variance analysis completed.")

    return(list(model = Models,
                protein = Proteins,
                promean = SampleMean,
                mu = GroupMean,
                sigma = GroupVar))
}

