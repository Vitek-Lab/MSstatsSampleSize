#' Sample size calculation plot for hypothesis testing
#'
#' Calculate sample size for future experiments based on intensity-based linear model.
#'
#' @details The function fits intensity-based linear model on the input `data'.
#' Then it uses the fitted models and the fold changes estimated from the models to calculate
#' sample size for hypothesis testing through `designSampleSize' function from MSstats package.
#' It outputs the minimal number of biological replciates per condition to acquire the expected
#' FDR and power under different fold changes.
#'
#' @param data Protein abundance data matrix.
#' Rows are proteins and columns are biological replicates (samples).
#' @param annotation Group information for samples in data.
#' `BioReplicate' for sample ID and `Condition' for group information are required.
#' `BioReplicate' information should match with column names of `data'.
#' @param desired_FC the range of a desired fold change.
#' The first option (Default) is "data",
#' indicating the range of the desired fold change is directly estimated from the input `data',
#' which are the minimal fold change and the maximal fold change in the input `data'.
#' The second option is a vector which includes
#' the lower and upper values of the desired fold change (For example, c(1.25,1.75)).
#' @param select_testing_proteins the standard to select the proteins for hypothesis testing and sample size calculation.
#' The variance (and the range of desired fold change if desiredFC = "data") for sample size calculation
#' will be estimated from the selected proteins.
#' It can be 1) "proportion" of total number of proteins in the input data or
#' 2) "number" to specify the number of proteins.
#' "proportion" indicates that user should provide the value for `protein_proportion' option.
#' "number" indicates that user should provide the value for `protein_number' option.
#' @param protein_proportion Proportion of total number of proteins in the input data to test.
#' For example, input data has 1,000 proteins and user selects `protein_proportion=0.1'.
#' Proteins are ranked in decreasing order based on their mean abundance across all the samples.
#' Then, 1,000 * 0.1 = 100 proteins will be selected from the top list to test.
#' Default is 1.0, which meaans that all the proteins will be used.
#' @param protein_number Number of proteins to test. For example, `protein_number=1000'.
#' Proteins are ranked in decreasing order based on their mean abundance across all the samples
#' and top `protein_number' proteins will be selected to test.
#' Default is 1000.
#' @param FDR a pre-specified false discovery ratio (FDR) to control the overall false positive.
#' Default is 0.05
#' @param power a pre-specified statistical power which defined as the probability of detecting a
#' true fold change. You should input the average of power you expect. Default is 0.9
#' @param width Width of the saved pdf file. Default is 5.
#' @param height Height of the saved pdf file. Default is 5.
#' @param address The name of folder that will store the results. Default folder is the current working directory.
#' The other assigned folder has to be existed under the current working directory.
#' An output pdf file is automatically created with the default name of `HypothesisTestingSampleSizePlot.pdf'.
#' The command address can help to specify where to store the file as well as how to modify the beginning of the file name.
#' If address=FALSE, plot will be not saved as pdf file but showed in window.
#'
#' @importFrom MSstats designSampleSize designSampleSizePlots
#' @importFrom stats quantile
#' @export
#'
#' @return sample size plot for hypothesis testing :
#' the plot for the minimal number of biological replciates per condition to acquire the expected
#' FDR and power under different fold changes.
#' @return data frame with columns desiredFC, numSample, FDR, power and CV
#' @author Ting Huang, Meena Choi, Olga Vitek
#' @examples
#' data(OV_SRM_train)
#' data(OV_SRM_train_annotation)
#'
#' # sample size plot for hypothesis testing
#' HT_res <- designSampleSizeHypothesisTestingPlot(data = OV_SRM_train,
#'                                                 annotation= OV_SRM_train_annotation,
#'                                                 desired_FC = "data",
#'                                                 select_testing_proteins = "proportion",
#'                                                 protein_proportion = 1.0,
#'                                                 protein_number = 1000,
#'                                                 FDR=0.05,
#'                                                 power=0.9)
#'
#' # data frame with columns desiredFC, numSample, FDR, power and CV
#' head(HT_res)
#'
designSampleSizeHypothesisTestingPlot <- function(data,
                                              annotation,
                                              desired_FC = "data",
                                              select_testing_proteins = "proportion",
                                              protein_proportion = 1.0,
                                              protein_number = 1000,
                                              FDR=0.05,
                                              power=0.9,
                                              height = 5,
                                              width = 5,
                                              address = "") {


    ## Estimate the mean abundance and variance of each protein in each phenotype group
    parameters <- estimateVar(data, annotation)

    ###############################################################################
    ## log file
    ## save process output in each step
    loginfo <- .logGeneration()
    finalfile <- loginfo$finalfile
    processout <- loginfo$processout

    processout <- rbind(processout,
                        as.matrix(c(" ", " ", "MSstatsSampleSize - designSampleSizeHypothesisTesting function", " "), ncol=1))

    ###############################################################################

    ## match between data and annotation
    data <- data[, annotation$BioReplicate]
    group <- as.factor(as.character(annotation$Condition))

    ## input checking
    ## 2.1 desired_FC
    if ( !is.element("data", desired_FC) & !length(desired_FC) == 2 ) {
        processout <- rbind(processout, c("ERROR : desired_FC should be `data` or a a vector with
                                          the lower and upper values of the desired fold change.
                                          Please check it."))
        write.table(processout, file=finalfile, row.names=FALSE)

        stop("ERROR : desired_FC should be `data` or a a vector with
                the lower and upper values of the desired fold change.
                Please check it. \n")

    }

    if(!is.element("data", desired_FC)){
        processout <- rbind(processout, paste0(paste0("desired_FC = ", desired_FC), collapse=","))
        write.table(processout, file=finalfile, row.names=FALSE)
    }

    ## 2.2 select_testing_proteins
    if ( !is.element("proportion", select_testing_proteins) & !is.element("number", select_testing_proteins) ) {
        processout <- rbind(processout, c("ERROR : select_testing_proteins should be either `proportion` or `number`. Please check it."))
        write.table(processout, file=finalfile, row.names=FALSE)

        stop("ERROR : select_testing_proteins should be either `proportion` or `number`. Please check it. \n")
    }

    processout <- rbind(processout, c(paste0("select_testing_proteins = ", select_testing_proteins)))
    write.table(processout, file=finalfile, row.names=FALSE)


    ## 2.3 protein_proportion
    if(is.element("proportion", select_testing_proteins)){
        if (is.null(protein_proportion) ) {
            processout <- rbind(processout, c("ERROR : protein_proportion is required for select_testing_proteins=`proportion`. Please provide the value for protein_proportion."))
            write.table(processout, file=finalfile, row.names=FALSE)

            stop("ERROR : protein_proportion is required for select_testing_proteins=`proportion`. Please provide the value for protein_proportion. \n")

        } else if ( protein_proportion < 0 | protein_proportion > 1 ) {
            processout <- rbind(processout, c("ERROR : protein_proportion should be between 0 and 1. Please check the value for protein_proportion."))
            write.table(processout, file=finalfile, row.names=FALSE)

            stop("ERROR : protein_proportion should be between 0 and 1. Please check the value for protein_proportion. \n")
        }

        processout <- rbind(processout, c(paste0("protein_proportion = ", protein_proportion)))
        write.table(processout, file=finalfile, row.names=FALSE)

    }


    ## 2.4 protein_number
    num_total_proteins <- nrow(data)

    if(is.element("number", select_testing_proteins)){
        if (is.null(protein_number) ) {
            processout <- rbind(processout, c("ERROR : protein_number is required for select_testing_proteins=`number`. Please provide the value for protein_number."))
            write.table(processout, file=finalfile, row.names=FALSE)

            stop("ERROR : protein_number is required for select_testing_proteins=`number`. Please provide the value for protein_number. \n")

        } else if ( protein_number < 0 | protein_number > num_total_proteins ) {
            processout <- rbind(processout,
                                c(paste0("ERROR : protein_number should be between 0 and the total number of protein(", num_total_proteins,
                                         "). Please check the value for protein_number.")))
            write.table(processout, file=finalfile, row.names=FALSE)

            stop(message(paste0("ERROR : protein_number should be between 0 and the total number of protein(", num_total_proteins,
                                "). Please check the value for protein_number. \n")))
        }

        processout <- rbind(processout, c(paste0("protein_number = ", protein_number)))
        write.table(processout, file=finalfile, row.names=FALSE)
    }


    ## 2.5 FDR
    if ( (!is.numeric(FDR)) | FDR < 0 | FDR > 1 ) {
        processout <- rbind(processout, c("ERROR : FDR should be a numeric value between 0 and 1."))
        write.table(processout, file=finalfile, row.names=FALSE)

        stop("ERROR : FDR should be a numeric value between 0 and 1. \n")
    }

    processout <- rbind(processout, c(paste0("FDR = ", FDR)))
    write.table(processout, file=finalfile, row.names=FALSE)


    ## 2.6 power
    if ( (!is.numeric(power)) | power < 0 | power > 1 ) {
        processout <- rbind(processout, c("ERROR : power should be a numeric value between 0 and 1."))
        write.table(processout, file=finalfile, row.names=FALSE)

        stop("ERROR : power should be a numeric value between 0 and 1. \n")
    }

    processout <- rbind(processout, c(paste0("power = ", power)))
    write.table(processout, file=finalfile, row.names=FALSE)


    ###############################################################################
    ## Prepare the parameters for hypothesis testing
    models <- parameters$model
    mu <- parameters$mu
    sigma <- parameters$sigma
    promean <- parameters$promean
    proteins <- parameters$protein

    ## Generate the protein number to test
    # 1. based on proportion
    if (select_testing_proteins == "proportion") {
        nproteins <- length(promean)
        protein_num <- round(nproteins*protein_proportion)

        # 2. based on number
    } else {
        protein_num <- protein_number
    }

    ## select proteins based on their mean abundance
    selectedPros <- order(promean, decreasing = TRUE)[seq_len(protein_num)]

    processout <- rbind(processout, c(paste0(" Number of proteins to teset: ", protein_num)))
    write.table(processout, file=finalfile, row.names=FALSE)
    message(" Number of proteins to test: ", protein_num)

    ## prepare the expected fold change and linear models for hypothesis testing
    models_2 <- models[selectedPros]
    mu_2 <- mu[selectedPros, ]
    # sigma_2 <- sigma[selectedPros, ]
    # promean <- parameters$promean[selectedPros]
    # proteins_2 <- proteins[proteins %in% selectedPros]

    ## prepare the expected fold change for hypothesis testing
    if (is.element("data", desired_FC)) {
        FCs <- NULL
        # 1. use the fold change estimated from input data
        for(i in seq(ncol(mu_2)-1)){
            for(j in seq(i+1, ncol(mu_2))){
                FCs <- c(FCs, 2^(mu_2[,i] - mu_2[,j]))

            }
        }

        # make all the fold changes not less than 1
        FCs <- ifelse(FCs > 1, FCs, 1/FCs)
        desired_FC <- round(quantile(FCs, probs = c(0.05, 0.95)), 4)
        if(desired_FC["5%"] < 1.1){
            processout <- rbind(processout, c(paste0(" The 1st quantile of fold changes from input data is ",
                                                     desired_FC["5%"],
                                                     ", which is too small. The minimal fold change for sample size estimation is set to 1.1")))
            write.table(processout, file=finalfile, row.names=FALSE)
            message(" The 1st quantile of fold changes from input data is ", desired_FC["5%"],
                    ", which is too small. The minimal fold change for sample size estimation is set to 1.1")

            # set the minimal fold change to 1.1
            desired_FC["5%"] <- 1.1

        }


    } else {
        # 2. use the predefined fold change range
        desired_FC = desired_FC

    }


    ## Calculate sample size for future experiments:
    # Minimal number of biological replicates per condition
    res <- designSampleSize(data = models_2,
                            numSample = TRUE,
                            desiredFC = desired_FC,
                            FDR = FDR,
                            power = power)

    if (address != FALSE) {
        allfiles <- list.files()

        num <- 0
        plotfilenaming <- paste0(address, "HypothesisTestingSampleSizePlot")
        plotfinalfile <- paste0(address, "HypothesisTestingSampleSizePlot.pdf")

        while (is.element(plotfinalfile, allfiles)) {
            num <- num + 1
            plotfinalfile <- paste0(paste(plotfilenaming, num, sep="-"), ".pdf")
        }

        pdf(plotfinalfile, width=width, height=height)
    }

    designSampleSizePlots(res)

    if (address != FALSE) {
        dev.off()
    }

    message(" Drew the sample size plot for hypothesis testing!")
    processout <- rbind(processout, as.matrix(c(" Drew the sample size plot for hypothesis testing!"), ncol=1))
    write.table(processout, file=finalfile, row.names=FALSE)

    return(res)

}
