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
#' `Run' for MS run, `BioReplicate' for biological subject ID and `Condition' for group information are required.
#' `Run' information should be the same with the column of `data'. Multiple `Run' may come from same `BioReplicate'.
#' @param log2Trans Default is FALSE. If TRUE, the input `data' is log-transformed with base 2.
#' @param desired_FC the range of a desired fold change.
#' The first option (Default) is "data",
#' indicating the range of the desired fold change is directly estimated from the input `data',
#' which are the minimal fold change and the maximal fold change in the input `data'.
#' The second option is a vector which includes
#' the lower and upper values of the desired fold change (For example, c(1.25,1.75)).
#' @param protein_rank The standard to rank the proteins in the input `data'.
#' It can be 1) "mean" of protein abundances over all the samples or
#' 2) "sd" (standard deviation) of protein abundances over all the samples or
#' 3) the "combined" of mean abundance and standard deviation.
#' The proteins in the input `data' are ranked based on `protein_rank'
#' and the user can select a subset of proteins for hypothesis testing and sample size calculation.
#' @param protein_select select proteins with "low" or "high" mean abundance or
#' standard deviation (variance) or their combination for hypothesis testing and sample size calculation.
#' The variance (and the range of desired fold change if desiredFC = "data")
#' will be estimated from the selected proteins.
#' If `protein_order = "mean"' or protein_order = "sd"', `protein_select' should be "low" or "high".
#' Default is "high", indicating high abundance or standard deviation proteins are selected.
#' If `protein_order = "combined"', `protein_select' has two elements.
#' The first element corrresponds to the mean abundance.
#' The second element corrresponds to the standard deviation (variance).
#' Default is c("high", "low") (select proteins with high abundance and low variance).
#' @param protein_quantile_cutoff Quantile cutoff(s) for selecting protiens
#' for hypothesis testing and sample size calculation.
#' For example, when `protein_rank="mean"', and
#' `protein_select="high"', `protein_quantile_cutoff=0.1'
#' Proteins are ranked based on their mean abundance across all the samples.
#' Then, the top 10% highest abundant proteins are selected.
#' Default is 0.0, which means that all the proteins are used.
#' If `protein_rank = "combined"', `protein_quantile_cutoff'` has two cutoffs.
#' The first element corrresponds to the cutoff for mean abundance.
#' The second element corrresponds to the cutoff for the standard deviation (variance).
#' Default is c(0.0, 1.0), which means that all the proteins will be used.
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
#'                                                 log2Trans = FALSE,
#'                                                 desired_FC = "data",
#'                                                 protein_rank = "mean",
#'                                                 protein_select = "high",
#'                                                 protein_quantile_cutoff = 0.0,
#'                                                 FDR=0.05,
#'                                                 power=0.9)
#'
#' # data frame with columns desiredFC, numSample, FDR, power and CV
#' head(HT_res)
#'
designSampleSizeHypothesisTestingPlot <- function(data,
                                              annotation,
                                              log2Trans = FALSE,
                                              desired_FC = "data",
                                              protein_rank = "mean",
                                              protein_select = "high",
                                              protein_quantile_cutoff = 0.0,
                                              FDR=0.05,
                                              power=0.9,
                                              height = 5,
                                              width = 5,
                                              address = "") {


    ## Estimate the mean abundance and variance of each protein in each phenotype group
    parameters <- estimateVar(data = data,
                              annotation = annotation,
                              log2Trans = log2Trans)

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
    data <- data[, annotation$Run]
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

    ## 2.2 protein_rank
    if ( !protein_rank %in% c("mean", "sd", "combined") | length(protein_rank) > 1) {
        processout <- rbind(processout, c("ERROR : protein_rank should be `mean` or `sd` or `combined`. Please check it."))
        write.table(processout, file=finalfile, row.names=FALSE)

        stop("protein_rank should be `mean` or `sd` or `combined`. Please check it. \n")
    }

    processout <- rbind(processout, c(paste0("protein_rank = ", protein_rank)))
    write.table(processout, file=finalfile, row.names=FALSE)

    ## 2.3 protein_select
    if(protein_rank == "combined"){
        if(length(protein_select) != 2 | !all(protein_select %in% c("high", "low"))){
            processout <- rbind(processout, c("ERROR : When protein_rank = `combined`, protein_select should have two values,
                 each of which should be either `low` or `high`. Please check it. "))
            write.table(processout, file=finalfile, row.names=FALSE)

            stop("When protein_rank = `combined`, protein_select should have two values, each of which should be either `low` or `high`.
                 Please check it. \n")
        }
    } else{
        if(length(protein_select) != 1 | !all(protein_select %in% c("high", "low"))){
            processout <- rbind(processout,
                                c("ERROR : protein_select should be either `low` or `high`. Please check it. "))
            write.table(processout, file=finalfile, row.names=FALSE)

            stop("protein_select should be either `low` or `high`. Please check it. \n")
        }
    }

    processout <- rbind(processout, paste0(paste0("protein_select = ", protein_select), collapse=","))
    write.table(processout, file=finalfile, row.names=FALSE)

    ## 2.4 protein_quantile_cutoff
    if(protein_rank == "combined"){
        if(length(protein_quantile_cutoff) != 2 |
           any(protein_quantile_cutoff < 0) |
           any(protein_quantile_cutoff > 1)){
            processout <- rbind(processout, c("ERROR : When protein_rank = `combined`, protein_quantile_cutoff should have two values,
                 each of which should be between 0 and 1. Please check it. "))
            write.table(processout, file=finalfile, row.names=FALSE)

            stop("When protein_rank = `combined`, protein_quantile_cutoff should have two values, each of which should be between 0 and 1.
                 Please check it. \n")
        }
    } else{
        if(length(protein_quantile_cutoff) != 1 |
           any(protein_quantile_cutoff < 0) |
           any(protein_quantile_cutoff > 1)){
            processout <- rbind(processout, c("ERROR : protein_quantile_cutoff should be between 0 and 1.
                                              Please check it. "))
            write.table(processout, file=finalfile, row.names=FALSE)

            stop("protein_quantile_cutoff should be between 0 and 1. Please check it. \n")
        }
    }

    processout <- rbind(processout, paste0(paste0("protein_quantile_cutoff = ", protein_quantile_cutoff), collapse=","))
    write.table(processout, file=finalfile, row.names=FALSE)


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
    proteins <- parameters$protein

    ## select the proteins to simulate
    temp <- data.frame(Protein = proteins,
                       Mean = parameters$promean[proteins],
                       SD = parameters$prosd[proteins])

    # 1. use both mean and standard deviation cutoff
    if(protein_rank == "combined"){
        mean_cutoff <- quantile(temp$Mean, probs = protein_quantile_cutoff[1], na.rm = TRUE)
        sd_cutoff <- quantile(temp$SD, probs = protein_quantile_cutoff[2], na.rm = TRUE)
        if(protein_select[1] == "low" & protein_select[2] == "low") { # low abundance and low variance
            selectedPros <- temp[temp$Mean <= mean_cutoff & temp$SD <= sd_cutoff, "Protein"]
        } else if(protein_select[1] == "low" & protein_select[2] == "high") { # low abundance and high variance
            selectedPros <- temp[temp$Mean <= mean_cutoff & temp$SD > sd_cutoff, "Protein"]
        } else if(protein_select[1] == "high" & protein_select[2] == "low") { # high abundance and low variance
            selectedPros <- temp[temp$Mean > mean_cutoff & temp$SD <= sd_cutoff, "Protein"]
        } else if(protein_select[1] == "high" & protein_select[2] == "high") { # high abundance and high variance
            selectedPros <- temp[temp$Mean > mean_cutoff & temp$SD > sd_cutoff, "Protein"]
        }

        # 2. use mean cutoff
    } else if(protein_rank == "mean"){
        mean_cutoff <- quantile(temp$Mean, probs = protein_quantile_cutoff, na.rm = TRUE)
        if(protein_select == "low") { # low protein abundance
            selectedPros <- temp[temp$Mean <= mean_cutoff, "Protein"]
        } else if(protein_select == "high") { # high protein abundacne
            selectedPros <- temp[temp$Mean > mean_cutoff, "Protein"]
        }

        # 3. use standard deviation cutoff
    } else if(protein_rank == "sd"){
        sd_cutoff <- quantile(temp$SD, probs = protein_quantile_cutoff, na.rm = TRUE)
        if(protein_select == "low") { # low variance
            selectedPros <- temp[temp$SD <= sd_cutoff, "Protein"]
        } else if(protein_select == "high") { # high variance
            selectedPros <- temp[temp$SD > sd_cutoff, "Protein"]
        }

    }

    protein_num <- length(selectedPros)
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
