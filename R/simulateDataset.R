#' Simulate datasets with the given numbers of biological replicates and proteins based on the preliminary \emph{data}
#' @details This function simulate datasets with the given numbers of biological replicates and proteins based on the preliminary data (input for this function).
#' The function fits intensity-based linear model on the input prelimiary data `data` in order to get variance and mean abundance, using `estimateVar` function.
#' Then it uses variance components and mean abundance to simulate new training data with the given sample size and protein number.
#' It outputs the number of simulated proteins, a vector with the number of simulated samples in a condition,
#' the list of simulated training datasets and the (simulated) validation dataset.

#'
#' @param data Protein abundance data matrix. Rows are proteins and columns are biological replicates(samples).
#' @param annotation Group information for samples in data. `BioReplicate` for sample ID and `Condition` for group information are required.
#' `BioReplicate` information should match with column names of `data`.
#' @param num_simulations Number of times to repeat simulation experiments (Number of simulated datasets). Default is 10.
#' @param expected_FC Expected fold change of proteins. The first option (Default) is `data`, indicating the fold changes are directly estimated from the input `data`.
#' The second option is a vector with predefined fold changes of listed proteins. The vector names must match with the unique information of Condition in `annotation`.
#' One group must be selected as a baseline and has fold change 1 in the vector.
#' The user should provide list_diff_proteins, which users expect to have the fold changes greater than 1.
#' Other proteins that are not available in `list_diff_proteins` will be expected to have fold change = 1
#' @param list_diff_proteins Vector of proteins names which are set to have fold changes greater than 1 between conditions.
#' If user selected `expected_FC= `"data `" `, this should be NULL.
#' @param select_simulated_proteins The standard to select the simulated proteins among data.
#' It can be 1) `proportion` of total number of proteins in the input data or 2) `number` to specify the number of proteins.
#' `proportion` indicates that user should provide the value for `protein_proportion` option.
#' `number` indicates that user should provide the value for `protein_number` option.
#' @param protein_proportion Proportion of total number of proteins in the input data to simulate.
#' For example, input data has 1,000 proteins and user selects `protein_proportion=0.1`.
#' Proteins are ranked in decreasing order based on their mean abundance across all the samples.
#' Then, 1,000 * 0.1 = 100 proteins will be selected from the top list to simulate.
#' Default is 1.0, which meaans that all the proteins will be used.
#' @param protein_number Number of proteins to simulate. For example, `protein_number=1000`.
#' Proteins are ranked in decreasing order based on their mean abundance across all the samples and top `protein_number` proteins will be selected to simulate.
#' Default is 1000.
#' @param samples_per_group Number of samples per group to simulate. Default is 50.
#' @param simulate_validation Default is FALSE. If TRUE, simulate the validation set; otherwise, the input `data` will be used as the validation set.
#' @param valid_samples_per_group Number of validation samples per group to simulate. This option works only when user selects `simulate_validation=TRUE`. Default is 50.
#'
#' @return \emph{num_proteins} the number of simulated proteins.
#' @return \emph{num_samples} a vector with the number of simulated samples in each condition.
#' @return \emph{simulation_train_Xs} is the list of simulated protein abundance matrices. Each element of the list represents one simulation.
#' @return \emph{simulation_train_Ys} is the list of simulated condition vectors. Each element of the list represents one simulation.
#' @return \emph{valid_X} is the validation protein abundance matrix, which is used for classification.
#' @return \emph{valid_Y} is the condition vector of validation samples.
#'
#' @author Ting Huang, Meena Choi, Olga Vitek.
#'
#' @examples
#' data(OV_SRM_train)
#' data(OV_SRM_train_annotation)
#'
#' # num_simulations = 10: simulate 10 times
#' # expected_FC = "data": fold change estimated from OV_SRM_train
#' # select_simulated_proteins = "proportion":
#' # select the simulated proteins based on the proportion of total proteins
#' # simulate_validation = FALSE: use input OV_SRM_train as validation set
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
#'                                       simulate_validation = FALSE,
#'                                       valid_samples_per_group = 50)
#'
#' # the number of simulated proteins
#' simulated_datasets$num_proteins
#'
#' # a vector with the number of simulated samples in each condition
#' simulated_datasets$num_samples
#'
#' # the list of simulated protein abundance matrices
#' # Each element of the list represents one simulation
#' head(simulated_datasets$simulation_train_Xs[[1]]) # first simulation
#'
#'# the list of simulated condition vectors
#'# Each element of the list represents one simulation
#'head(simulated_datasets$simulation_train_Ys[[1]]) # first simulation
#'
#' @export
#'
#' @importFrom utils sessionInfo read.table write.table
#'
simulateDataset <- function(data,
                            annotation,
                            num_simulations = 10,
                            expected_FC = "data",
                            list_diff_proteins =  NULL,
                            select_simulated_proteins = "proportion",
                            protein_proportion = 1.0,
                            protein_number = 1000,
                            samples_per_group = 50,
                            simulate_validation = FALSE,
                            valid_samples_per_group = 50) {

    ## Estimate the mean abundance and variance of each protein in each phenotype group
    parameters <- estimateVar(data, annotation)

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
                        as.matrix(c(" ", " ", "MSstatsSampleSize - simulateDataset function", " "), ncol=1))


    ###############################################################################
    ## Input and option checking

    data <- data[, annotation$BioReplicate]
    group <- annotation$Condition


    ## 2. check input for option
    ## 2.1 number of simulation : any lower limit?
    if (num_simulations < 10) {
        processout <- rbind(processout, c(paste0("ERROR : num_simulation =", num_simulations, ", which is smaller than 10. - stop")))
        write.table(processout, file=finalfile, row.names=FALSE)

        stop("Please use more than 10 simulations. \n")
    }

    processout <- rbind(processout, c(paste0("Number of simulation = ", num_simulations)))
    write.table(processout, file=finalfile, row.names=FALSE)


    ## 2.2 expected_FC and list_diff_proteins
    if ( !is.element("data", expected_FC) & !is.element(1, expected_FC) ) {
        processout <- rbind(processout, c("ERROR : expected_FC should be `data` or a vector including 1. Please check it."))
        write.table(processout, file=finalfile, row.names=FALSE)

        stop("ERROR : expected_FC should be `data` or a vector including 1. Please check it. \n")
    }

    if(!is.element("data", expected_FC)){
        processout <- rbind(processout, paste0(paste0("expected_FC = ", expected_FC), collapse=","))
        write.table(processout, file=finalfile, row.names=FALSE)
    }


    ## 2.3 list_diff_proteins
    if ( is.numeric(expected_FC) & is.null(list_diff_proteins) ) {
        processout <- rbind(processout, c("ERROR : list_diff_proteins are required for predefined expected_FC. Please provide the vector for list_diff_proteins."))
        write.table(processout, file=finalfile, row.names=FALSE)

        stop("ERROR : list_diff_proteins are required for predefined expected_FC. Please provide the vector for list_diff_proteins. \n")
    }

    if(!is.element("data", expected_FC)){
        processout <- rbind(processout, paste0("list_diff_proteins = ", paste0(list_diff_proteins, collapse=",")))
        write.table(processout, file=finalfile, row.names=FALSE)
    }


    ## 2.4 select_simulated_proteins
    if ( !is.element("proportion", select_simulated_proteins) & !is.element("number", select_simulated_proteins) ) {
        processout <- rbind(processout, c("ERROR : select_simulated_protein should be either `proportion` or `number`. Please check it."))
        write.table(processout, file=finalfile, row.names=FALSE)

        stop("ERROR : select_simulated_protein should be either `proportion` or `number`. Please check it. \n")
    }

    processout <- rbind(processout, c(paste0("select_simulated_proteins = ", select_simulated_proteins)))
    write.table(processout, file=finalfile, row.names=FALSE)


    ## 2.5 protein_proportion
    if(is.element("proportion", select_simulated_proteins)){
        if (is.null(protein_proportion) ) {
            processout <- rbind(processout, c("ERROR : protein_proportion is required for select_simulated_protein=`proportion`. Please provide the value for protein_proportion."))
            write.table(processout, file=finalfile, row.names=FALSE)

            stop("ERROR : protein_proportion is required for select_simulated_protein=`proportion`. Please provide the value for protein_proportion. \n")

        } else if ( protein_proportion < 0 | protein_proportion > 1 ) {
            processout <- rbind(processout, c("ERROR : protein_proportion should be between 0 and 1. Please check the value for protein_proportion."))
            write.table(processout, file=finalfile, row.names=FALSE)

            stop("ERROR : protein_proportion should be between 0 and 1. Please check the value for protein_proportion. \n")
        }

        processout <- rbind(processout, c(paste0("protein_proportion = ", protein_proportion)))
        write.table(processout, file=finalfile, row.names=FALSE)

    }



    ## 2.6 protein_number
    num_total_proteins <- nrow(data)

    if(is.element("number", select_simulated_proteins)){
        if (is.null(protein_number) ) {
            processout <- rbind(processout, c("ERROR : protein_number is required for select_simulated_protein=`number`. Please provide the value for protein_number."))
            write.table(processout, file=finalfile, row.names=FALSE)

            stop("ERROR : protein_number is required for select_simulated_protein=`number`. Please provide the value for protein_number. \n")

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

    ## 2.7 samples_per_group
    if ( !is.numeric(samples_per_group) ) {
        processout <- rbind(processout, c("ERROR : sample_per_group should be numeric. Please provide the numeric value for samples_per_group."))
        write.table(processout, file=finalfile, row.names=FALSE)

        stop("ERROR : sample_per_group should be numeric. Please provide the numeric value for samples_per_group. \n")

    } else if ( samples_per_group%%1 != 0 ) { ## not integer, then round

        samples_per_group <- round(samples_per_group)

        processout <- rbind(processout, c("NOTE : samples_per_group should be integer. Rounded samples_per_group will be used."))
        write.table(processout, file=finalfile, row.names=FALSE)

        message("NOTE : samples_per_group should be integer. Rounded samples_per_group will be used.")
    }

    processout <- rbind(processout, c(paste0("samples_per_group = ", samples_per_group)))
    write.table(processout, file=finalfile, row.names=FALSE)


    ## 2.8 simulated value
    if ( !is.logical(simulate_validation) ) {
        processout <- rbind(processout, c("ERROR : simulate_validation should be logical. Please provide either TRUE or FALSE for simulate_validation."))
        write.table(processout, file=finalfile, row.names=FALSE)

        stop("ERROR : simulate_validation should be logical. Please provide either TRUE or FALSE for simulate_validation. \n")
    }

    processout <- rbind(processout, c(paste0("simulate_validation = ", simulate_validation)))
    write.table(processout, file=finalfile, row.names=FALSE)


    ## 2.9 valid_samples_per_group
    if ( simulate_validation & (is.null(valid_samples_per_group) | !is.numeric(valid_samples_per_group)) ) {
        processout <- rbind(processout, c("ERROR : valid_samples_per_group is required for simulate_validation=TRUE. Please provide the numeric value for valid_samples_per_group."))
        write.table(processout, file=finalfile, row.names=FALSE)

        stop("ERROR : valid_samples_per_group is required for simulate_validation=TRUE. Please provide the numeric value for valid_samples_per_group. \n")
    }

    processout <- rbind(processout, c(paste0("valid_samples_per_group = ", valid_samples_per_group)))
    write.table(processout, file=finalfile, row.names=FALSE)


    ###############################################################################
    processout <- rbind(processout, c("Preparing simulation...."))
    write.table(processout, file=finalfile, row.names=FALSE)
    message(" Preparing simulation...")



    ## Prepare the parameters for simulation experiment
    mu <- parameters$mu
    sigma <- parameters$sigma
    promean <- parameters$promean
    proteins <- parameters$protein

    ## prepare the mean and variance matrix for simulation
    # 1. use the fold change estimated from input data
    if (is.element("data", expected_FC)) {
        sim_mu <- mu
        sim_sigma <- sigma

    # 2. use the user-defined fold change
    } else {
        sim_mu <- mu
        sim_sigma <- sigma

        baseline <- names(expected_FC)[expected_FC == 1] # baseline group
        otherlines <- names(expected_FC)[expected_FC != 1] # compared groups

        for(i in 1:length(otherlines)){
            # set group mean for differential proteins based on predefined fold change
            sim_mu[rownames(sim_mu) %in% list_diff_proteins, otherlines[i]] <-
            sim_mu[rownames(sim_mu) %in% list_diff_proteins, baseline] * expected_FC[otherlines[i]]
            # set group mean for non-differential proteins (fold change = 1)
            sim_mu[!rownames(sim_mu) %in% list_diff_proteins, otherlines[i]] <-
            sim_mu[!rownames(sim_mu) %in% list_diff_proteins, baseline]
        }
    }

    ## Generate the training sample size to simulate
    ngroup <- length(unique(group)) # Number of phenotype groups
    num_samples <- rep(samples_per_group, ngroup)
    names(num_samples) <- unique(group)
    train_size <- samples_per_group * ngroup

    processout <- rbind(processout, c(paste0(" Size of training data to simulate: ", train_size)))
    write.table(processout, file=finalfile, row.names=FALSE)
    message(" Size of training data to simulate: ", train_size)

    ## Generate the protein number to simulate
    # 1. based on proportion
    if (select_simulated_proteins == "proportion") {
        nproteins <- nrow(mu)
        protein_num <- round(nproteins*protein_proportion)

    # 2. based on number
    } else {
        protein_num <- protein_number
    }

    ## select proteins based on their mean abundance
    selectedPros <- order(promean, decreasing = TRUE)[1:protein_num]

    ## prepare the mean and variance for simulation
    mu_2 <- mu[selectedPros, ]
    sigma_2 <- sigma[selectedPros, ]

    ## prepare the validation set
    # 1.  simulate samples in the validation data
    if (simulate_validation) {

        valid <- .sampleSimulation(m = valid_samples_per_group,
                                   mu = mu_2,
                                   sigma = sigma_2)
        valid_X <- as.data.frame(valid$X)
        valid_Y <- as.factor(valid$Y)
    } else{ # use input data as validation set
        ## !! impute the missing values by randomly selecting values for each protein
        valid_X <- as.data.frame(apply(data[selectedPros, ],1, function(x) .random.imp(x))) # impute missing values
        valid_Y <- as.factor(group)
    }

    processout <- rbind(processout, c(paste0(" Number of proteins to simulate: ", protein_num)))
    write.table(processout, file=finalfile, row.names=FALSE)
    message(" Number of proteins to simulate: ", protein_num)

    processout <- rbind(processout, c(" Start to run the simulation..."))
    write.table(processout, file=finalfile, row.names=FALSE)
    message(" Start to run the simulation...")

    simulation_train_Xs <- list()
    simulation_train_Ys <- list()

    for (i in 1:num_simulations) { ## Number of simulations
        message("  Simulation: ", i)

        ## simulate samples in the training data
        train <- .sampleSimulation(m = samples_per_group,
                                   mu = mu_2,
                                   sigma = sigma_2)
        X <- as.data.frame(train$X)
        Y <- as.factor(train$Y)

        simulation_train_Xs[[paste("Simulation", i, sep="")]] <- X
        simulation_train_Ys[[paste("Simulation", i, sep="")]] <- Y
    }

    processout <- rbind(processout, c(" Simulation completed."))
    write.table(processout, file=finalfile, row.names=FALSE)
    message(" Simulation completed.")

    return(list(num_proteins = protein_num, # number of proteins
                num_samples = num_samples, # number of samples per group
                simulation_train_Xs = simulation_train_Xs,
                simulation_train_Ys = simulation_train_Ys,
                valid_X = valid_X,
                valid_Y = valid_Y))

}

