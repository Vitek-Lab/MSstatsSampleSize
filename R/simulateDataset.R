#' Simulate datasets with the given number of biological replicates and
#' proteins based on the input \emph{data}
#' @details This function simulate datasets with
#' the given numbers of biological replicates and
#' proteins based on the input dataset (input for this function).
#' The function fits intensity-based linear model on the input \emph{data}
#' in order to get variance and mean abundance, using \code{\link{estimateVar}} function.
#' Then it uses variance components and mean abundance to simulate new training data
#' with the given sample size and protein number.
#' It outputs the number of simulated proteins,
#' a vector with the number of simulated samples in a condition,
#' the list of simulated training datasets,
#' the input preliminary dataset and
#' the (simulated) validation dataset.
#'
#' @param data Protein abundance data matrix.
#' Rows are proteins and columns are biological replicates (samples).
#' @param annotation Group information for samples in data.
#' `Run' for MS run, `BioReplicate' for biological subject ID and `Condition' for group information are required.
#' `Run' information should be the same with the column of `data'. Multiple `Run' may come from same `BioReplicate'.
#' @param log2Trans Default is FALSE. If TRUE, the input `data' is log-transformed with base 2.
#' @param num_simulations Number of times to repeat simulation experiments
#' (Number of simulated datasets). Default is 10.
#' @param samples_per_group Number of samples per group to simulate. Default is 50.
#' @param protein_rank The standard to rank the proteins in the input `data'.
#' It can be 1) "mean" of protein abundances over all the samples or
#' 2) "sd" (standard deviation) of protein abundances over all the samples or
#' 3) the "combined" of mean abundance and standard deviation.
#' The proteins in the input `data' are ranked based on `protein_rank'
#' and the user can select a subset of proteins to simulate.
#' @param protein_select select proteins with "low" or "high" mean abundance or
#' standard deviation (variance) or their combination to simulate.
#' If `protein_order = "mean"' or protein_order = "sd"', `protein_select' should be "low" or "high".
#' Default is "high", indicating high abundance or standard deviation proteins are selected to simulate.
#' If `protein_order = "combined"', `protein_select' has two elements.
#' The first element corrresponds to the mean abundance.
#' The second element corrresponds to the standard deviation (variance).
#' Default is c("high", "low") (select proteins with high abundance and low variance).
#' @param protein_quantile_cutoff Quantile cutoff(s) for selecting protiens to simulate.
#' For example, when `protein_rank="mean"', and
#' `protein_select="high"', `protein_quantile_cutoff=0.1'
#' Proteins are ranked based on their mean abundance across all the samples.
#' Then, the top 10% highest abundant proteins are selected to simulate.
#' Default is 0.0, which means that all the proteins are used.
#' If `protein_rank = "combined"', `protein_quantile_cutoff'` has two cutoffs.
#' The first element corrresponds to the cutoff for mean abundance.
#' The second element corrresponds to the cutoff for the standard deviation (variance).
#' Default is c(0.0, 1.0), which means that all the proteins will be used.
#' @param expected_FC Expected fold change of proteins.
#' The first option (Default) is "data",
#' indicating the fold changes are directly estimated from the input `data'.
#' The second option is a vector with predefined fold changes of listed proteins.
#' The vector names must match with the unique information of Condition in `annotation'.
#' One group must be selected as a baseline and has fold change 1 in the vector.
#' The user should provide list_diff_proteins, which users expect to have the fold changes greater than 1.
#' Other proteins that are not available in `list_diff_proteins' will be expected to have fold change = 1
#' @param list_diff_proteins Vector of proteins names
#' which are set to have fold changes greater than 1 between conditions.
#' If user selected `expected_FC= "data" ', this should be NULL.
#' @param simulate_validation Default is FALSE. If TRUE, simulate the validation set;
#' otherwise, the input `data' will be used as the validation set.
#' @param valid_samples_per_group Number of validation samples per group to simulate.
#' This option works only when user selects `simulate_validation=TRUE'. Default is 50.
#'
#' @return \emph{num_proteins} is the number of simulated proteins.
#'  It should be set up by parameters, named \emph{protein_proportion} or \emph{protein_number}
#' @return \emph{num_samples} is a vector with the number of simulated samples in each condition.
#' It should be same as the parameter, \emph{samples_per_group}
#' @return \emph{input_X} is the input protein abundance matrix `data'.
#' @return \emph{input_Y} is the condition vector for the input `data.
#' @return \emph{simulation_train_Xs} is the list of simulated protein abundance matrices.
#' Each element of the list represents one simulation.
#' @return \emph{simulation_train_Ys} is the list of simulated condition vectors.
#' Each element of the list represents one simulation.
#' @return \emph{valid_X} is the validation protein abundance matrix, which is used for classification.
#' @return \emph{valid_Y} is the condition vector of validation samples.
#'
#' @author Ting Huang, Meena Choi, Sumedh Sankhe, Olga Vitek.
#'
#' @examples
#' data(OV_SRM_train)
#' data(OV_SRM_train_annotation)
#'
#' # num_simulations = 10: simulate 10 times
#' # protein_rank = "mean", protein_select = "high", and protein_quantile_cutoff = 0.0:
#' # select the proteins with high mean abundance based on the protein_quantile_cutoff
#' # expected_FC = "data": fold change estimated from OV_SRM_train
#' # simulate_validation = FALSE: use input OV_SRM_train as validation set
#' # valid_samples_per_group = 50: 50 samples per condition
#' simulated_datasets <- simulateDataset(data = OV_SRM_train,
#'                                       annotation = OV_SRM_train_annotation,
#'                                       log2Trans = FALSE,
#'                                       num_simulations = 10,
#'                                       samples_per_group = 50,
#'                                       protein_rank = "mean",
#'                                       protein_select = "high",
#'                                       protein_quantile_cutoff = 0.0,
#'                                       expected_FC = "data",
#'                                       list_diff_proteins =  NULL,
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
#' head(simulated_datasets$simulation_train_Ys[[1]]) # first simulation
#'
#' @export
#'
#' @importFrom utils sessionInfo read.table write.table
#'
simulateDataset <- function(data, annotation, log2Trans = FALSE,
                            num_simulations = 10, samples_per_group = 50,
                            protein_select = list(Mean = "high", SD = "high"),
                            protein_quantile_cutoff = list(Mean = 0, SD = 0),
                            expected_FC = "data", list_diff_proteins =  NULL,
                            simulate_validation = FALSE,
                            valid_samples_per_group = 50, ...) {
    
    #### 1. log file save process output in each step ####
    dots <- list(...)
    session <- dots$session
    func <- as.list(sys.call())[[1]]
    conn <- .find_log(...)
    
    ##### Simulation Process ####
    res <- .catch_faults({
        
        if(is.null(dots$parameters)){
            parameters <- estimateVar(data = data, annotation = annotation, 
                                      log2Trans = log2Trans) 
        }else{
            .status(detail = "Using provided estimated means and variance",
                    log = conn$con, ...)
            parameters <- dots$parameters
        }
        
        
        data <- data[, annotation$Run]
        group <- as.factor(as.character(annotation$Condition))
        
        ##### 2. check input for option ####
        #### 2.1 number of simulation : any lower limit? ####
        if(num_simulations < 10){
            stop("CALL_", func, "_Please use more than 10 simulations.")
        }
        .status(detail = sprintf("Number of Simulations = %s", num_simulations),
                log = conn$con, ...)
        
        #### 2.2-2.5 samples_per_group ####
        if(!is.numeric(samples_per_group)){
            stop("CALL_",func,"_sample_per_group should be numeric. Please",
                 "provide the numericvalue for samples_per_group.")
        } else if (samples_per_group%%1 != 0){ ## not integer, then round
            samples_per_group <- round(samples_per_group)
            .status(detail = "samples_per_group should be integer. 
                Rounded samples_per_group will be used.", log = conn$con,
                    func = func, ...)
        }
        .status(detail = sprintf("samples_per_group = %s", samples_per_group),
                log = conn$con, ...)
        
        
        if(class(protein_quantile_cutoff) != "list"){
          stop("CALL_",func,"_protein_quantile_cutoff should be provided as :",
               " list(Mean = 0, SD = 0)")  
        }
        
        sapply(protein_quantile_cutoff, function(x){
          if(x<0 | x>1){
              stop("Quantile cutoff values should be between 0 and 1")
          }
        })
        
        if(class(protein_select) != "list"){
            stop("CALL_",func,"_protein_select should be provided as :",
                 " list(Mean = 'high', SD = 'low')")
        }
        
        sapply(protein_select, function(x){
            if(!x %in% c('high','low')){
                stop("CALL_",func,"_protein_select values should be either",
                     " 'high' or 'low'")
            }
        })
        
        
        if(all(c("Mean", "SD") %in% names(protein_quantile_cutoff))){
            .status(detail = "Combined quantile cutoff's provided", log = conn$con)
            .status(detail = sprintf("Mean quantile cutoff = %s", 
                                    protein_quantile_cutoff[['Mean']]),
                    log = conn$con, ...)
            .status(detail = sprintf("SD quantile cutoff = %s",
                                     protein_quantile_cutoff[['SD']]),
                    log = conn$con, ...)
        }else if(names(protein_quantile_cutoff) == "Mean"){
            .status(detail = sprintf("Mean quantile cutoff = %s", 
                                     protein_quantile_cutoff[['Mean']]),
                    log = conn$con, ...)
        }else if(names(protein_quantile_cutoff) == "SD"){
            .status(detail = sprintf("SD quantile cutoff = %s",
                                     protein_quantile_cutoff[['SD']]),
                    log = conn$con, ...)
        }
        .status(detail = sprintf("protien_select = (%s)",
                                 paste(protein_select, collapse = ", ")),
                log = conn$con, ...)
       
        #### 2.6 expected_FC and list_diff_proteins ####
        if ( !is.element("data", expected_FC) & !is.element(1, expected_FC) ) {
            stop("CALL_", func,"_expected_FC should be `data` or a vector including 1.",
                 "Please check it.")
        }
        
        if(!is.element("data", expected_FC)){
            .status(detail = sprintf("expected_FC = %s", expected_FC),
                    log = conn$con, ...)
        }
        
        #### 2.7 list_diff_proteins ####
        if(is.numeric(expected_FC) & is.null(list_diff_proteins)) {
            stop("CALL_",func,"list_diff_proteins are required for predefined",
                 "expected_FC. Please provide the vector for list_diff_proteins.")
        }
        
        if(!is.element("data", expected_FC)){
            .status(detail = sprintf("list_diff_proteins = %s", 
                                     paste(list_diff_proteins, collapse = ",")),
                    log = conn$con, ...)
        }
        
        #### 2.8 simulated value ####
        if(!is.logical(simulate_validation)) {
            stop("CALL_",func,"simulate_validation should be logical.",
                 "Please provide either TRUE or FALSE for simulate_validation.")
        }
        .status(detail = sprintf("simulate_validaiton = %s", simulate_validation),
                log = conn$con, ...)
        
        #### 2.9 valid_samples_per_group ####
        if(simulate_validation & (is.null(valid_samples_per_group) | 
                                  !is.numeric(valid_samples_per_group))) {
            stop("CALL_",func,"valid_samples_per_group is required for",
                 "simulate_validation=TRUE. Please provide the numeric value",
                 "for valid_samples_per_group.")
        }
        .status(detail = sprintf("valid_samples_per_group = %s",
                                 valid_samples_per_group),
                log = conn$con, ...)
        .status("Preparing simulation...", log = conn$con, ...)
        #### 3. Prepare the parameters for simulation experiment ####
        mu <- parameters$mu
        sigma <- parameters$sigma
        proteins <- parameters$protein
        
        ## select the proteins to simulate
        temp <- data.frame(Protein = proteins,
                           Mean = parameters$promean[proteins],
                           SD = parameters$prosd[proteins])
        
        selectedPros <- .protein_select(temp, cutoff = protein_quantile_cutoff,
                                        equality = protein_select, conn = conn)
        
        
        ## prepare the mean and variance matrix for simulation
        # 1. use the fold change estimated from input data
        if (is.element("data", expected_FC)) {
            sim_mu <- mu[selectedPros, ]
            sim_sigma <- sigma[selectedPros, ]
            
            # 2. use the user-defined fold change
        } else {
            sim_mu <- mu[selectedPros, ]
            sim_sigma <- sigma[selectedPros, ]
            
            baseline <- names(expected_FC)[expected_FC == 1] # baseline group
            otherlines <- names(expected_FC)[expected_FC != 1] # compared groups
            
            if(!all(list_diff_proteins %in%  selectedPros)){
                selected_diff_proteins <- intersect(list_diff_proteins, selectedPros)
                
                if (length(selected_diff_proteins) == 0) {
                    stop("None of list_diff_proteins are selected to simulate based",
                         "on protein_quantile_cutoff. Please check it.")
                }
                
                .status(detail= sprintf("NOTE : %s%% list_diff_proteins are selected to simulate based on protein_quantile_cutoff", 
                                        round(length(selected_diff_proteins)/length(list_diff_proteins), 4)*100),
                        log = conn$con, ...)
            }
            
            for(i in seq_along(otherlines)){
                # set group mean for differential proteins based on predefined fold change
                sim_mu[rownames(sim_mu) %in% selected_diff_proteins, otherlines[i]] <-
                    sim_mu[rownames(sim_mu) %in% selected_diff_proteins, baseline] * expected_FC[otherlines[i]]
                # set group mean for non-differential proteins (fold change = 1)
                sim_mu[!rownames(sim_mu) %in% selected_diff_proteins, otherlines[i]] <-
                    sim_mu[!rownames(sim_mu) %in% selected_diff_proteins, baseline]
            }
        }
        
        ## Generate the training sample size to simulate
        ngroup <- length(unique(group)) # Number of phenotype groups
        num_samples <- rep(samples_per_group, ngroup)
        names(num_samples) <- unique(group)
        train_size <- samples_per_group * ngroup
        
        .status(detail = sprintf("Size of training data to simulate: %s", train_size),
                log = conn$con, ...)
        
        ## prepare the validation set
        # 1.  simulate samples in the validation data
        if (simulate_validation) {
            valid <- .sampleSimulation(m = valid_samples_per_group,
                                       mu = sim_mu,
                                       sigma = sim_sigma)
            valid_X <- as.data.frame(valid$X)
            valid_Y <- as.factor(valid$Y)
            
        } else{ # use input data as validation set
            ## !! impute the missing values by randomly selecting values for each protein
            valid_X <- as.data.frame(apply(data[selectedPros, ],1, function(x){
                .randomImputation(x)  
            })) # impute missing values
            valid_Y <- as.factor(group)
        }
        
        protein_num <- length(selectedPros)
        .status(detail = sprintf("Number of proteins to simulate: %s", protein_num),
                log = conn$con, ...)
        .status(detail = "Start to run the simulation", log = conn$con, ...)
        
        simulation_train_Xs <- list()
        simulation_train_Ys <- list()
        
        for (i in seq_len(num_simulations)) { ## Number of simulations
            .status(detail = sprintf("Simulation: %s", i), log = conn$con,
                    func = func, ...)
            
            ## simulate samples in the training data
            train <- .sampleSimulation(m = samples_per_group,
                                       mu = sim_mu,
                                       sigma = sim_sigma)
            X <- as.data.frame(train$X)
            Y <- as.factor(train$Y)
            
            simulation_train_Xs[[paste("Simulation", i, sep="")]] <- X
            simulation_train_Ys[[paste("Simulation", i, sep="")]] <- Y
        }
        
        .status(detail = "Simulation completed.", log = conn$con, ...)
        
        #TODO create structure of a fixed class, see ?structure
        list(num_proteins = protein_num, # number of proteins
             num_samples = num_samples, # number of samples per group
             simulation_train_Xs = simulation_train_Xs,
             simulation_train_Ys = simulation_train_Ys,
             input_X = t(data),
             input_Y = group,
             valid_X = valid_X,
             valid_Y = valid_Y)
    }, conn = conn, session = session)
    #### Return ####
    return(res)
}
