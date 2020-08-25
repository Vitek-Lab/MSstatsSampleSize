#' PCA plot for each simulation
#' @details This function draws PCA plot for the whole input dataset and each simulated dataset
#' in `simulations' (input for this function).
#' It outputs the number of simulations plus 1 of PCA plots.
#' The first page shows a PCA plot for the input preliminary dataset.
#' Each of the following pages shows a PCA plot under one simulation.
#' x-axis of PCA plot is the first component and y-axis is the second component.
#' This function can be used to validate whether the simulated dataset looks consistent with the input dataset.
#'
#' @param simulations A list of simulated datasets. It should be the output of \code{\link{simulateDataset}} function.
#' @param which.PCA Select one PCA plot to show. It can be "all", "allonly", or "simulationX".
#' X should be index of simulation, such as "simulation1" or "simulation5".
#' Default is "all", which generates all the plots.
#' "allonly" generates the PCA plot for the whole input dataset.
#' "simulationX" generates the PCA plot for a specific simulated dataset (given by index).
#' @param x.axis.size size of x-axis labeling in PCA Plot. Default is 10.
#' @param y.axis.size size of y-axis labels. Default is 10.
#' @param dot.size size of dots in PCA plot. Default is 3.
#' @param legend.size size of legend above Profile plot. Default is 7.
#' @param width width of the saved pdf file. Default is 6.
#' @param height height of the saved pdf file. Default is 5.
#' @param address the name of folder that will store the results. Default folder is the current working directory.
#' The other assigned folder has to be existed under the current working directory.
#' An output pdf file is automatically created with the default name of `PCAPlot.pdf'.
#' The command address can help to specify where to store the file as well as how to modify the beginning of the file name.
#' If address=FALSE, plot will be not saved as pdf file but showed in window.
#'
#' @return PCA plot : x-axis of PCA plot is the first component and y-axis is the second component.
#' @author Ting Huang, Meena Choi, Olga Vitek
#' @examples
#' data(OV_SRM_train)
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
#' # output a PDF file with multiple PCA plots
#' designSampleSizePCAplot(simulated_datasets)
#'
#' @import ggplot2
#' @import graphics
#' @importFrom utils sessionInfo read.table write.table
#' @importFrom grDevices pdf dev.off
#' @importFrom stats prcomp
#' @export
#'
designSampleSizePCAplot <- function(simulations, which.PCA = "all", save.pdf = T,
                                    ...) {
    
    ###############################################################################
    ## log file
    ## save process output in each step
    dots <- list(...)
    session <- dots$session
    func <- as.list(sys.call())[[1]]
    pc_plot <- list()
    conn <- .find_log(...)
    
    ## parameter checking: which.PCA
    if (!(length(which.PCA) == 1 & (which.PCA %in% c("all", "data") |
                                    (startsWith(which.PCA, "simulation"))))) {
        stop("CALL_",func,"_which.PCA should be one of 'all', 'data' or ",
             "'simulation + index of simulation' (such as 'simulation1').")
    }
    
    .status(sprintf("which.PCA = %s", which.PCA), log = conn$con, func = func)
    
    if(grepl("simulation", which.PCA) | which.PCA == "all"){
        index <- length(simulations$simulation_train_Xs)
        iter <- as.numeric(gsub("[[:alpha:]]","",which.PCA))
        if(!is.na(iter)){
            iters <- iter
        }else{
            iters <- seq_len(index)
        }
        
        for(i in iters){
            title <- ifelse(is.na(iter),
                            sprintf("Simulation: %s", i),
                            sprintf("Simulation: %s", iter))
            pr_comp <- .do_prcomp(sim = simulations$simulation_train_Xs[[i]],
                                  sim_y = simulations$simulation_train_Ys[[i]])
            pc_plot[[i]] <- .pca_plot(data = pr_comp$pc.result, exp_var = pr_comp$exp.var
                                      ,title = title, ...)
            .status(detail = sprintf("Plotted %s of %s", i, length(iters)),
                    log = conn$con, func = func, ...)
        }
        
        if(which.PCA == 'all'){
            pr_comp <- .do_prcomp(simulations$input_X, simulations$input_Y)
            pc_plot_ip <- .pca_plot(data = pr_comp$pc.result, exp_var = pr_comp$exp.var,
                                    title = "Input Dataset", ...)
            pc_plot <- append(list(pc_plot_ip), pc_plot)
        }
    }else if(which.PCA == 'data'){
        pr_comp <- .do_prcomp(simulations$input_X, simulations$input_Y)
        pc_plot <- .pca_plot(data = pr_comp$pc.result, exp_var = pr_comp$exp.var,
                             title = "Input Dataset", ...)
    }else{
        stop("Improper which arguement provided, should be either 'all', 'data' or",
             "'simulation'+index, example ='simulation1'")
    }
    
    
    if(save.pdf | !grepl("simulation", which.PCA)){
        if(!save.pdf){
            warning("CALL_",func,"_All PCA's requested, forcing file save")
        }
        file_name <- file.path(getwd(), sprintf("PCA_Plot_%s.pdf",
                                                format(Sys.time(), "%Y%m%d%H%M")))
        pdf(file = file_name, width = 8.5, height = 8.5)
        for(i in seq_along(pc_plot)){
            print(pc_plot[[i]])
        }
        dev.off()
        .status(detail = sprintf("Plots saved at %s", file_name), log = conn$con, 
                func =func)
    } else{
        return(pc_plot[[iter]])
    }
}
