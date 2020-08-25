#############################################
## designSampleSizeClassificationPlots
#############################################
#' Visualization for sample size calculation in classification
#'
#' @description To illustrate the mean classification accuracy and protein 
#' importance under different sample sizes
#' through predictive accuracy plot and protein importance plot.
#'
#' @details This function visualizes for sample size calculation in classification.
#' Mean predictive accuracy and mean protein importance under each sample size 
#' is from the input `data', which is the output from function 
#' \code{\link{designSampleSizeClassification}}.
#'
#' To illustrate the mean predictive accuracy and protein importance under 
#' different sample sizes, it generates two types of plots in pdf files as 
#' output: (1) The predictive accuracy plot, The X-axis represents different 
#' sample sizes and y-axis represents the mean predictive accuracy.
#' The reported sample size per condition can be used to design future experiment
#'
#' (2) The protein importance plot includes multiple subplots.
#' The number of subplots is equal to `list_samples_per_group'.
#' Each subplot shows the top `num_important_proteins_show` most important 
#' proteins under each sample size. The Y-axis of each subplot is the protein 
#' name and X-axis is the mean protein importance under the sample size.
#'
#' @param data A list of outputs from function \code{\link{designSampleSizeClassification}}. Each element represents the results under a specific sample size.
#' The input should include at least two simulation results with different sample sizes.
#' @param optimal_threshold The maximal cutoff for deciding the optimal sample 
#' size. Default is 0.0001. Large cutoff can lead to smaller optimal sample size
#' whereas small cutoff produces large optimal sample size.
#' @param num_important_proteins_show The number of proteins to show in protein 
#' importance plot.
#' @param protein_importance_plot TRUE(default) draws protein importance plot.
#' @param predictive_accuracy_plot TRUE(default) draws predictive accuracy plot.
#' @param ... Arguements that can be passed to ggplot2::theme functions to alter 
#' the visuals
#' @return predictive accuracy plot is the mean predictive accuracy under 
#' different sample sizes.
#' The X-axis represents different sample sizes and y-axis represents the mean 
#' predictive accuracy.
#' @return protein importance plot includes multiple subplots. The number of 
#' subplots is equal to `list_samples_per_group`. Each subplot shows the top 
#' `num_important_proteins_show` most important proteins under each sample size.
#' The Y-axis of each subplot is the protein name and X-axis is the mean protein
#' importance under the sample size.
#' @return a numeric value which is the estimated optimal sample size per group
#' for the input dataset for classification problem.
#' @author Ting Huang, Meena Choi, Olga Vitek.
#' @examples
#' data(OV_SRM_train)
#' data(OV_SRM_train_annotation)
#'
#' # simulate different sample sizes
#' # 1) 10 biological replicats per group
#' # 2) 25 biological replicats per group
#' # 3) 50 biological replicats per group
#' # 4) 100 biological replicats per group
#' list_samples_per_group <- c(10, 25, 50, 100)
#'
#' # save the simulation results under each sample size
#' multiple_sample_sizes <- list()
#' for(i in seq_along(list_samples_per_group)){
#'     # run simulation for each sample size
#'     simulated_datasets <- simulateDataset(data = OV_SRM_train,
#'                                           annotation = OV_SRM_train_annotation,
#'                                           num_simulations = 10, # simulate 10 times
#'                                           expected_FC = "data",
#'                                           list_diff_proteins =  NULL,
#'                                           select_simulated_proteins = "proportion",
#'                                           protein_proportion = 1.0,
#'                                           protein_number = 1000,
#'                                           samples_per_group = list_samples_per_group[i],
#'                                           simulate_valid = FALSE,
#'                                           valid_samples_per_group = 50)
#'
#'     # run classification performance estimation for each sample size
#'     res <- designSampleSizeClassification(simulations = simulated_datasets,
#'                                           parallel = TRUE)
#'
#'     # save results
#'     multiple_sample_sizes[[i]] <- res
#' }
#'
#' ## make the plots and save them to disk
#' designSampleSizeClassificationPlots(data = multiple_sample_sizes, save.pdf = T)
#'
#'## make accuracy plot print in the Plots panes
#' designSampleSizeClassificationPlots(data = multiple_sample_sizes, =predictive_accuracy_plot = T)
#' @importFrom reshape2 melt
#' @import ggplot2
#' @importFrom gridExtra grid.arrange
#' @importFrom utils sessionInfo read.table write.table
#' @importFrom grDevices pdf dev.off
#' @export
#'
designSampleSizeClassificationPlots <- function(data,
                                                optimal_threshold = 0.001,
                                                num_important_proteins_show = 10,
                                                protein_importance_plot = TRUE,
                                                predictive_accuracy_plot = TRUE,
                                                save.pdf = F, ...) {

    ###############################################################################
    ## log file
    ## save process output in each step
    dots <- list(...)
    session <- dots$session
    func <- as.list(sys.call())[[1]]
    conn <- .find_log(...)
    ################################################################################
    
    res <- .catch_faults({
        if('ssclassification' %in% class(data)){
            f_imp <- .format_df(dat = data$mean_feature_importance,
                                sample = unique(data$num_samples),
                                top_n = num_important_proteins_show)
            
            acc_tbl <- .format_df(dat = data$predictive_accuracy,
                                  sample = unique(data$num_samples))
            
        }else{
            f_imp <- do.call('rbind', lapply(data, function(x){
                .format_df(dat = x$mean_feature_importance,
                           sample = unique(x$num_samples),
                           top_n = num_important_proteins_show)
            }))
            
            acc_tbl <- do.call('rbind', lapply(data, function(x){
                .format_df(dat = x$predictive_accuracy,
                           sample = unique(x$num_samples))
            }))
            
        }
        
        names(f_imp) <- c('frequency', 'protein', 'sample')
        names(acc_tbl) <- c('acc', 'simulation', 'sample')
        ylim_imp <- c(0,length(unique(acc_tbl$simulation)))
        
        if(length(unique(acc_tbl$sample))>1){
            opt_obj <- .identify_optimal(df = acc_tbl, cutoff = optimal_threshold)
            opt_val <- opt_obj$opt
            acc_plot <- .plot_acc(df = opt_obj$df, y_lim = opt_obj$y_lim,
                                  optimal_ss = opt_val)
        }else{
            opt_val <- unique(acc_tbl$sample)
            y_lim <- c(min(acc_tbl$acc)-0.1, 1)
            acc_tbl$fill_col <- 'red'
            acc_tbl$mean_acc <- mean(acc_tbl$acc, na.rm = T)
            acc_plot <- .plot_acc(df = acc_tbl, y_lim = y_lim,
                                  optimal_ss = opt_val)
        }
        
        p <- NULL
        if(save.pdf | (protein_importance_plot && predictive_accuracy_plot)){
            
            if(predictive_accuracy_plot){
                file <- sprintf("Accuracy_plot_%s.pdf",
                                format(Sys.time(), "%Y%m%d%H%M%S"))
                .status(detail = 'Plotting Accuracy plots', log = conn$con, ...)
                pdf(file = file)
                print(acc_plot)
                dev.off()
                .status(detail = sprintf("Accuracy Plot stored at %s", file),
                        log = conn$con, ...)
            }
            
            if(protein_importance_plot){
                file <- sprintf("Protein_importance_plot_%s.pdf",
                                format(Sys.time(), "%Y%m%d%H%M%S"))
                plots <- list()
                pdf(file = file)
                for(i in unique(f_imp$sample)){
                    df <- subset(f_imp, sample == i)
                    print(.plot_imp(df = df, samp = i,
                                    ylim = ylim_imp,
                                    x.axis.size = 6, y.axis.size = 6,
                                    margin = 0.5))
                }
                dev.off()
                .status(detail = sprintf("Protein Importance plots stored at %s",
                                         file), log = conn$con, ...)
            }
            
            p <- NULL
        } else if (predictive_accuracy_plot){
            p <- acc_plot
        } else if(protein_importance_plot && !is.null(session)){
            opt_val <- ifelse(is.null(dots$samp), opt_val, dots$samp)
            p <- .plot_imp(df = f_imp, samp = opt_val, ylim = ylim_imp,
                           facet = FALSE)
        } else if(protein_importance_plot && is.null(session)){
            p <- .plot_imp(df = f_imp, ylim = ylim_imp, facet = TRUE)
        }
        .status(detail = sprintf("Estimated optimal sample size is %s", opt_val),
                log = conn$con, ...)
        
        return(list('optimal_size'=opt_val, 'plot' = p))
    }, conn = conn, session = session)
    return(res)
}
