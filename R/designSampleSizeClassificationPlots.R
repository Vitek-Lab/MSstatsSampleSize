#############################################
## designSampleSizeClassificationPlots
#############################################
#' Visualization for sample size calculation in classification
#'
#' @description To illustrate the mean classification accuracy and protein importance under different sample sizes
#' through predictive accuracy plot and protein importance plot.
#'
#' @details This function visualizes for sample size calculation in classification.
#' Mean predictive accuracy and mean protein importance under each sample size is from the input `data',
#' which is the output from function \code{\link{designSampleSizeClassification}}.
#'
#' To illustrate the mean predictive accuracy and protein importance under different sample sizes,
#' it generates two types of plots in pdf files as output: (1) The predictive accuracy plot,
#' The X-axis represents different sample sizes and y-axis represents the mean predictive accuracy.
#' The reported sample size per condition can be used to design future experiment
#'
#' (2) The protein importance plot includes multiple subplots.
#' The number of subplots is equal to `list_samples_per_group'.
#' Each subplot shows the top `num_important_proteins_show` most important proteins under each sample size.
#' The Y-axis of each subplot is the protein name and X-axis is the mean protein importance under the sample size.
#'
#' @param data A list of outputs from function \code{\link{designSampleSizeClassification}}. Each element represents the results under a specific sample size.
#' The input should include at least two simulation results with different sample sizes.
#' @param list_samples_per_group A vector includes the different sample sizes simulated. This is required.
#' The number of simulated sample sizes in the input `data' should be equal to the length of list_samples_per_group
#' @param num_important_proteins_show The number of proteins to show in protein importance plot.
#' @param protein_importance_plot TRUE(default) draws protein importance plot.
#' @param predictive_accuracy_plot TRUE(default) draws predictive accuracy plot.
#' @param x.axis.size Size of x-axis labeling in predictive accuracy plot and protein importance plot. Default is 10.
#' @param y.axis.size Size of y-axis labels in predictive accuracy plot and protein importance plot. Default is 10.
#' @param protein_importance_plot_width Width of the saved pdf file for protein importance plot. Default is 3.
#' @param protein_importance_plot_height Height of the saved pdf file for protein importance plot. Default is 3.
#' @param predictive_accuracy_plot_width Width of the saved pdf file for predictive accuracy plot. Default is 4.
#' @param predictive_accuracy_plot_height Height of the saved pdf file for predictive accuracy plot. Default is 4.
#' @param ylimUp_predictive_accuracy The upper limit of y-axis for predictive accuracy plot. Default is 1. The range should be 0 to 1.
#' @param ylimDown_predictive_accuracy The lower limit of y-axis for predictive accuracy plot. Default is 0.0. The range should be 0 to 1.
#' @param address the name of folder that will store the results. Default folder is the current working directory.
#' The other assigned folder has to be existed under the current working directory.
#' An output pdf file is automatically created with the default name of `PredictiveAccuracyPlot.pdf' and `ProteinImportancePlot.pdf'.
#' The command address can help to specify where to store the file as well as how to modify the beginning of the file name.
#' If address=FALSE, plot will be not saved as pdf file but showed in window.
#'
#' @return predictive accuracy plot is the mean predictive accuracy under different sample sizes.
#' The X-axis represents different sample sizes and y-axis represents the mean predictive accuracy.
#' @return protein importance plot includes multiple subplots. The number of subplots is equal to `list_samples_per_group'.
#' Each subplot shows the top `num_important_proteins_show' most important proteins under each sample size.
#' The Y-axis of each subplot is the protein name and X-axis is the mean protein importance under the sample size.
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
#' ## make the plots
#' designSampleSizeClassificationPlots(data = multiple_sample_sizes,
#'                                     list_samples_per_group = list_samples_per_group)
#'
#' @importFrom reshape2 melt
#' @import ggplot2
#' @importFrom gridExtra grid.arrange
#' @importFrom utils sessionInfo read.table write.table
#' @importFrom grDevices pdf dev.off
#' @export
#'
designSampleSizeClassificationPlots <- function(data,
                                                list_samples_per_group,
                                                num_important_proteins_show = 10,
                                                protein_importance_plot = TRUE,
                                                predictive_accuracy_plot = TRUE,
                                                x.axis.size = 10,
                                                y.axis.size = 10,
                                                protein_importance_plot_width = 3,
                                                protein_importance_plot_height = 3,
                                                predictive_accuracy_plot_width = 4,
                                                predictive_accuracy_plot_height = 4,
                                                ylimUp_predictive_accuracy = 1,
                                                ylimDown_predictive_accuracy = 0.0,
                                                address = "") {

    prots <- freq <- NULL

    ###############################################################################
    ## log file
    ## save process output in each step
    loginfo <- .logGeneration()
    finalfile <- loginfo$finalfile
    processout <- loginfo$processout

    processout <- rbind(processout,
                        as.matrix(c(" ", " ", "MSstatsSampleSize - designSampleSizeClassificationPlot", " "), ncol=1))

    ################################################################################
    ## need to simulate at least two different sample sizes
    num_sample_size <- length(data)

    if ( num_sample_size < 2 ) {

        processout <- rbind(processout,
                            "The required input should include at least two simulation results with different sample sizes from simulateDataset and designSampleSizeClassifications. - stop")
        write.table(processout, file=finalfile, row.names=FALSE)

        stop("The required input should include at least two simulation results with different sample sizes from simulateDataset anddesignSampleSizeClassifications.")
    }

    ## the number of simulation in data should be equal to the number of list_samples_per_group
    if ( num_sample_size != length(list_samples_per_group) ) {

        processout <- rbind(processout,
                            "The number of simulation in the input of designSampleSizeClassificationPlot should be equal to the number of list_samples_per_group. - stop")
        write.table(processout, file=finalfile, row.names=FALSE)

        stop("The number of simulation in the input of designSampleSizeClassificationPlot should be equal to the number of list_samples_per_group.")
    }


    ## !! assumption : multiple simulations came from only different sample sizes
    ## !! If there are multiple simulation from different protein numbers, the function should be changed.

    ################################################################################
    ## 1. make the plot for protein importance

    sample_size <- NULL
    mean_PA <- NULL
    FI_plots <- list()

    ## if there is list protein numbers. we need rows and columns.

    for(i in seq_len(num_sample_size)){

        sample_size <- c(sample_size, list_samples_per_group[i])
        mean_PA <- c(mean_PA, data[[i]]$mean_predictive_accuracy)
        num_simulations <- ncol(data[[i]]$feature_importance)

        FI <- data[[i]]$mean_feature_importance
        # select the top most important proteins
        FI <- FI[order(FI, decreasing=TRUE)[num_important_proteins_show:1]]
        FI <- data.frame(freq=FI, prots=names(FI), row.names = NULL)
        FI$prots <- factor(FI$prots, levels = FI$prots)
        FI_plots[[i]] <- ggplot(data=FI, aes(x = prots, y = freq))+
            geom_bar(stat="identity") +
            coord_flip()+
            labs(title = paste(list_samples_per_group[i], "samples/group"), size=10,
                 x = "",
                 y = "") +
            ylim(0, num_simulations)+
            theme(
                panel.background = element_rect(fill = 'white', colour = "black"),
                panel.grid.major = element_line(colour = 'gray95'),
                panel.grid.minor = element_blank(),
                strip.background = element_rect(fill = 'gray95'),
                strip.text.x = element_text(colour = c("#00B0F6"), size=14),
                axis.text.x = element_text(size = x.axis.size, colour="black"),
                axis.text.y = element_text(size = y.axis.size, colour="black"),
                axis.ticks = element_line(colour = "black"),
                axis.title.x = element_text(size = x.axis.size+5, vjust= -0.4),
                axis.title.y = element_text(size = y.axis.size+5, vjust=0.3),
                title = element_text(size = x.axis.size+2, vjust=1.5))
    }


    ## size of width : * n
    height <- protein_importance_plot_height
    width <- protein_importance_plot_width * num_sample_size

    ## print out Protein Importance Plot
    if(protein_importance_plot){
        if (address != FALSE) {
            allfiles <- list.files()

            num <- 0
            plotfilenaming <- paste0(address, "ProteinImportancePlot")
            plotfinalfile <- paste0(address, "ProteinImportancePlot.pdf")

            while (is.element(plotfinalfile, allfiles)) {
                num <- num + 1
                plotfinalfile <- paste0(paste(plotfilenaming, num, sep="-"), ".pdf")
            }

            pdf(plotfinalfile, width=width, height=height)
        }

        do.call(grid.arrange, c(FI_plots, list(ncol=length(list_samples_per_group))))


        if (address != FALSE) {
            dev.off()
        }

        processout <- rbind(processout, as.matrix(c(" Drew Protein Importance Plot."), ncol=1))
        write.table(processout, file=finalfile, row.names=FALSE)
        message(" Drew Protein Importance Plot.")

    }


    ##############################################################################
    ## 2. make the plot for predictive accuracy

    ## get the mean accuracy
    plotdata1 <- data.frame(meanPA = mean_PA,
                            samplesize = sample_size)

    if(predictive_accuracy_plot){
        if (address != FALSE) {
            allfiles <- list.files()

            num <- 0
            plotfilenaming <- paste0(address, "PredictiveAccuracyPlot")
            plotfinalfile <- paste0(address, "PredictiveAccuracyPlot.pdf")

            while (is.element(plotfinalfile, allfiles)) {
                num <- num + 1
                plotfinalfile <- paste0(paste(plotfilenaming, num, sep="-"), ".pdf")
            }

            pdf(plotfinalfile, width=predictive_accuracy_plot_width, height=predictive_accuracy_plot_height)
        }

        ## need to update for multiple protein numbers

        p1 <- ggplot(data = plotdata1, aes(x= sample_size, y= mean_PA)) +
            geom_point() +
            geom_line() +
            scale_x_continuous(breaks = list_samples_per_group) + ## specify the defined sample size
            scale_y_continuous(limits = c(ylimDown_predictive_accuracy, ylimUp_predictive_accuracy)) +
            labs(x="Pre-defined sample size", y='Mean accuracy') +
            theme(
                panel.background = element_rect(fill = 'white', colour = "black"),
                panel.grid.major = element_line(colour = 'gray95'),
                panel.grid.minor = element_blank(),
                strip.background=element_rect(fill = 'gray95'),
                strip.text.x = element_text(colour = c("#00B0F6"), size=14),
                axis.text.x = element_text(size = x.axis.size, colour="black"),
                axis.text.y = element_text(size = y.axis.size, colour="black"),
                axis.ticks = element_line(colour = "black"),
                axis.title.x = element_text(size = x.axis.size+4, vjust=-0.4),
                axis.title.y = element_text(size = y.axis.size+4, vjust=0.3),
                title = element_text(size = x.axis.size+3, vjust=1.5))

        print(p1)

        if (address != FALSE) {
            dev.off()
        }

        processout <- rbind(processout, as.matrix(c(" Drew Predictive Accuracy Plot."), ncol=1))
        write.table(processout, file=finalfile, row.names=FALSE)
        message(" Drew Predictive Accuracy Plot.")

    }

}

