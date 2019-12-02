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
designSampleSizePCAplot <- function(simulations,
                                    which.PCA = "all",
                                    x.axis.size = 10,
                                    y.axis.size = 10,
                                    dot.size = 3,
                                    legend.size = 7,
                                    width = 6,
                                    height = 5,
                                    address = "") {

    ###############################################################################
    ## log file
    ## save process output in each step
    loginfo <- .logGeneration()
    finalfile <- loginfo$finalfile
    processout <- loginfo$processout

    processout <- rbind(processout,
                        as.matrix(c(" ", " ", "MSstatsSampleSize - designSampleSizePCAplot", " "), ncol=1))

    ## parameter checking: which.PCA
    if (!(length(which.PCA) == 1 &
        (which.PCA %in% c("all", "data") |
         (startsWith(which.PCA, "simulation"))))) {

        processout <- rbind(processout, c("ERROR : which.PCA should be one of \"all\", \"data\" or \"simulation + index of simulation\" (such as \"simulation1\")."))
        write.table(processout, file=finalfile, row.names=FALSE)

        stop("which.PCA should be one of \"all\", \"data\" or \"simulation + index of simulation\" (such as \"simulation1\"). \n")

    }

    processout <- rbind(processout, c(paste0("which.PCA = ", which.PCA)))
    write.table(processout, file=finalfile, row.names=FALSE)

    ###############################################################################
    ## number of simulations
    iter <- length(simulations$simulation_train_Xs)

    if (address != FALSE) {
        allfiles <- list.files()

        num <- 0
        plotfilenaming <- paste0(address, "PCAPlot")
        plotfinalfile <- paste0(address, "PCAPlot.pdf")

        while (is.element(plotfinalfile, allfiles)) {
            num <- num + 1
            plotfinalfile <- paste0(paste(plotfilenaming, num, sep="-"), ".pdf")
        }

        pdf(plotfinalfile, width=width, height=height)
    }

    #########################################
    # draw the first PCA plot on the input preliminary dataset
    input <- simulations$input_X
    input.group <- simulations$input_Y
    input.group <- factor(input.group, levels=sort(levels(input.group)))

    if(which.PCA == "all" | which.PCA == "allonly"){
        ## PCA
        result.pca <- prcomp(input, scale. = TRUE)
        summary.pca <- summary(result.pca)
        # $x : the values of each sample in terms of the principal components

        # extract PC1 and PC2
        important.pc <- result.pca$x[, 1:2]
        pc.result <- data.frame(important.pc, group=input.group)

        # get explained var. : 'Proportion of variance'
        exp.var <- summary.pca$importance
        exp.var <- format(exp.var[2, 1:2]*100, digits=3)

        ## PC1 and PC2 are not scaled like biplot
        pcaplot <- ggplot(data = pc.result,
                          aes_string(x = 'PC1', y = 'PC2', color = 'group')) +
            geom_point(size = dot.size) +
            stat_ellipse() +
            labs(title = "Input dataset",
                 x = paste0("PC1 (", exp.var[1], "% explained var.)"),
                 y = paste0("PC2 (", exp.var[2], "% explained var.)")) +
            theme(
                panel.background = element_rect(fill = 'white', colour = "black"),
                panel.grid.major = element_line(colour = 'gray95'),
                panel.grid.minor = element_blank(),
                strip.background=element_rect(fill = 'gray95'),
                strip.text.x=element_text(colour = c("#00B0F6"), size=14),
                axis.text.x=element_text(size=x.axis.size, colour="black"),
                axis.text.y=element_text(size=y.axis.size, colour="black"),
                axis.ticks=element_line(colour="black"),
                axis.title.x=element_text(size=x.axis.size+5, vjust=-0.4),
                axis.title.y=element_text(size=y.axis.size+5, vjust=0.3),
                title=element_text(size=x.axis.size+8, vjust=1.5),
                legend.key = element_rect(fill='white', colour='white'),
                legend.position="right",
                legend.text=element_text(size=legend.size),
                legend.title = element_blank())

        print(pcaplot)

        processout <- rbind(processout, as.matrix(c(" Drew the PCA plot for input preliminary dataset "), ncol=1))
        write.table(processout, file=finalfile, row.names=FALSE)
        message(" Drew the PCA plot for input preliminary dataset ")

    }


    #########################################
    # draw the PCA plot for simulation data
    if(which.PCA == "all" | startsWith(which.PCA, "simulation")){
        if(which.PCA == "all") {
            index <- seq_len(iter) # draw plots for all the simulations

        } else{
            which.PCA <- gsub("simulation", "", which.PCA)
            index <- as.integer(which.PCA)

            if(index > iter){
                processout <- rbind(processout, c(paste0("ERROR: which.PCA should be one of \"all\", \"data\" or \"simulation + index of simulation\".
                                                         index must be no bigger than ", iter)))
                write.table(processout, file=finalfile, row.names=FALSE)

                stop("which.PCA should be one of \"all\", \"data\" or \"simulation + index of simulation\". index must be no bigger than ", iter)
            }

        }

        for(i in seq_along(index)) {

            # prepare the data for PCA analysis
            input <- simulations$simulation_train_Xs[[index[i]]]
            input.group <- simulations$simulation_train_Ys[[index[i]]]
            input.group <- factor(input.group, levels=sort(levels(input.group)))

            ## PCA
            result.pca <- prcomp(input, scale. = TRUE)
            ## There are two different scaling ( in the prcomp function, in the autoplot function)
            ## 1. http://www.sthda.com/english/wiki/ggfortify-extension-to-ggplot2-to-handle-some-popular-packages-r-software-and-data-visualization
            ## 2. https://cran.r-project.org/web/packages/ggfortify/vignettes/plot_pca.html
            summary.pca <- summary(result.pca)
            # $x : the values of each sample in terms of the principal components

            # extract PC1 and PC2
            important.pc <- result.pca$x[, 1:2]

            pc.result <- data.frame(important.pc, group=input.group)

            # get explained var. : 'Proportion of variance'
            exp.var <- summary.pca$importance
            exp.var <- format(exp.var[2, 1:2]*100, digits=3)

            ## PC1 and PC2 are not scaled like biplot
            pcaplot <- ggplot(data = pc.result,
                              aes_string(x = 'PC1', y = 'PC2', color = 'group')) +
                geom_point(size = dot.size) +
                stat_ellipse() +
                labs(title = paste0("Simulation:", index[i]),
                     x = paste0("PC1 (", exp.var[1], "% explained var.)"),
                     y = paste0("PC2 (", exp.var[2], "% explained var.)")) +
                theme(
                    panel.background = element_rect(fill = 'white', colour = "black"),
                    panel.grid.major = element_line(colour = 'gray95'),
                    panel.grid.minor = element_blank(),
                    strip.background=element_rect(fill = 'gray95'),
                    strip.text.x=element_text(colour = c("#00B0F6"), size=14),
                    axis.text.x=element_text(size=x.axis.size, colour="black"),
                    axis.text.y=element_text(size=y.axis.size, colour="black"),
                    axis.ticks=element_line(colour="black"),
                    axis.title.x=element_text(size=x.axis.size+5, vjust=-0.4),
                    axis.title.y=element_text(size=y.axis.size+5, vjust=0.3),
                    title=element_text(size=x.axis.size+8, vjust=1.5),
                    legend.key = element_rect(fill='white', colour='white'),
                    legend.position="right",
                    legend.text=element_text(size=legend.size),
                    legend.title = element_blank())

            print(pcaplot)

            processout <- rbind(processout, as.matrix(c(paste0(" Drew the PCA plot for simulation ", index[i])), ncol=1))
            write.table(processout, file=finalfile, row.names=FALSE)
            message(paste0(" Drew the PCA plot for simulation ", index[i]))
        }

    }

    if (address != FALSE) {
        dev.off()
    }

}
