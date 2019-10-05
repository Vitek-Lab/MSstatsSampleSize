#' PCA plot for each simulation
#' @details This function draws PCA plot for each simulated dataset.
#' It outputs a pdf file where the number of pages is equal to the number of simulations in `simulations' (input for this function).
#' Each page presents a PCA plot under one simulation. x-axis of PCA plot is the first component and y-axis is the second component.
#' This function can be used to validate whether the simulated dataset looks as expected.
#'
#' @param simulations A list of simulated datasets. It should be the output of \code{\link{simulateDataset}} function.
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

    allfiles <- list.files()

    filenaming <- "MSstatsSampleSize-ProgressReport"

    if (length(grep(filenaming,allfiles)) == 0) {

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
                        as.matrix(c(" ", " ", "MSstatsSampleSize - designSampleSizePCAplot", " "), ncol=1))
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

    for(i in 1:iter) {

        # prepare the data for PCA analysis
        input <- simulations$simulation_train_Xs[[i]]
        input.group <- simulations$simulation_train_Ys[[i]]

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
            labs(title = paste0("Simulation:", i),
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

        processout <- rbind(processout, as.matrix(c(paste0(" Drew the PCA plot for simulation ", i)), ncol=1))
        write.table(processout, file=finalfile, row.names=FALSE)
        message(paste0("Drew the PCA plot for simulation ", i))
    }

    if (address != FALSE) {
        dev.off()
    }

}
