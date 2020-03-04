#' Mean-SD plot
#' @description Draw the plot for the mean protein abundance vs standard deviation in each condition.
#' The `lowess' function is used to fit the LOWESS smoother between mean protein abundance and standard deviation.
#'
#' @param data A list with mean protein abundance matrix and standard deviation matrix.
#' It should be the output of \code{\link{estimateVar}} function.
#' @param x.axis.size Size of x-axis labeling in Mean-SD Plot. Default is 10.
#' @param y.axis.size Size of y-axis labels. Default is 10.
#' @param smoother_size Size of lowess smoother. Default is 1.
#' @param xlimUp The upper limit of x-axis for mean-SD plot. Default is 30.
#' @param ylimUp The upper limit of y-axis for mean-SD plot. Default is 3.
#' @param width Width of the saved pdf file. Default is 4.
#' @param height Height of the saved pdf file. Default is 4.
#' @param address The name of folder that will store the results. Default folder is the current working directory.
#' The other assigned folder has to be existed under the current working directory.
#' An output pdf file is automatically created with the default name of `MeanSDPlot.pdf'.
#' The command address can help to specify where to store the file as well as how to modify the beginning of the file name.
#' If address=FALSE, plot will be not saved as pdf file but showed in window.
#'
#' @return \emph{meanSDplot} is the plot for the mean protein abundance (X-axis) vs standard deviation (Y-axis) in each condition.
#' @author Ting Huang, Meena Choi, Olga Vitek
#' @examples
#' data(OV_SRM_train)
#' data(OV_SRM_train_annotation)
#'
#' variance_estimation <- estimateVar(data = OV_SRM_train,
#'                                    annotation = OV_SRM_train_annotation)
#'
#' meanSDplot(variance_estimation)
#'
#' @import ggplot2
#' @importFrom utils sessionInfo read.table write.table
#' @importFrom grDevices pdf dev.off
#' @importFrom stats lowess sd
#' @export
#'
meanSDplot <- function(data,
                       x.axis.size = 10,
                       y.axis.size = 10,
                       smoother_size = 1,
                       xlimUp = 30,
                       ylimUp = 3,
                       height = 4,
                       width = 4,
                       address = "") {


    ..density.. <-  freq <- x <- y <- NULL

    ###############################################################################
    ## log file
    ## save process output in each step
    loginfo <- .logGeneration()
    finalfile <- loginfo$finalfile
    processout <- loginfo$processout

    processout <- rbind(processout,
                        as.matrix(c(" ", " ", "MSstatsSampleSize - meanSDplot", " "), ncol=1))


    message("Drew the Mean-SD plot for simulation!")
    processout <- rbind(processout, as.matrix(c("Drew the Mean-SD plot for simulation."), ncol=1))
    write.table(processout, file=finalfile, row.names=FALSE)

    ###############################################################################

    if (address != FALSE) {
        allfiles <- list.files()

        num <- 0
        plotfilenaming <- paste0(address, "MeanSDPlot")
        plotfinalfile <- paste0(address, "MeanSDPlot.pdf")

        while (is.element(plotfinalfile, allfiles)) {
            num <- num + 1
            plotfinalfile <- paste0(paste(plotfilenaming, num, sep="-"), ".pdf")
        }

        pdf(plotfinalfile, width=width, height=height)
    }

    plotdata <- data.frame(mean=as.vector(data$mu), sd=as.vector(data$sigma))
    plot.lowess <- lowess(cbind(plotdata$mean,plotdata$sd))
    plot.lowess <- data.frame(x = plot.lowess$x, y = plot.lowess$y)

    meansdplot <-  ggplot(data = plotdata, aes(mean, sd)) +
        stat_density2d(aes(fill = ..density..^0.25), geom = "tile", contour = FALSE, n = xlimUp*10) +
        scale_fill_continuous(low = "white", high = "#0072B2")+
        geom_point(alpha = 0.02, shape = 20)+
        geom_line(data = plot.lowess, aes(x, y), color="orange", size = smoother_size) +
        labs(x = "Mean protein abundance per condition",
             y = "Standard deviation per condition") +
        scale_y_continuous(expand = c(0,0), limits = c(0,ylimUp)) +
        scale_x_continuous(expand = c(0,0), limits = c(0,xlimUp)) +
        theme(
            panel.background = element_rect(fill = 'white', colour = "black"),
            panel.grid.major = element_line(colour = 'black'),
            panel.grid.minor = element_blank(),
            strip.background = element_rect(fill = 'gray95'),
            strip.text.x = element_text(colour = c("#00B0F6"), size = 14),
            axis.text.x = element_text(size = x.axis.size, colour="black"),
            axis.text.y = element_text(size = y.axis.size, colour="black"),
            axis.ticks = element_line(colour = "black"),
            axis.title.x = element_text(size = x.axis.size+5, vjust = -0.4),
            axis.title.y = element_text(size = y.axis.size+5, vjust = 0.3),
            legend.position = "none")

    print(meansdplot)

    if (address != FALSE) {
        dev.off()
    }


}

