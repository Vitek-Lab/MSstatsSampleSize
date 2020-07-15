#' Mean-SD plot
#' @description Draw the plot for the mean protein abundance vs standard deviation
#' (square root of variance). The `lowess` function is used to fit the LOWESS 
#' smoother between mean protein abundance and standard deviation.
#'
#' @param data A list with mean protein abundance vector and standard deviation vector.
#' It should be the output of \code{\link{estimateVar}} function.
#' @param smoother_size Size of lowess smoother. Default is 1.
#' @param xlimUp The upper limit of x-axis for mean-SD plot. Default is 30.
#' @param ylimUp The upper limit of y-axis for mean-SD plot. Default is 3.
#' @param save.pdf A logical input, determines to save the plots as a pdf or not,
#' the pdf plot is saved in the current working directory, name of the created
#' file is displayed on the console and logged for easier access
#' @return \emph{meanSDplot} is the plot for the mean protein abundance (X-axis)
#' vs standard deviation (Y-axis).
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
meanSDplot <- function(data, smoother_size = 1, xlimUp = 30, ylimUp = 3, 
                       save.pdf = F, ...){
    
    ..density.. <-  freq <- x <- y <- NULL
    ###############################################################################
    ## log file
    ## save process output in each step
    dots <- list(...)
    session <- dots$session
    func <- as.list(sys.call())[[1]]
    if(is.null(dots$log_conn)){
        conn = mget("LOG_FILE", envir = .GlobalEnv,
                    ifnotfound = NA)
        if(is.na(conn)){
            rm(conn)
            conn <- .logGeneration()
        } else{
            conn <- .logGeneration(file = conn$LOG_FILE)   
        }
    }else{
        conn <- dots$log_conn
    }
    
    ###############################################################################
    
    plotdata <- data.frame(mean=data$promean, sd=data$prosd)
    plot.lowess <- lowess(cbind(plotdata$mean, plotdata$sd))
    plot.lowess <- data.frame(x = plot.lowess$x, y = plot.lowess$y)
    .status(detail = "Lowess data calculated", log = conn$con, func =func, ...)
    
    meansdplot <-  ggplot(data = plotdata, aes(mean, sd)) +
        stat_density2d(aes(fill = ..density..^0.25), geom = "tile", 
                       contour = FALSE, n = xlimUp*10) +
        scale_fill_continuous(low = "white", high = "#0072B2")+
        geom_point(alpha = 0.02, shape = 20)+
        geom_line(data = plot.lowess, aes(x, y), color="orange", size = smoother_size) +
        labs(x = "Mean protein abundance",
             y = "Standard deviation") +
        scale_y_continuous(expand = c(0,0), limits = c(0,ylimUp)) +
        scale_x_continuous(expand = c(0,0), limits = c(0,xlimUp)) +
        theme_MSstats(legend.position = "none", download = T)
    
    if(save.pdf){
        file_name <- sprintf("meandSD_plot_%s.pdf", format(Sys.time(),'%H%M%s'))
        file <- file.path(getwd(), file_name)
        pdf(file, height = 5, width = 5)
        print(meansdplot)
        dev.off()
        .status(detail = sprintf("File saved at %s", file), log = conn$con,
                func = func, ...)
    }
    .status(detail = "MeanSD plot printed", log = conn$con, func = func, ...)
    
    return(meansdplot)
}