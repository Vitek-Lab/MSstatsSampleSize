#' @export
MSstatsSSE_shiny <- function() {
  appDir <- system.file("sample_size_estimator", package = "MSstatsSampleSize")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `MSstatsSampleSize`.", call. = FALSE)
  }
  
  shiny::runApp(appDir, display.mode = "normal")
}