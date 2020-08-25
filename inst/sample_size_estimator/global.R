suppressMessages({
  library(shiny)
  library(shinydashboard)
  library(shinyWidgets)
  library(shinycssloaders)
  library(ggplot2)
  library(data.table)
  library(dplyr)
  library(MSstatsSampleSize)
  #### Caret package dependencies ####
  library(e1071)
  library(kernlab)
  library(naivebayes)
  library(randomForest)
})

wd <- getwd()
src_dir <- file.path(wd,"scripts")
req_files <- list.files(src_dir, full.names = T)
  
for (i in req_files) {
  tryCatch({
    source(i)
  }, error = function(e){
    print(conditionMessage(e))
    print(sprintf("Error in sourcing : %s", i))
  })
  
}
conn <<- .logGeneration()
message("Log File saved at ", conn$file)

#### GLOBAL VARS ####
FORMATS_LIST <- list("Protein-level quantification" = "standard", 
                     "Example from MSstatsSampleSize" = "examples")
FORMATS <- c("examples", "standard")
EXTENSTIONS <- c("text/csv",
                 "text/comma-separated-values,text/plain",
                 ".csv", "text/tab-separated-values", ".tsv")

MODELS <- c('rf','nnet','svmLinear','logreg','naive_bayes')
names(MODELS) <- c("Random Forest", "Neural Network",
                   "Support Vector Machines with Linear Kernel",
                   "Logistic Regression", "Naive Bayes")

STOPPING_METRIC <- c("AUTO", "deviance", "logloss", "MSE", "RMSE", "MAE", "RMSLE",
                     "AUC", "lift_top_group", "misclassification", "AUCPR",
                     "mean_per_class_error")

FOLD_ASSIGNMENT <- c("AUTO", "Random", "Modulo", "Stratified")


FAMILY <-  c("gaussian", "binomial", "quasibinomial", "ordinal", "multinomial",
           "poisson", "gamma", "tweedie", "negativebinomial")

SOLVER <- c("AUTO", "IRLSM", "L_BFGS", "COORDINATE_DESCENT_NAIVE", 
            "COORDINATE_DESCENT", "GRADIENT_DESCENT_LH", "GRADIENT_DESCENT_SQERR")

LINK <- c("family_default", "identity", "logit", "log", "inverse", "tweedie",
          "ologit")

B_GROUP <- ""
CURRMODEL <- ""
SIM_CHOICES <- 0

CSS_BUTTON <- "margin-top: 25px;
              display: inline-block;
              color: #fff;
              background-color: orange;
              border-color: black;
              font-size : 20px;"

CSS_BUTTON_REG <- "background-color: orange;
                  border-color: black;
                  color: #fff;"

CSS_BUTTON_RUN <- "background-color: orange;
                  border-color: black;
                  color: #fff;
                  margin-top:25px;"


#### Shiny functions ####
tuning_params <- function(...){
  dots <- list(...)
  params <- readRDS("params.rds")
  op <- list("enable" = F, "ids" = params$id)
  
  if("Parameter Tuning" %in% dots$checkbox){
    if("Use h2o Package" %in% dots$checkbox){
      show <- subset(params, classifier==dots$alg & use_h2o ==T)
      hide <- subset(params, classifier!=dots$alg | use_h2o == F)
    } else{
      show <- subset(params, classifier==dots$alg & use_h2o ==F)
      hide <- subset(params, classifier!=dots$alg | use_h2o ==T)
    }
    op <- list("enable" = T, "show_ids" = show$id, "hide_ids" = hide$id)
  }
  return(op)
}  



format_data <- function(format, count = NULL, annot = NULL, transform = NULL,
                        session = NULL, conn){
  
  shiny::validate(shiny::need(format %in% FORMATS, "Undefined Format"))
  if(format == "standard"){
    .status(detail = "Importing Protein Abundance file", value = 0.4,
           session = session, log = conn$con)
    #read abundance data from the file path provided
    wide <- read.csv(count$datapath, stringsAsFactors = FALSE)
    #No column names expected for the protein columns
    rownames(wide) <- wide[,1]
     wide[,1] <- NULL
    
    if(transform != "No Transform"){
      .status(detail = sprintf("Applying transformation of %s", transform),
              value = 0.45, session = session, log = conn$con)
      wide <- rapply(wide, f = get(transform), classes = c("numeric", "integer"),
                     how = "replace")
    }
    
    name <- count$name
    .status(detail = "Importing Annotations file", value = 0.5,
           session = session)
    #read annotations from the file path provided
    annot <- read.csv(annot$datapath, stringsAsFactors = FALSE)
    name <- count$name
    
  }else if(format == "examples"){
    .status(detail = "Importing Data from MSstatsSampleSize Package", value = 0.5,
           session = session, log = conn$con)
    #example data from the package
    wide <- MSstatsSampleSize::OV_SRM_train
    #examples data from the package
    annot <- MSstatsSampleSize::OV_SRM_train_annotation
    name <- "Ovarian Cancer SRM study"
  }else{
    stop("Not Defined")
  }
  #get summary about the bioreplicates and runs of the data
  data_summary <- .format_summary_table(data = annot)
  
  .status(detail = "Estimating Variance", value = 0.8, session = session,
         log = conn$con)
  
  var <- estimateVar(data = wide, annotation = annot, session = session,
                     log_conn = conn) 
  
  
  .status(detail = "Creating Summary Table", value = 0.9, session = session,
         log = conn$con)
  sum_table <- data.frame(Parameter =  c("Number of Proteins", "Number of Groups"),
                          Values = c(nrow(wide), length(annot$Condition)))
  .status(detail = "Creating Box Plots", value = 0.95, session = session,
         log = conn$con)
  .status(detail = sprintf("Dataset Name: %s", name), log = conn$con)
  return(list("wide_data" = wide, "annot_data" = annot,
              "n_prot" = nrow(wide), 
              "cond_sum_table" = data_summary,
              "n_group" = length(unique(annot$Condition)), "dataset_name" = name,
              "sum_table" = as.matrix(sum_table), "est_var" = var))
}
