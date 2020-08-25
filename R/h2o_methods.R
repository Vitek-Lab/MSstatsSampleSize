#' @title H2o configuration
#' @description Get configuration details of the h2o instance to start, the details
#' can be specified in the Rprofile.site file in the R installation path
#' @importFrom parallel detectCores
#' @keywords internal
#' @return A named list containing required h2o configuration details, if none 
#' provided, defaults are used
h2o_config <- function(...){
  options("h2o.use.data.table" = T)
  options('h2o.fwrite' = T)
  config <- list()
  func <- as.list(sys.call())[[1]]
  config$threads <- max(parallel::detectCores()-1,
                        1)
  
  mem <- as.numeric(system("awk '/MemFree/ {print $2}' /proc/meminfo", 
                           intern=TRUE))/(1024*1024)
  if(mem<2){
    stop("CALL_",func,"_Atleast 2Gb of RAM required to use h2o")
  }
  
  if(mem>2 && mem <4){
    warning("CALL_",func,"_Atleast 4Gb of RAM recommended to be used with h2o")
  }

  config$max_mem <- paste0(format(mem,digits = 0),"G")
  log_dir <- file.path(getwd(),"h2o_logs")
  dir.create(log_dir, showWarnings = F)
  config$log_dir <- log_dir
  config$log_level <- "DEBUG"
  return(config) 
}



#' @title Classify using h2o classification algorithms
#' @description This function enables the shiny application to apply various 
#' algorithms as defined by the user
#' @param n_samp A comma separated chracter vector with number of sample sizes 
#' @param sim_data A list object containing simulated data
#' @param classifier Classifier to use as selected by the user can use one of the
#' following `Logistic Regression`, `Random Forest`, `Support Vector Machine`,
#' `Naive Bayes`, `Neural Network`
#' @param stopping_metric 
#' @param seed
#' @param nfolds
#' @param fold_assignment
#' @param iters
#' @param alpha
#' @param family
#' @param solver
#' @param link
#' @param min_sdev
#' @param laplace
#' @param eps
#' @param session 
ss_classify_h2o <- function(n_samp, sim_data, classifier, stopping_metric = "AUTO",
                            seed = -1, nfolds = 0, fold_assignment = "AUTO", iters = 200,
                            alpha = 0, family, solver, link, min_sdev, laplace, eps,
                            session = NULL){
  
  samp <- unlist(strsplit(n_samp,","))
  config <- h2o_config()
  h2o::h2o.init(nthreads = config$threads, max_mem_size = config$max_mem,
                log_dir = config$log_dir, log_level = config$log_level)
  max_val <- 0
  iter <- 0
  modelz <- list()
  f_imp <- data.table()
  
  for(i in samp){
    
    valid_x <- sim_data[[i]]$valid_X
    valid_y <- sim_data[[i]]$valid_Y
    valid <- h2o::as.h2o(cbind("condition" = valid_y, valid_x),
                         destination_frame = "valid")
    
    train_x_list <- sim_data[[i]]$simulation_train_Xs 
    train_y_list = sim_data[[i]]$simulation_train_Ys
    for(index in seq_along(train_x_list)){
      if(max_val == 0){
        max_val <- length(train_x_list) * length(samp)
      }
      iter = iter + 1/max_val
      
      status(detail = sprintf("Classifying Sample Size %s of %s, Simulation %s of %s",
                              which(samp == i), length(samp), index, length(train_x_list)),
             session = session, value = iter)
      
      x <- train_x_list[[index]]
      y <- train_y_list[[index]]
      train <- h2o::as.h2o(cbind("condition" = y, x), destination_frame = "train")
      
      if(classifier == "rf"){
        model <- h2o::h2o.randomForest(y = 1, training_frame = train,
                                       stopping_metric = stopping_metric, seed = seed, 
                                       balance_classes = FALSE, nfolds = nfolds,
                                       fold_assignment = fold_assignment,
                                       build_tree_one_node = T)
        var_imp <- h2o::h2o.varimp(model)
        sel_imp <- var_imp$variable[1:10]
        h2o::h2o.rm(train)
        train <- h2o::as.h2o(cbind("condition"=y,x[,sel_imp]), destination_frame = "train")
        model <- h2o::h2o.randomForest(y = 1, training_frame = train,
                                       validation_frame = valid,
                                       stopping_metric = stopping_metric, seed = seed, 
                                       balance_classes = FALSE, nfolds = nfolds,
                                       fold_assignment = fold_assignment,
                                       build_tree_one_node = T)
        
      } else if (classifier == "nnet"){
        l1 = 0
        l2 = 0
        rho = 0.99
        epochs = 10
        hidden = c(250,250)
        activation = "Rectifier"
        model <- h2o::h2o.deeplearning(y = 1, 
                                       training_frame = train,
                                       l1= l1, l2=l2,
                                       activation = activation,
                                       hidden = hidden,
                                       epochs = epochs)
        found_imp <- h2o::h2o.varimp(model)
        sel_imp <- found_imp$variable[1:10]
        train <- h2o::as.h2o(cbind("condition"=y,x[,sel_imp]), destination_frame = "train")
        model <- h2o::h2o.deeplearning(y = 1, training_frame = train,
                                       validation_frame = valid, l1 = l1, l2 = l2,
                                       activation = activation, hidden = hidden,
                                       epochs = epochs)
        
      } else if (classifier == "svmLinear"){
        model <- h2o::h2o.psvm(y = 1, training_frame = train, 
                               max_iterations = iters,
                               seed = seed, validation_frame = valid,
                               disable_training_metrics = F)
      } else if (classifier == "logreg"){
        if(length(unique(as.vector(train[,1]))) > 2){
          family <- "multinomial"
        }
        model <- h2o::h2o.glm(y = 1, training_frame = train,
                              seed = seed, family = family, alpha = alpha, 
                              nfolds = nfolds, solver = solver,
                              link = link)
        
        var_imp <- h2o::h2o.varimp(model)
        sel_imp <- var_imp$variable[1:10]
        h2o::h2o.rm(train)
        train <- h2o::as.h2o(cbind("condition"=y,x[,sel_imp]), destination_frame = "train")
        model <- h2o::h2o.glm(y = 1, training_frame = train, seed = seed,
                              family = family, alpha = alpha, nfolds = nfolds, 
                              solver = solver, link = link, validation_frame = valid)
      } else if (classifier == "naive_bayes"){
        model <- h2o::h2o.naiveBayes(y = 1, training_frame = train,
                                     ignore_const_cols = TRUE,
                                     nfolds = nfolds,
                                     fold_assignment = fold_assignment,
                                     seed = seed, laplace = laplace,
                                     min_sdev = min_sdev,
                                     eps_sdev = eps)
      } else{
        stop("Not defined")
      }
      perf <- h2o::h2o.performance(model = model, newdata = valid)
      cm <- perf@metrics$cm$table[1:length(unique(y)),
                                  1:length(unique(y))]
      cm <- as.matrix(sapply(cm, as.numeric))
      #accuracy of current simulation
      acc <- sum(diag(cm))/sum(cm)
      
      name_val <- sprintf("Sample%s %s", i,  names(train_x_list)[index])
      var_imp <-  h2o::h2o.varimp(model)
      rm(train, perf, cm, model)
      
      modelz[[name_val]] <- list("acc" = acc, "var_imp" = var_imp)
    }
  }
  rm(valid, valid_x, valid_y, train_x_list, train_y_list)
  gc()
  h2o:::.h2o.garbageCollect()
  h2o::h2o.shutdown(prompt = F)
  
  return(list("models" = modelz))
}



.classification_performance_h2o <- function(index, classifier, train_x_list,
                                            train_y_list, valid_x, valid_y, 
                                            top_K, ...){
  
  train <- train_x_list[[index]]
  predictors <- names(train)
  train <- cbind(train, "condition" = train_y_list[[index]])
  train <- as.h2o(train, destination_frame = 'training_frame')
  train["condition"] <- as.factor(train["condition"])
  
  model <- .classification_model_h2o(x = predictors, y = 'condition', 
                                     training_frame = train,
                                     classifier = classifier, ...)
  
  fimp <- h2o.varimp(model)
  predictors <- fimp[seq_len(top_K), 'variable']
  
  pred.model <- .classification_model_h2o(x = predictors, y = 'condition',
                                          training_frame = train,
                                          classifier = classifier, ...)
  
  valid <- valid_x[, (colnames(valid_x) %in% predictors)]
  valid <- cbind(valid, "condition" = valid_y)
  valid <- as.h2o(valid)
  valid["condition"] <- as.factor(valid["condition"])
  
  perf <- h2o::h2o.performance(model = model, newdata = valid)
  cm <- perf@metrics$cm$table[1:length(unique(valid_y)),
                              1:length(unique(valid_y))]
  cm <- as.matrix(sapply(cm, as.numeric))
  #accuracy of current simulation
  acc <- sum(diag(cm))/sum(cm)
  
  run <- list(acc = acc, f_imp = predictors)
  return(run)
}




.classification_model_h2o <- function(x, y, training_frame, classifier, ...){
  models <- c("h2o::h2o.randomForest", "h2o::h2o.deeplearning",
              "h2o::h2o.naiveBayes", "h2o::h2o.glm", "h2o::h2o.psvm")
  names(models) <- c("rf", "nnet", "naive_bayes", "logreg", "svmLinear")
  
  call <- match.call(expand.dots = T)
  call[[5]] <- NULL
  #call[[6]] <- NULL
  func <- unname(models[classifier])
  call[[1]] <- func
  str_func <- deparse(call)
  str_func <- gsub('"',"", str_func)
  str_func <- gsub("= train", "= training_frame", str_func)
  str_func <- gsub("= predictors", "= x", str_func)
  str_func <- gsub("condition", "'condition'", str_func)
  call <- parse(text = str_func)
  output <- eval(call)
  
  return(output)
}
# res3 <- list()
# 
# 
# for(i in seq_along(sim1$simulation_train_Xs)){
# 
#   res3[[i]] <- .classification_performance_h2o(index = i, classifier = 'naive_bayes',
#                                          train_x_list = sim1$simulation_train_Xs,
#                                          train_y_list = sim1$simulation_train_Ys,
#                                          valid_x = sim1$valid_X, valid_y = sim1$valid_Y,
#                                          top_K = 10)
# 
# }
# 
# 
