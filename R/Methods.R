#' Simulate extended datasets for sample size estimation
#'
#' @param m number of samples per group to simulate
#' @param mu a matrix of mean abundance in each phenotype group and each protein
#' @param sigma a matrix of variance in each phenotype group and each protein
#' @return \emph{X} A simulated matrix with required sample size
#' @return \emph{Y} Group information corresponding with \emph{X}
#' @keywords internal
.sampleSimulation <- function(m, mu, sigma) {

    nproteins <- nrow(mu)
    ngroup <- ncol(mu)
    ## Determine the size of each phenotype group
    samplesize <- rep(m, ngroup)

    ## Simulate the data matrix
    sim_matrix <- matrix(rep(0, nproteins * sum(samplesize)), ncol=sum(samplesize))
    for (i in 1:nproteins) {
        abun <- NULL
        for (j in 1:ngroup) {
            abun <- c(abun, rnorm(samplesize[j], mu[i, j], sigma[i, j]))
        }
        sim_matrix[i, ] <- abun
    }

    sim_matrix <- t(sim_matrix)
    colnames(sim_matrix) <- rownames(mu)
    #Simulate the phenotype information
    group <- rep(colnames(mu), times=samplesize)

    index <- sample(length(group), length(group))
    sim_matrix <- sim_matrix[index, ]
    group <- group[index]

    return(list(X=sim_matrix,
                Y=as.factor(group)))
}



#' For each protein, impute the missing values based on the observed values
#'
#' @param data protein abundance data for one protein.
#' @return Imputed protein abundance data
#' @keywords internal
.random.imp <- function (data){
    missing <- is.na(data) # count missing values
    n.missing <- sum(missing)
    data.obs <- data[!missing] # keep the observed values
    imputed <- data
    # impute missing values by randomly selecting observed values
    imputed[missing] <- sample (data.obs, n.missing, replace=TRUE)
    return (imputed)
}

