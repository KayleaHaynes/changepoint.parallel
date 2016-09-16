#' Parallel PELT SM2
#' @description Splits the data into chunks by sending the first point to the first core, the second point to the second core etc, and runs parallel changepoint detection on each chunk then merges the results and runs parallel detection using the detected changepoints in step 1 as the candidate changepoint locations.
#' @param data A vector of data-points within which you with to find changepoints.
#' @param penalty The value of the penalty.
#' @param pruning If true PELT is used, if false Optimal Partitioning is used.
#' @param sum.stat This can be ``norm.sum" or ``exp.sum".
#' @param cost This can be ``norm.mean.cost", ``norm.var.cost", ``norm.meanvar.cost" or ``exp.cost".
#' @param ncores Number of cores to use.
#'
#' @return The detected changepoints.
#'
#' @author Kaylea Haynes

#'@examples
#'### Set up parallel environment ###
#'library(doParallel)
#'library(foreach)
#'
#'ncores <- c(10)
#'cl <- makeCluster(ncores)
#'registerDoParallel(cl)
#
#'### Generate some data from the blocks data set ####
#'cpts <- round(c(0,0.1, 0.13, 0.15, 0.23, 0.25, 0.40, 0.44, 0.65, 0.76, 0.78, 0.81,1)*10000)
#'segment_param <- c(0, 4, -1, 2, -2, 3, -1.2, 0.9, 5.2, 2.1, 4.2, 0)
#'data <- NULL
#'
#'for (i in 1:(length(cpts)-1)){
#'   data_new <- rep(segment_param[i],cpts[i+1] - cpts[i])
#'  data <- c(data, data_new)
#'}
#'data <- data + rnorm(10000,0,1)
#'
#'Parallel_PELTSM2(data, 2*log(length(data)), TRUE, sum.stat = norm.sum, cost = norm.mean.cost, minseglen =1, ncores=10)
#' @export

Parallel_PELTSM2 <- function(data, penalty, pruning, sum.stat = norm.sum, cost, minseglen, ncores){
  sumx <- sum.stat(data)
  n <- length(data)
  changepoints <- (foreach(i = 1:(ncores)) %dopar% PELTSM2_step1(sumx, penalty,PRUNE = pruning,cost = cost,coreID = i, ncores = ncores))
  changepointslocs <- sort(unlist(sapply(1:ncores, function(x) changepoints[[x]][[2]])))
  changepointslocs <- unique(changepointslocs)
  ##### STEP 2 #####
  changepointslocs <- changepointslocs[changepointslocs <= n]
  changepointslocs <- changepointslocs[changepointslocs >= 0]
  full_changepoints <- PELTSM2_step2(sumx, penalty,PRUNE = pruning,cost = cost, minseglen = minseglen, cptslocs = c(changepointslocs))[[2]]
  full_changepoints
  return(full_changepoints)
}






