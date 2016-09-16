#' Parallel PELT SM1
#' @description Splits the data into chunks and runs parallel changepoint detection on each chunk then merges the results and runs parallel detection using the detected changepoints in step 1 as the candidate changepoint locations.
#' @param data A vector of data-points within which you with to find changepoints.
#' @param penalty The value of the penalty.
#' @param pruning If true PELT is used, if false Optimal Partitioning is used.
#' @param sum.stat This can be ``norm.sum" or ``exp.sum".
#' @param cost This can be ``norm.mean.cost", ``norm.var.cost", ``norm.meanvar.cost" or ``exp.cost".
#' @param ncores Number of cores to use.
#' @param boundary This can either be ``fixed" or ``adaptive".
#' @param boundary_value If boundary is fixed then this is the number of points to use around the boundary.
#' @param minseglen Minimim length a segment can be
#' @return The detected changepoints.
#'
#' @author Kaylea Haynes

#'@examples
#'
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
#'Parallel_PELTSM1(data, 2*log(length(data)), TRUE, sum.stat = norm.sum, cost = norm.mean.cost, ncores=10, boundary = "fixed",boundary_value = 20, 1)
#' @export

Parallel_PELTSM1 <- function(data, penalty, pruning, sum.stat = norm.sum, cost = norm.mean.cost, ncores=2, boundary = "fixed",boundary_value = 0,
                             minseglen){

  sumx <- sum.stat(data)
  n <- length(data)
  changepoints <- (foreach(i = 1:(ncores)) %dopar% PELTSM1_step1(sumx, pen= penalty,PRUNE = pruning,cost = cost, minseglen = minseglen, coreID = i, ncores = ncores))
  changepointslocs <- (sort(unlist(sapply(1:ncores, function(x) floor(c((n/ncores*(x-1)),changepoints[[x]][[2]],(n/ncores*x)))))))
  changepointslocs <- unique(changepointslocs)

  if (boundary == "fixed"){
    cuts <- c(0,floor((length(data)/ncores)*c(1:ncores)))
    changepointslocs <- sort(unique(c(changepointslocs, unlist(lapply(1:(length(cuts)), function(x) (cuts[x] - boundary_value):(cuts[x]+boundary_value))))))
  }

  else if (boundary == "adaptive"){
    cuts <- floor((length(data)/ncores)*c(1:ncores))

    changepointslocsnew <- unique(unlist(lapply(1:(ncores-1), function(x){c(cuts[x]:min(changepointslocs[changepointslocs >cuts[x] & changepointslocs <= cuts[x+1]]),
                                                                            max(changepointslocs[changepointslocs >= cuts[x] & changepointslocs <cuts[x+1]]):cuts[x+1])})))
    changepointslocsnew <- unique(sort(c(max(changepointslocs[changepointslocs < cuts[1]]):cuts[1],
                                         changepointslocsnew, cuts[ncores-1]:min(changepointslocs[changepointslocs > cuts[ncores-1]]), changepointslocs)))
    changepointslocs <- changepointslocsnew
  }

  ##### STEP 2 #####
  changepointslocs <- changepointslocs[changepointslocs <= n]
  changepointslocs <- changepointslocs[changepointslocs >= 0]
  full_changepoints <- PELTSM1_step2(sumx, penalty,PRUNE = pruning,cost = cost, minseglen = minseglen, cptslocs = changepointslocs)[[2]]
  full_changepoints
  return(full_changepoints)
}


