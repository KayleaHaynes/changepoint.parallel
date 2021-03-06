#' @keywords internal
#' @export
PELTSM1_step1<- function(sumx,pen = 2*log(dim(sumx[2]-1)),PRUNE = pruning,cost = cost, minseglen = minseglen, coreID=i, ncores = ncores, lastchangelike = lastchangelike,numchangecpts = numchangecpts, lastchangecpts =lastchangecpts){

  nwhole <- dim(sumx)[2]-1
  if (coreID < ncores){
    sumx <- sumx[,floor((((nwhole/ncores)*(coreID-1)+1):((nwhole/ncores)*(coreID)+1)))]
  }
  else{
    sumx <- sumx[,(floor((nwhole/ncores)*(ncores-1))+1):(nwhole+1)]
  }

  n <- dim(sumx)[2]-1

  if(missing(lastchangelike)) {lastchangelike = array(0,dim = n+1)}
  if(missing(lastchangecpts)) {lastchangecpts = array(0,dim = n+1)}
  if(missing(numchangecpts)) {numchangecpts = array(0,dim = n+1)}

  lastchangelike[1] <- -pen
  numchangecpts[1] <- 0
  lastchangecpts[1] <- 0

  for (j in (minseglen:(2*minseglen-1))){
    lastchangecpts[j+1] = 0
  }

  for(j in (minseglen:(2*minseglen - 1))){
    lastchangelike[j+1] = cost(j,0,sumx)
  }

  for(j in (minseglen:(2*minseglen - 1))){
    numchangecpts[j+1] = 1
  }
  #  if(minseglen>1){
  #  lastchangelike[2:minseglen]=Inf ##cannot have changepoint at 1,minseglen-1.
  #  }

  checklist <- array() #stores the candidate changepoint positions
  checklist[1] <- 0
  checklist[2] <- minseglen



  for (tstar in (2*minseglen):n) {
    if(lastchangelike[tstar+1] == 0){
      tmplike <- lastchangelike[checklist+1]+cost(tstar,checklist,sumx) + pen

      #### Store changepoints and cost function for each tstar ###
      lastchangecpts[tstar+1] <- checklist[min(which(tmplike == min(tmplike[!is.na(tmplike)])))]
      lastchangelike[tstar+1] <- min(tmplike[!is.na(tmplike)])
      numchangecpts[tstar +1] <- numchangecpts[lastchangecpts[tstar+1]+1]+1

      if(PRUNE){
        checklist <- checklist[(tmplike - pen) <= lastchangelike[tstar+1]]
      }
    }
    checklist <- c(checklist,tstar-(minseglen-1))
  }
  cp=n
  while(cp[1]>0){
    cp=c(lastchangecpts[cp[1]],cp) }
  cp <- floor((nwhole/ncores)*(coreID-1)) + cp
  return(list(lastchangecpts, cp, lastchangelike,numchangecpts))
}


################################################################################################
###############################cost functions###################################################
################################################################################################
#' @keywords internal
#' @export
exp.sum=function(data){
  return(cumsum(c(0,data)))
}
#' @keywords internal
#' @export
exp.cost=function(tau,R,sumx){
  return(- (tau-R) * (log(tau-R) -log(sumx[tau+1] - sumx[R+1])) )
}
#' @keywords internal
#' @export
meandat <- function(data){
  x <- mean(data)
  return(x)
}
#' @keywords internal
#' @export
norm.sum = function(data){
  return(matrix(data = c(cumsum(c(0,data)),cumsum(c(0,data^2)), cumsum(c(0,(data-meandat(data))^2))),nrow = 3,byrow = T))
}
#' @keywords internal
#' @export
norm.mean.cost = function(tau,R,sumx){
  return((sumx[2,tau +1] - sumx[2,R+1])-(sumx[1,tau +1] - sumx[1,R+1])^2/(tau - (R)))
}
#' @keywords internal
#' @export
norm.var.cost = function(tau,R,sumx){
  mll.var <- function(x,n){
    neg = x<=0
    x[neg==TRUE] = 10000
    return(n*(log(2*pi) + log(x/n) + 1))
  }
  cost <- mll.var(sumx[3,tau+1] - sumx[3,R+1], tau - R)
  return(cost)

}
#' @keywords internal
#' @export
norm.meanvar.cost = function(tau,R,sumx){

  cost <- array()
  cost[which((tau-R) <= 1)] = 0
  x <- R[which((tau-R) > 1)]
  t2 = sumx[2,tau+1] - sumx[2,x+1]
  t = sumx[1,tau+1] - sumx[1,x+1]
  # cost[which((tau-R) > 1)]= (tau-x) * (log(2*pi) + log((((sumx[2,tau +1] - sumx[2,x+1])-(tau-x) *((sumx[1,tau +1] - sumx[1,x+1])/(tau - (x)))^2))/(tau-x)) + 1)
  sigsq = (t2 - (t^2)/(tau-x))/(tau-x)
  sigsq[which(sigsq <= 0)] = 0.00000000001
  cost[which((tau-R) > 1)] = (tau-x) * (log(2*pi) + log(sigsq) + 1)
  return(cost)
}


