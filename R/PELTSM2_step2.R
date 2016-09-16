#' @keywords internal
#' @export
PELTSM2_step2<- function(sumx, pen = 2*log(dim(sumx)[2]-1),PRUNE = pruning,cost = norm.mean.cost, minseglen = 1, cptslocs = changepoints){

  n <- dim(sumx)[2]-1
  lastchangelike = array(0,dim = n+1)
  lastchangecpts = array(0,dim = n+1)
  numchangecpts = array(0,dim = n+1)
  lastchangelike[1] <- -pen

  checklist <- array() #stores the candidate changepoint positions
  checklist[1] <- 0
  for (tstar in cptslocs[-1]){
    checklist_remove <- which((tstar - checklist) < minseglen)
    if (length(checklist_remove > 0)){
      checklist1 <- checklist[-checklist_remove]
    }
    else {checklist1 <- checklist}
    tmplike <- lastchangelike[checklist1+1]+cost(tstar,checklist1,sumx) + pen
    #### Store changepoints and cost function for each tstar ###
    lastchangecpts[tstar + 1] <- checklist[min(which(tmplike == min(tmplike[!is.na(tmplike)])))]
    lastchangelike[tstar+1] <- min(tmplike[!is.na(tmplike)])

    if(PRUNE){
      checklist1 <- checklist1[(tmplike - pen) <= lastchangelike[tstar+1]]
    }
    if (length(checklist_remove > 0)){
      checklist <- c(checklist1,checklist[checklist_remove],tstar)
    }
    else{checklist <- c(checklist1,tstar)}
  }
  cp=tstar
  lastchangecpts2 <- lastchangecpts[-1]
  while(cp[1]>0){
    cp=c(lastchangecpts2[cp[1]],cp) }
  return(list(lastchangecpts,cp , lastchangelike))
}
