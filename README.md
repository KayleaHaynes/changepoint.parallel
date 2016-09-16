# changepoint.parallel
This package is a very basic package to run changepoint detection.  The code here runs algorithms discussed in my PhD thesis.  Please email me for more information.
There are 2 different methods used to split the data over multiple cores.  These are SM1 and SM2 which are depicted below. 

##### Example 
```
### Set up parallel environment ###
library(doParallel)
library(foreach)

ncores <- c(10)
cl <- makeCluster(ncores)
registerDoParallel(cl)

### Generate some data from the blocks data set ####
cpts <- round(c(0,0.1, 0.13, 0.15, 0.23, 0.25, 0.40, 0.44, 0.65, 0.76, 0.78, 0.81,1)*10000)
segment_param <- c(0, 4, -1, 2, -2, 3, -1.2, 0.9, 5.2, 2.1, 4.2, 0)
data <- NULL

for (i in 1:(length(cpts)-1)){
   data_new <- rep(segment_param[i],cpts[i+1] - cpts[i])
  data <- c(data, data_new)
}
data <- data + rnorm(10000,0,1)

Parallel_PELTSM1(data, 2*log(length(data)), TRUE, sum.stat = norm.sum, cost = norm.mean.cost, ncores=10, boundary = "fixed",boundary_value = 20, 1)
Parallel_PELTSM2(data, 2*log(length(data)), TRUE, sum.stat = norm.sum, cost = norm.mean.cost, minseglen =1, ncores=10)
```
