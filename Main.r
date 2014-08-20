
require(e1071)


# read sample data

setwd("S:/SMOSSVR/South_Pacific/aggregation/out_emiss")

datafile = "SPacific_June_emiss_fil.csv"
data = read.csv(datafile, header = FALSE) 

names(data) = c("ARGO","15H","15V","30H","30V","45H","45V","60H","60V","SST","WS","SSS")








  # create subsets for training / validation






  #------------------------------------------

# create all combinations of features, omitting in-situ ARGO data and retrieved SSS
n <- 10
l <- rep(list(0:1), n)

comb = expand.grid(l)
names(comb) = c("15H","15V","30H","30V","45H","45V","60H","60V","SST","WS")

# ommit row with no features selected
comb = comb[2:nrow(comb),]

# dataframe with svr performences for feature combinations
comb2   = comb
comb2$R = NA

SMOS = as.logical(0)
L3   = as.logical(0)

# loop over feature combinations
for (i in nrow(comb)) {
  
  # training feature selection  
  selstring = as.logical(comb[i,])
  # number of features used for training
  nf = sum(selstring)
  
  selstring = c(SMOS, selstring, L3)
  x = data[,selstring]
  
  # SMOS TBs / Emissivities
  y = data[,1]
  
  
  # write selected feature combination to file for .pat conversion
  # .... requires correct formatting
  
  
  # conversion to .pat / sigma-finder 
  setwd("C:/Projects/6_SVR/svr/Amazon")  
  
  callconverter = paste("converter -r -l -d ", datafile,".txt -o", datafile,".pat -f", nf)
  system(callconverter)
  
  callsigma = paste("sigma-finder -r -d ", datafile, ".pat -o", datafile, "_sigma.txt -s 10")
  system(callsigma)
  
  # read data from sigma finder output  
  sigmas<-rep(rep(NA),12)
  
  n_sigmas=paste(data, "_sigma.txt")
  f_sigmas=file(n_sigmas, open="r")
  linn=readLines(f_sigmas)  
  
  sigmas[1:10] = linn[5:14]
  
  #C and epsilon parameter
  sigmas[11]   = as.numeric(substr(linn[17], nchar(linn[17])-14+1, nchar(linn[17])))
  sigmas[12]   = as.numeric(substr(linn[20], nchar(linn[20])-9+1,  nchar(linn[20])))
  close(f_sigmas)      
  
  # create dataframe to save performence of individual runs
  out = matrix(data = NA, nrow = 10, ncol = 6, byrow = FALSE, dimnames = NULL)
  #out = as.data.frame(out)
  
  # n-fold cross validation (currently n=5)  
  for (j in 1:5) {
   
    
   for (k in 1:10) {
   
    # svm
    m = svm(x, y, kernel="radial", gamma = sigmas[k], cost = sigmas[11], epsilon = sigmas[12])
    new = predict(m, x)
    
    out[k,j] = cor(y, new)
   
    #plot(y, new)
    #points(x, log(x), col = 2)
    #points(x, new, col = 4)
   }
   
   
   # ....
  
  }
  
  # read best svr parameters for training
  for (l in 1:10) {
    out[l,6] = mean(out[l,1:5])  
  }
  
  # location of maximum mean
  bfit = which.max(out[,6])
  
  # use parameters with best performence for validation 
  m = svm(x, y, kernel="radial", gamma = sigmas[bfit], cost = sigmas[11], epsilon = sigmas[12])
  new = predict(m, x) 
  
  # write results to table
  comb2[i,11] = cor(y, new)

}