
require(e1071)


# read sample data

setwd("S:/SMOSSVR/South_Pacific/aggregation/out_emiss")

file = "SPacific_June_emiss_fil.csv"
data = read.csv(file, header = FALSE) 

names(data) = c("ARGO","15H","15V","30H","30V","45H","45V","60H","60V","SST","WS","SSS")


  # create subsets for training / validation  
  lim = nrow(data)/10
  lim = floor(lim) 

  training = data[1:(lim*5),]
  testing  = data[((lim*5)+1):nrow(data),]

  t1 = training[1:lim,]
  t2 = training[(lim+1):(lim*2),]
  t3 = training[((lim*2)+1):(lim*3),]
  t4 = training[((lim*3)+1):(lim*4),]
  t5 = training[((lim*4)+1):(lim*5),]

  SMOS_training = data[1:(lim*5),1]
  SMOS_testing  = data[((lim*5)+1):nrow(data),1]

  SMOS_t1 = training[1:lim,1]
  SMOS_t2 = training[(lim+1):(lim*2),1]
  SMOS_t3 = training[((lim*2)+1):(lim*3),1]
  SMOS_t4 = training[((lim*3)+1):(lim*4),1]
  SMOS_t5 = training[((lim*4)+1):(lim*5),1]


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
for (i in 1:nrow(comb)) {
  
  # training feature selection  
  selstring = as.logical(comb[i,])
  # number of features used for training
  nf = sum(selstring)
  
  selstring = c(SMOS, selstring, L3)
  tr_sub = training[,selstring]
  tr_sub = testing[,selstring]
  
  tr1_sub = t1[,selstring]
  tr2_sub = t2[,selstring]
  tr3_sub = t3[,selstring]
  tr4_sub = t4[,selstring]
  tr5_sub = t5[,selstring]
  
  test_sub = testing[,selstring]
  
  
  # write selected feature combination to file for .pat conversion
  # .... requires correct formatting
  
  # SMOS TBs / Emissivities
  yy = data[,1]
  xx = data[,selstring]
  xx = as.data.frame(xx)
  
  for (m in 1:nf) {
   for(n in 1:nrow(xx)) {  
    xx[n,m] = paste(m,":",xx[n,m], sep="")
   }
  }
  
  xx = cbind(yy,xx)
  
  # ... format does not seem to work ...
  
  datafile  = "sigma_output"
  datafile2 = "sigma_output.txt"
  write.table(format(xx, digits=5), datafile2, quote=FALSE, row.names=FALSE, col.names=FALSE, sep = "\t",)      
  
  
  
  # conversion to .pat / sigma-finder 
  #setwd("C:/Projects/6_SVR/svr/Amazon")  
  
  callconverter = paste("converter -r -l -d ", datafile, ".txt -o ", datafile,".pat -f ", nf, sep="")
  system(callconverter)
  
  callsigma = paste("sigma-finder -r -d ", datafile, ".pat -o ", datafile, "_sigma.txt -s 10", sep="")
  system(callsigma)
  
  # read data from sigma finder output  
  sigmas<-rep(rep(NA),12)
  
  n_sigmas=paste(datafile, "_sigma.txt", sep="")
  f_sigmas=file(n_sigmas, open="r")
  linn=readLines(f_sigmas)  
  
  sigmas[1:10] = as.numeric(linn[5:14])
  
  
  #C and epsilon parameter
  sigmas[11]   = as.numeric(substr(linn[17], nchar(linn[17])-14+1, nchar(linn[17])))
  sigmas[12]   = as.numeric(substr(linn[20], nchar(linn[20])-9+1,  nchar(linn[20])))
  sigmas = as.numeric(sigmas)
  sigmas[1:10] = 1 / sigmas[1:10]
  
  close(f_sigmas)      
  
  # create dataframe to save performence of individual runs
  out = matrix(data = NA, nrow = 10, ncol = 6, byrow = FALSE, dimnames = NULL)
  #out = as.data.frame(out)
  
  # n-fold cross validation (currently n=5)  
  for (j in 1:5) {
    
    x = paste("tr",j,"_sub",sep="")
    y = paste("SMOS_t",j,sep="")
   
    
   for (k in 1:10) {
   
    # svm
    m = svm(get(x), get(y), kernel="radial", gamma = sigmas[k], cost = sigmas[11], epsilon = sigmas[12])
    new = predict(m, get(x))
    
    out[k,j] = cor(get(y), new)
   
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
  m = svm(test_sub, SMOS_testing, kernel="radial", gamma = sigmas[bfit], cost = sigmas[11], epsilon = sigmas[12])
  new = predict(m, test_sub)
  
  # write results to table
  comb2[i,11] = cor(SMOS_testing, new)
}

comb2[,12] = comb2[,11] * comb2[,11]