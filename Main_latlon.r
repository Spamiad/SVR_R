
require(e1071)


# read sample data

setwd("C:/ESA_GTP/svr6months")

file = "out_latlon_emiss_format.csv"
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


  train1 = rbind(t1,t2,t3,t4)
  test1  = t5

  train2 = rbind(t1,t2,t3,t5)
  test2  = t4

  train3 = rbind(t1,t2,t4,t5)
  test3  = t3

  train4 = rbind(t1,t3,t4,t5)
  test4  = t2
 
  train5 = rbind(t2,t3,t4,t5)
  test5  = t1

  ARGO_train1 = rbind(t1,t2,t3,t4)
  ARGO_train1 = ARGO_train1[,1]

  ARGO_test1  = t5
  ARGO_test1 = ARGO_test1[,1]

  ARGO_train2 = rbind(t1,t2,t3,t5)
  ARGO_train2 = ARGO_train2[,1]

  ARGO_test2  = t4
  ARGO_test2 = ARGO_test2[,1]

  ARGO_train3 = rbind(t1,t2,t4,t5)
  ARGO_train3 = ARGO_train3[,1]

  ARGO_test3  = t3
  ARGO_test3 = ARGO_test3[,1]

  ARGO_train4 = rbind(t1,t3,t4,t5)
  ARGO_train4 = ARGO_train4[,1]

  ARGO_test4  = t2
  ARGO_test4 = ARGO_test4[,1]

  ARGO_train5 = rbind(t2,t3,t4,t5)
  ARGO_train5 = ARGO_train5[,1]

  ARGO_test5  = t1
  ARGO_test5  = ARGO_test5[,1]

  #------------------------------------------

# create all combinations of features, omitting in-situ ARGO data and retrieved SSS
n <- 10
l <- rep(list(0:1), n)

comb = expand.grid(l)
names(comb) = c("15H","15V","30H","30V","45H","45V","60H","60V","SST","WS")

# ommit row with no features selected
comb = comb[2:nrow(comb),]


    # manual feature selection to save processing time
r1 = c(1,1,0,0,0,0,0,0,0,0,0,0)
r2 = c(1,1,1,1,1,1,1,1,1,1,1,1)
r3 = c(1,1,0,0,0,0,0,0,0,0,1,0)


comb = rbind(r1,r2,r3)
comb = as.data.frame(comb)
names(comb) = c("lat","lon","15H","15V","30H","30V","45H","45V","60H","60V","SST","WS")

#comb = comb[2:nrow(comb),]


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
  
  tr1_sub = train1[,selstring]
  tr2_sub = train2[,selstring]
  tr3_sub = train3[,selstring]
  tr4_sub = train4[,selstring]
  tr5_sub = train5[,selstring]
  
  ts1_sub = test1[,selstring]
  ts2_sub = test2[,selstring]
  ts3_sub = test3[,selstring]
  ts4_sub = test4[,selstring]
  ts5_sub = test5[,selstring]
  
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
  sigmas[12]   = as.numeric(substr(linn[20], nchar(linn[20])-7+1,  nchar(linn[20])))
  sigmas = as.numeric(sigmas)
  sigmas[1:10] = 1 / sigmas[1:10]
  
  close(f_sigmas)      
  
  # create dataframe to save performence of individual runs
  out = matrix(data = NA, nrow = 10, ncol = 6, byrow = FALSE, dimnames = NULL)
  #out = as.data.frame(out)
  
  # n-fold cross validation (currently n=5)  
  for (j in 1:5) {
    
    x = paste("tr",j,"_sub",sep="")
    y = paste("ARGO_train",j,sep="")
   
    z = paste("ts",j,"_sub",sep="") 
    zz= paste("ARGO_test",j,sep="")
      
   for (k in 1:10) {
   
    # svm
    m = svm(get(x), get(y), kernel="radial", gamma = sigmas[k], cost = sigmas[11], epsilon = 0.1, cachesize = 250, tolerance = 0.000010)
    new = predict(m, get(z))
    
    out[k,j] = cor(get(zz), new)
   
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
  m = svm(tr_sub, SMOS_training, kernel="radial", gamma = sigmas[bfit], cost = sigmas[11], epsilon = 0.1, cachesize = 250, tolerance = 0.000010)
  new = predict(m, test_sub)
  
  # write results to table
  comb2[i,11] = cor(SMOS_testing, new)
}

comb2[,12] = comb2[,11] * comb2[,11]
