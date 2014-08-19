



# read sample data

setwd("S:/SMOSSVR/South_Pacific/aggregation/out_emiss")
data = read.csv("SPacific_June_emiss_fil.csv", header = FALSE) 

names(data) = c("ARGO","15H","15V","30H","30V","45H","45V","60H","60V","SST","WS","SSS")



# loop over feature combinations

# create all combinations of features, omitting in-situ ARGO data and retrieved SSS
n <- 10
l <- rep(list(0:1), n)


  # create subsets for training / validation






  #------------------------------------------


comb = expand.grid(l)
names(comb) = c("15H","15V","30H","30V","45H","45V","60H","60V","SST","WS")

# ommit row with no features selected
comb = comb[2:nrow(comb),]

SMOS = as.logical(0)
L3   = as.logical(0)

# loop over feature combinations
for (i in nrow(comb)) {
  
  # training feature selection  
  selstring = as.logical(comb[i,])
  selstring = c(SMOS, selstring, L3)
  x = data[,selstring]
  
  # SMOS TBs / Emissivities
  y = data[,1]
  
  
  # conversion to .pat / sigma-finder  
  # write runner_converter.bat
  # ....  
  
  # write sigma-finder.bat
  # ....
  
  setwd("C:/Projects/6_SVR/svr/Amazon")
  system("runner_converter.bat")
  system("runner_sigma-finder.bat")
  
  # read data from sigma finder output
  
  # n-fold cross validation (currently n=5)  
  for (i in 5) {
   
    
   for (i in 10) {
   
   # svm
   m = svm(x, y, kernel="radial")
   new = predict(m, x)
   
   #plot(y, new)
   #points(x, log(x), col = 2)
   #points(x, new, col = 4)
   }
   
   
   # ....
  
  }
  
   # read best svr parameters for training
   # ....
  
  
   # validation 
   m = svm(x, y, kernel="radial")
   new = predict(m, x)  
   #plot(y, new)
   #points(x, log(x), col = 2)
   #points(x, new, col = 4)
  
   # write results to table
   # ....
}






 


