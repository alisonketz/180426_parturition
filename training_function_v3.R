loo.train.v3 = function(d.train,part.window=126,ph,vdrop){
  
  ###
  ### testing the training evaluation function
  ###
  source("training_detector.R")
  source("training_test_detector.R")
  source("evaluate.R")
  
  individs = unique(d.train[,1]) # vector of id's for individuals in the training set
  nInd.train = length(unique(d.train[,1])) #number of individuals in training set
  nEps=15 #number of epsilon quantiles tested
  nCovs=dim(d.train)[2]-3 # numver of features used ofr predictions
  loo.eval=rep(list(),3*nInd.train)# there are 3 criteria we can use to tune epsilon for anomaly detection
  loo.fittest=rep(list(),3*nInd.train)# there are 3 criteria we can use to tune epsilon for anomaly detection
  nLoo=1 #number in leave one out training
  decide=matrix(NA,nr=nCovs,nc=3)
  decide.indx=list()

  eval.train=array(NA,c(nInd.train,nCovs,nEps,3))
    for(j in 1:nInd.train){
      d.ind=d.train[d.train[,1]==individs[j],]
      for(h in 1:nCovs){
        for(k in 1:nEps){# loop over epsilon values
          fit.train = training(d=d.ind[,c(1:3,h+3)],pw=part.window,eps=k/100,vd=vdrop[j])# fit model
          eval.temp=evaluate(alarm=fit.train$alarm,possible.hits=ph[j],nInd=nLoo,vitdropday = vdrop[j])
          eval.train[j,h,k,1]=eval.temp$out.prec
          eval.train[j,h,k,2]=eval.temp$out.recall
          eval.train[j,h,k,3]=eval.temp$out.F1# calculate eval = recall
        }#end k
      }#end h
    }#end j
  
  browser()
  for(m in 1:3){
    eval.train[,,,m]
  }
  
  eval.train[,,,2]
  for(j in 1:nInd.train){
    df=d.train[d.train[,1]!=individs[j],]
    fit.train = training(d=df,pw=part.window,eps=.1,vd=vdrop[j])
    fit.test = evaluate(alarm=fit.train$alarm,possible.hits=ph[-j],nInd=nLoo,vitdropday = vdrop[-j])
  }




apply(eval.train[,,,2],c(1,2),mean)


  return(list(loo.eval = loo.eval,decide=decide,eval.train=eval.train,compare = compare,loo.fittest=loo.fittest))
}