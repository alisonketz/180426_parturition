loo.train.recall = function(d.train,part.window=126,ph,vdrop){
  
  ###
  ### testing the training evaluation function
  ###
  source("training_detector.R")
  source("training_detector.R")
  source("evaluate.R")
  
  individs = unique(d.train[,1])
  nInd.train = length(unique(d.train[,1]))
  nEps=15 #number of epsilon quantiles tested
  nCovs=dim(d.train)[2]-3
  loo.eval=rep(list(),nInd.train)# there are 3 criteria we can use to tune epsilon for anomaly detection
  loo.fittest=rep(list(),nInd.train)# there are 3 criteria we can use to tune epsilon for anomaly detection
  n.loo=nInd.train-1 #number in leave one out training
  decide=matrix(NA,nr=nCovs,nc=nEps)
  decide.indx=list()
  
  train.cov=array(NA,c(nCovs,nEps,nInd.train))
  for(j in 1:nInd.train){
      df=d.train[d.train[,1]!=individs[j],]# subset data to leave out one individual
      d.ind=d.train[d.train[,1]==individs[j],] # subset data for the left out individual
      for(h in 1:nCovs){#loop over features
        for(k in 1:nEps){# loop over epsilon values
          fit.train = training(d=df[,c(1:3,h+3)],pw=part.window,eps=k/100,vd=vdrop[-j])# fit model leaving one individual out
          eval.temp=evaluate(alarm=fit.train$alarm,possible.hits=ph[-j],nInd=n.loo,vitdropday = vdrop[-j]) # evaluate the fitted model per individual
          train.cov[h,k,j]=eval.temp$out.recall # record results of training
          fit.test=training_test(d=d.ind,pw=part.window,eps=k/100,vd=vdrop[j])# fit model to the left out individual to test
          eval.test=evaluate(alarm=fit.test$alarm,possible.hits=ph[j],nInd=1,vitdropday = vdrop[j])# evaluate the test
          eval.cov[h,k,j]=eval.test$out.recall # record results of test
        }
      }
      # loo.fittest[[j]]=fit.test #save the fit
  }    
  return(list(eval.cov=eval.cov,train.cov=train.cov))
}
