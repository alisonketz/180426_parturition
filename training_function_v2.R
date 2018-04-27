loo.train.v2 = function(d.train,part.window=126,ph,vdrop){

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
      nLoo=nInd.train #number in leave one out training
      decide=matrix(NA,nr=nCovs,nc=3)
      eval=array(NA,c(nInd.train,3,3))
      
      for(m in 1:3){
          train.cov=array(NA,c(nCovs,nEps,3))
          for(h in 1:nCovs){
            for(k in 1:nEps){# loop over epsilon values
              fit.train = training(d=d.train[,c(1:3,h+3)],pw=part.window,eps=k/100,vd=vdrop)# fit model
              eval.temp=evaluate(alarm=fit.train$alarm,possible.hits=ph,nInd=nLoo,vitdropday = vdrop)
              train.cov[h,k,1]=eval.temp$out.prec
              train.cov[h,k,2]=eval.temp$out.recall
              train.cov[h,k,3]=eval.temp$out.F1# calculate eval = recall
            }
          }
          decide=apply(train.cov,c(1,3),which.max)/100

          for(j in 1:nInd.train){
              d.ind=d.train[d.train[,1]==individs[j],]
              fit.test=training_test(d=d.ind,pw=part.window,eps=decide[,m],vd=vdrop[j])# fit model
              eval.test=evaluate(alarm=fit.test$alarm,possible.hits=ph[j],nInd=1,vitdropday = vdrop[j])
              eval[j,m,1]=eval.test$out.prec
              eval[j,m,2]=eval.test$out.recall
              eval[j,m,3]=eval.test$out.F1
          }#end j
      }#end m  

return(list(decide=decide,train.cov=train.cov,eval=eval))
}