loo.train = function(d.train,part.window=126,ph,vdrop){

      ###
      ### testing the training evaluation function
      ###
      source("training_test_detector.R")
      source("training_detector.R")
      source("evaluate.R")
  
      individs = unique(d.train[,1]) # vector of id's for individuals in the training set
      nInd.train = length(unique(d.train[,1])) #number of individuals in training set
      nEps=15 #number of epsilon quantiles tested
      nCovs=dim(d.train)[2]-3 # numver of features used ofr predictions
      loo.eval=rep(list(),3*nInd.train)# there are 3 criteria we can use to tune epsilon for anomaly detection
      loo.fittest=rep(list(),3*nInd.train)# there are 3 criteria we can use to tune epsilon for anomaly detection
      nLoo=nInd.train-1 #number in leave one out training
      decide=matrix(NA,nr=nCovs,nc=3)
      decide.indx=list()
      
      for(m in 1:3){
        for(j in 1:nInd.train){
          df=d.train[d.train[,1]!=individs[j],]
          train.cov=array(NA,c(nCovs,nEps,3))
          for(h in 1:nCovs){
            for(k in 1:nEps){# loop over epsilon values
              fit.train = training(d=df[,c(1:3,h+3)],pw=part.window,eps=k/100,vd=vdrop[-j])# fit model
              eval.temp=evaluate(alarm=fit.train$alarm,possible.hits=ph[-j],nInd=nLoo,vitdropday = vdrop[-j])
              train.cov[h,k,1]=eval.temp$out.prec
              train.cov[h,k,2]=eval.temp$out.recall
              train.cov[h,k,3]=eval.temp$out.F1# calculate eval = recall
            }
          }
          compare=apply(train.cov,c(1,3),max,na.rm=TRUE)
          for(h in 1:nCovs){
            for(i in 1:3){
              decide[h,i]=min(which(train.cov[h,,i]==compare[h,i]))
            }
          }

          decide=decide/100
          d.ind=d.train[d.train[,1]==individs[j],]
          fit.test= training_test(d=d.ind,pw=part.window,eps=decide[,m],vd=vdrop[j])# fit model
          eval.test=evaluate(alarm=fit.test$alarm,possible.hits=ph[j],nInd=1,vitdropday = vdrop[j])
          loo.indx=j+(m-1)*nInd.train
          loo.eval[[loo.indx]]=eval.test
          loo.fittest[[loo.indx]]=fit.test
      }    
      }
return(list(loo.eval = loo.eval,decide=decide,train.cov=train.cov,compare = compare,loo.fittest=loo.fittest))
}