loo.train.v1 = function(d.train,part.window=126,ph,vdrop){

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
      decide=rep(list(matrix(NA,nr=nCovs,nc=3)),nInd.train)
      eval=array(NA,c(nInd.train,3,3))

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
              train.cov[h,k,3]=eval.temp$out.F1
            }
          }
          decide[[j]]=apply(train.cov,c(1,3),which.max)/100
          d.ind=d.train[d.train[,1]==individs[j],]
          fit.test=training_test(d=d.ind,pw=part.window,eps=decide[[j]][,m],vd=vdrop[j])# fit model
          eval.test=evaluate(alarm=fit.test$alarm,possible.hits=ph[j],nInd=1,vitdropday = vdrop[j])
          eval[j,m,1]=eval.test$out.prec
          eval[j,m,2]=eval.test$out.recall
          eval[j,m,3]=eval.test$out.F1
        }    
      }
      #epsilon averaged over maximum of results based on precision
      
      prec.max = max(eval[,1,1])
      prec.indx = which(eval[,1,1]==prec.max)
      prec.save=rep(NA,nCovs)
      for(i in prec.indx){
        prec.save=cbind(prec.save,decide[[i]][,1])
      }
      prec.save=prec.save[,-1]

      #epsilon averaged over maximum of results based on recall
      
      recall.max = max(eval[,2,2])
      recall.indx = which(eval[,2,2]==recall.max)
      recall.save=rep(NA,nCovs)
      for(i in recall.indx){
        recall.save=cbind(recall.save,decide[[i]][,2])
      }
      recall.save=recall.save[,-1]

      #epsilon averaged over maximum of results based on F1
      F1.max = max(eval[,3,3])
      F1.indx = which(eval[,3,3]==F1.max)
      F1.save=rep(NA,nCovs)
      for(i in F1.indx){
          F1.save=cbind(F1.save,decide[[7]][,3])
      }
      F1.save=F1.save[,-1]

      if(!is.matrix(prec.save)){
        prec.eps = prec.save
      }else{
        prec.eps = apply(prec.save,1,median)
      }
      if(!is.matrix(recall.save)){
        recall.eps = recall.save
      }else{
        recall.eps = apply(recall.save,1,median)
      }
      if(!is.matrix(F1.save)){
        F1.eps = F1.save
      }else{
        F1.eps = apply(F1.save,1,median)
      }
      epsilon.star = cbind(prec.eps,recall.eps,F1.eps)
      
      
      ###
      ### average across all example of decide
      ###
      decide.prec=apply(sapply(decide,function(x){x[,1]}),1,mean)
      decide.recall=apply(sapply(decide,function(x){x[,2]}),1,mean)
      decide.F1=apply(sapply(decide,function(x){x[,3]}),1,mean)
      
      decide.star = cbind(decide.prec,decide.recall,decide.F1)

      return(list(eval=eval,decide=decide,train.cov=train.cov,epsilon.star=epsilon.star,decide.star=decide.star))
}