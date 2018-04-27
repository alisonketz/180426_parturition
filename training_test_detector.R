training_test = function(d,eps,pw=126,vd){
  
  #d is dataframe of features, 
  #the first column is the id
  #the second column is the julian day of observation
  #the third column is the labeled anomaly (1) for vit drop event
  
  nCovs = dim(d)[2]-3
  id = unique(d[,1])
  nInd = length(id)
  ci = 4:dim(d)[2]
  
  if(length(eps)!=nCovs){cat("length(epsilon) != number features, try again \n");return}
  if(length(vd)!=nInd){cat("length(vitdropday) != number individuals, try again \n");return}
  
  #For individual case, run anomaly dection on single individual
  if(nInd==1){
    d.temp=d[d[,2]>=pw,]
    n.temp=dim(d.temp)[1]
    
    if(nCovs==1){#nInd = 1, ncovs=1
      ind.mean = mean(d[d[,2]<pw,ci],na.rm=TRUE)
      ind.sd =  sd(d[d[,2]<pw,ci],na.rm=TRUE)
      
      detect.quant = quantile(d[d[,2]<(vd-2),ci],probs=eps,na.rm=TRUE)
      # detect.quant = quantile(d[,ci],probs=eps,na.rm=TRUE)
      # detect.quant = quantile(d[d[,2]<pw,ci],probs=eps,na.rm=TRUE)
      # detect.quant = quantile(d.temp[,ci],probs=eps,na.rm=TRUE)
      
      threshold.density = dnorm(detect.quant,ind.mean,ind.sd)
      threshold.p = pnorm(detect.quant,mean=ind.mean,sd=ind.sd)
      lower.prob=pnorm(d.temp[,ci],mean=ind.mean,sd=ind.sd)
      
      detect.density = dnorm(d.temp[,ci],ind.mean,ind.sd)
      
      hits.indx = which(is.finite(detect.density) & detect.density<=threshold.density)
      # hits.indx = which(is.finite(detect.density) & detect.density<=threshold.density & lower.prob <= threshold.p)
      outs.indx = which(!(1:n.temp %in% hits.indx))
      alarm = list(rbind(d.temp[hits.indx,]))
      
      results.prob=rep(NA,n.temp)
      results.prob[hits.indx]=1-(lower.prob[hits.indx]/threshold.p)
      results.prob[outs.indx] =1-(1-lower.prob[outs.indx])/(1-threshold.p)
      
    }#end nInd=1,nCovs=1
    else {# if nInd=1 nCovs >1
      
      detect.quant=rep(NA,nCovs)
      for(i in 1:nCovs){
        detect.quant[i] = quantile(d[d[,2]<(vd-2),i+3],probs=eps[i],na.rm=TRUE)
        # detect.quant[i] = quantile(d[,i+3],probs=eps[i],na.rm=TRUE)
        # detect.quant[i] = quantile(d[d[,2]<pw,i+3],probs=eps[i],na.rm=TRUE)
        # detect.quant[i] = quantile(d.temp[,i+3],probs=eps[i],na.rm=TRUE)
      }
      ind.mean = apply(d[d[,2]<pw,ci],2,mean,na.rm=TRUE)
      ind.sd = apply(d[d[,2]<pw,ci],2,sd,na.rm=TRUE)
      Sigma=diag(ind.sd)
      threshold.density = dmvnorm(detect.quant,ind.mean,Sigma)
      threshold.p = pmvnorm(lower=rep(-Inf,nCovs),upper=detect.quant,mean=ind.mean,sigma=Sigma)[1]
      
      detect.density = rep(NA,n.temp)
      lower.prob=rep(NA,n.temp)
      for(i in 1:n.temp){
        up=as.numeric(d.temp[i,ci])
        up.na=is.na(up)
        detect.density[i] = dmvnorm(up[!up.na],ind.mean[!up.na],diag(ind.sd[!up.na]))
        lower.prob[i]=pmvnorm(lower=rep(-Inf,nCovs-sum(up.na)),upper = up[!up.na],mean=ind.mean[!up.na],sigma=diag(ind.sd[!up.na]))[1]
      }
      
      # hits.indx = which(is.finite(detect.density) & detect.density<=threshold.density & lower.prob <= threshold.p)
      hits.indx = which(is.finite(detect.density) & detect.density<=threshold.density)
      
      outs.indx = which(!(1:n.temp %in% hits.indx))
      alarm=list(rbind(d.temp[hits.indx,]))
      
      results.prob=rep(NA,n.temp)
      results.prob[hits.indx]=1-(lower.prob[hits.indx]/threshold.p)
      results.prob[outs.indx] =1-(1-lower.prob[outs.indx])/(1-threshold.p)
    }#end else nCovs>1,nInd = 1
  }#end nInd =1
  else{ # nInd>1       
    n.temp=rep(NA,nInd)
    for(j in 1:nInd){
      d.temp1=d[d[,1]==id[j],]
      d.temp=d.temp1[d.temp1[,2]>=pw,]
      n.temp[j]=dim(d.temp)[1]
    }
    n.temp.max=max(n.temp)
    
    alarm=rep(list(),nInd)
    hits.indx=rep(list(),nInd)
    outs.indx=rep(list(),nInd)
    threshold.p=rep(NA,nInd)
    threshold.density=rep(NA,nInd)
    detect.density=matrix(NA,nr=n.temp.max,nc=nInd)
    lower.prob=matrix(NA,nr=n.temp.max,nc=nInd)
    results.prob=matrix(NA,nr=n.temp.max,nc=nInd)
    
    if(nCovs==1){ #nInd>1, nCovs = 1
      detect.quant=rep(NA,nInd)
      ind.mean=rep(NA,nInd)
      ind.sd=rep(NA,nInd)
      for(j in 1:nInd){
        d.temp1=d[d[,1]==id[j],]
        
        ind.mean[j] = mean(d.temp1[d.temp1[,2]<pw,ci],na.rm=TRUE)
        ind.sd[j] = sd(d.temp1[d.temp1[,2]<pw,ci],na.rm=TRUE)
        
        detect.quant[j]=quantile(as.numeric(d.temp1[d.temp1[,2]<(vd[j]-2),ci]),eps,na.rm=TRUE)
        # detect.quant[j]=quantile(as.numeric(d.temp1[,ci]),eps,na.rm=TRUE)
        # detect.quant[j]=quantile(as.numeric(d.temp1[d.temp1[,2]<pw,ci]),eps,na.rm=TRUE)
        # detect.quant[j]=quantile(as.numeric(d.temp1[d.temp1[,2]>=pw,ci]),eps,na.rm=TRUE)
        
        threshold.density[j] = dnorm(detect.quant[j],ind.mean[j],ind.sd[j])
        threshold.p[j] = pnorm(detect.quant[j],mean=ind.mean[j],sd=ind.sd[j])
        
        d.temp=d.temp1[d.temp1[,2]>=pw,]
        
        #probability of hit/nothit
        for(i in 1:n.temp[j]){
          lower.prob[i,j]=pnorm(as.numeric(d.temp[i,ci]),mean=ind.mean[j],ind.sd[j])
          detect.density[i,j] = dnorm(as.numeric(d.temp[i,ci]),ind.mean[j],ind.sd[j])
        }
        hits.indx[[j]] = which(is.finite(detect.density[,j]) & detect.density[,j]<=threshold.density[j])
        # hits.indx[[j]] = which(is.finite(detect.density[,j]) & detect.density[,j]<=threshold.density[j] & lower.prob[,j] <= threshold.p[j])
        outs.indx[[j]] = which(!(1:n.temp[j] %in% hits.indx[[j]]))
        
        alarm[[j]]=d.temp[hits.indx[[j]],]
        
        results.prob[hits.indx[[j]],j] = 1-(lower.prob[hits.indx[[j]],j]/threshold.p[j])
        results.prob[outs.indx[[j]],j] =1-(1-lower.prob[outs.indx[[j]],j])/(1-threshold.p[j])
        
      }#endfor nInd
    }#endif nCovs==1
    else{#nind>1,nCovs >1
      
      detect.quant=matrix(NA,nr=nCovs,nc=nInd)
      ind.mean=matrix(NA,nr=nCovs,nc=nInd)
      ind.sd=matrix(NA,nr=nCovs,nc=nInd)
      
      for(j in 1:nInd){
        
        d.temp1=d[d[,1]==id[j],]
        
        ind.mean[,j] = apply(d.temp1[d.temp1[,2]<pw,ci],2,mean,na.rm=TRUE)
        ind.sd[,j] = apply(d.temp1[d.temp1[,2]<pw,ci],2,sd,na.rm=TRUE)
        
        for(i in 1:nCovs){
          detect.quant[i,j] = quantile(d.temp1[d.temp1[,2]<(vd[j]-2),i+3],probs=eps[i],na.rm=TRUE)
          # detect.quant[i,j] = quantile(d.temp1[,i+3],probs=eps[i],na.rm=TRUE)
          # detect.quant[i,j] = quantile(d.temp1[d.temp1[,2]<pw,i+3],probs=eps[i],na.rm=TRUE)
          # detect.quant[i,j] = quantile(d.temp1[d.temp1[,2]>=pw,i+3],probs=eps[i],na.rm=TRUE)
        }
        
        threshold.density[j] = dmvnorm(detect.quant[,j],ind.mean[,j],diag(ind.sd[,j]))
        threshold.p[j] = pmvnorm(lower=rep(-Inf,nCovs),upper=detect.quant[,j],mean=ind.mean[,j],sigma=diag(ind.sd[,j]))[1]
        
        d.temp=d.temp1[d.temp1[,2]>=pw,]
        
        for(i in 1:n.temp[j]){
          up=as.numeric(d.temp[i,ci])
          up.na=is.na(up)
          detect.density[i,j] = dmvnorm(up[!up.na],mean=ind.mean[!up.na,j],sigma=diag(ind.sd[!up.na,j]))
          lower.prob[i,j]=pmvnorm(lower=rep(-Inf,nCovs-sum(up.na)),upper = up[!up.na],mean=ind.mean[!up.na,j],sigma=diag(ind.sd[!up.na,j]))[1]
        }
        hits.indx[[j]] = which(is.finite(detect.density[,j]) & detect.density[,j]<=threshold.density[j])
        # hits.indx[[j]] = which(is.finite(detect.density[,j]) & detect.density[,j]<=threshold.density[j] & lower.prob[,j] <= threshold.p[j])
        outs.indx[[j]] = which(!(1:n.temp[j] %in% hits.indx[[j]]))
        
        alarm[[j]]=d.temp[hits.indx[[j]],]
        
        results.prob[hits.indx[[j]],j] = 1-(lower.prob[hits.indx[[j]],j]/threshold.p[j])
        results.prob[outs.indx[[j]],j] =1-(1-lower.prob[outs.indx[[j]],j])/(1-threshold.p[j])
        
      }#endfor
    }#end else
  }#end group nInd>1 else
  
  return(list(alarm=alarm,
              detect.density=detect.density,
              threshold.density=threshold.density,
              detect.quant=detect.quant,
              results.prob=results.prob,
              hits.indx=hits.indx,
              outs.indx=outs.indx,
              lower.prob=lower.prob,
              ind.mean=ind.mean,
              ind.sd=ind.sd
  ))
  
}
