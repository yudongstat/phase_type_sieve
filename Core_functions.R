beta_min <- function(mypara,y=numeric(0),censor=numeric(0),delta=numeric(0),
                   m=numeric(0),m1=numeric(0),mm1=numeric(0),
                   long_Hs=numeric(0),long_Hs_transit=numeric(0),
                   long_Hs1=numeric(0),long_Hs1_transit=numeric(0),
                   long_Hs3=numeric(0),long_Hs3_transit=numeric(0),
                   p0ps = numeric(0), p1ps = numeric(0), p2ps = numeric(0),
                   covariate_reported = numeric(0), covariate_unreported = numeric(0)){
  s1 <- abs(diag(p0ps))[1]
  s2 <- abs(diag(p1ps))[1]
  s3 <- abs(diag(p2ps))[1]
  res<-0
  mybeta<-mypara[1:covariate_dim]
  myalpha<-mypara[(covariate_dim+1):(2*covariate_dim)]
  for(sub1 in 1:length(y)){
    res<-res-sum(diag(long_Hs[sub1,,])[1:m])*exp(-sum(covariate_reported[sub1,]*mybeta))*s1
    res<-res-sum(diag(long_Hs1[sub1,,])[1:m1])*exp(-sum(covariate_reported[sub1,]*myalpha))*s2
  }
  for(sub2 in 1:length(censor)){
    res<-res-exp(-sum(covariate_unreported[sub2,]*mybeta))*sum(diag(long_Hs3[sub2,,])[1:(m*m1)])*s1
    res<-res-exp(-sum(covariate_unreported[sub2,]*myalpha))*sum(diag(long_Hs3[sub2,,])[1:(m*m1)])*s2
  }
  for(sub1 in 1:length(y)){
    for(jj in 1:(m-1)){
      res<-res-sum(covariate_reported[sub1,]*mybeta)*long_Hs_transit[sub1,jj,jj+1]
    }
    res<-res-sum(covariate_reported[sub1,]*mybeta)*delta[sub1]
    for(jj in 1:(m1-1)){
      res<-res-sum(covariate_reported[sub1,]*myalpha)*long_Hs1_transit[sub1,jj,jj+1]
    }
    res<-res-sum(covariate_reported[sub1,]*myalpha)*(1-delta[sub1])
  }
  for(sub2 in 1:length(censor)){
    for(qqq in 1:(m-1)){
      res<-res-sum(covariate_unreported[sub2,]*mybeta)*sum(diag(long_Hs3_transit[sub2,((qqq-1)*m1+1):(qqq*m1),(qqq*m1+1):((qqq+1)*m1)]))
      res<-res-sum(covariate_unreported[sub2,]*myalpha)*sum(diag(long_Hs3_transit[sub2,((qqq-1)*m1+1):(qqq*m1-1),((qqq-1)*m1+2):(qqq*m1)]))
    }
    qqq<-m
    res<-res-sum(covariate_unreported[sub2,]*myalpha)*sum(diag(long_Hs3_transit[sub2,((qqq-1)*m1+1):(qqq*m1-1),((qqq-1)*m1+2):(qqq*m1)]))
    res<-res-sum(covariate_unreported[sub2,]*mybeta)*sum(long_Hs3_transit[sub2,((m-1)*m1+1):(m*m1-1),(m*m1+1):mm1])
    res<-res-sum(covariate_unreported[sub2,]*myalpha)*sum(long_Hs3_transit[sub2,(1:(m-1))*m1,(m*m1+1):mm1])
    res<-res+log(s1*exp(-sum(covariate_unreported[sub2,]*mybeta))+s2*exp(-sum(covariate_unreported[sub2,]*myalpha)))*sum(long_Hs3_transit[sub2,m*m1,(m*m1+1):mm1])
  }
  return(-res)
}

emcore <- function(y,delta,y1,censor,pips,pi1ps,pi2ps,est_beta,est_alpha,Rmax,p0ps,p1ps,p2ps,nups,nu1ps,nu2ps,covariate_reported,covariate_unreported){
  mu2ps=max(abs(diag(p2ps)))+0.1
  old_y<-y
  old_y1<-y1
  m=length(pips)
  Hs<-matrix(0,m,m)
  newHs<-matrix(0,m,m)
  long_Hs<-array(0,dim = c(length(y),m,m))
  long_Hs_transit<-array(0,dim = c(length(y),m,m))
  Bs<-rep(0,m)
  for(i in 1:length(y)){
    mups<-max(abs(diag(exp(-sum(covariate_reported[i,]*est_beta))*p0ps)))+0.1
    re_D0=exp(-sum(covariate_reported[i,]*est_beta))*p0ps/mups+diag(m)
    if(delta[i]==1){
      vec1<-nups*exp(-sum(covariate_reported[i,]*est_beta))
      hc<-pips*expAtv(A=exp(-sum(covariate_reported[i,]*est_beta))*p0ps,v=vec1,t=y[i])$eAtv
      sumhc<-sum(hc)
      Bs<-Bs+hc/sumhc
    }else{
      vec1<-rep(1,m)
      hc<-pips*expAtv(A=exp(-sum(covariate_reported[i,]*est_beta))*p0ps,v=vec1,t=y[i])$eAtv
      sumhc<-sum(hc)
      Bs<-Bs+hc/sumhc
    }
    CHks=matrix(0,m,m)
    calps=matrix(vec1,m,Rmax+1)
    for(j in 2:(Rmax+1)) calps[,j]=re_D0%*%calps[,j-1]
    cbets=matrix(dpois(Rmax+1,mups*y[i])*pips,Rmax+1,m,byrow=T)
    for(j in Rmax:1){
      cbets[j,]=cbets[j+1,]%*%re_D0+dpois(j,mups*y[i])*pips
    }
    for(j in 1:(Rmax+1)){
      CHks=CHks+calps[,j]%*%t(cbets[j,])
    }
    CHks<-t(CHks)/sumhc/mups
    long_Hs[i,,]<-CHks
    Hs<-Hs+CHks
    newHs<-newHs+CHks*exp(-sum(covariate_reported[i,]*est_beta))*p0ps
    long_Hs_transit[i,,]<-CHks*exp(-sum(covariate_reported[i,]*est_beta))*p0ps
  }
  
  ###consider the drop-out time
  m1=length(pi1ps)
  Hs1<-matrix(0,m1,m1)
  newHs1<-matrix(0,m1,m1)
  long_Hs1<-array(0,dim = c(length(y),m1,m1))
  long_Hs1_transit<-array(0,dim = c(length(y),m1,m1))
  Bs1<-rep(0,m1)
  for(i in 1:length(y)){
    mu1ps<-max(abs(diag(exp(-sum(covariate_reported[i,]*est_alpha))*p1ps)))+0.1
    re_D0=exp(-sum(covariate_reported[i,]*est_alpha))*p1ps/mu1ps+diag(m1)
    if(delta[i]==0){
      vec1<-nu1ps*exp(-sum(covariate_reported[i,]*est_alpha))
      hc<-pi1ps*expAtv(A=exp(-sum(covariate_reported[i,]*est_alpha))*p1ps,v=vec1,t=y[i])$eAtv
      sumhc<-sum(hc)
      Bs1<-Bs1+hc/sumhc
    }else{
      vec1<-rep(1,m1)
      hc<-pi1ps*expAtv(A=exp(-sum(covariate_reported[i,]*est_alpha))*p1ps,v=vec1,t=y[i])$eAtv
      sumhc<-sum(hc)
      Bs1<-Bs1+hc/sumhc
    }
    CHks1=matrix(0,m1,m1)
    calps=matrix(vec1,m1,Rmax+1)
    for(j in 2:(Rmax+1)) calps[,j]=re_D0%*%calps[,j-1]
    cbets=matrix(dpois(Rmax+1,mu1ps*y[i])*pi1ps,Rmax+1,m1,byrow=T)
    for(j in Rmax:1){
      cbets[j,]=cbets[j+1,]%*%re_D0+dpois(j,mu1ps*y[i])*pi1ps
    }
    for(j in 1:(Rmax+1)){
      CHks1=CHks1+calps[,j]%*%t(cbets[j,])
    }
    CHks1<-t(CHks1)/sumhc/mu1ps
    long_Hs1[i,,]<-CHks1
    Hs1<-Hs1+CHks1
    newHs1<-newHs1+CHks1*exp(-sum(covariate_reported[i,]*est_alpha))*p1ps
    long_Hs1_transit[i,,]<-CHks1*exp(-sum(covariate_reported[i,]*est_alpha))*p1ps
  }
  
  ###afterwards, we consider the observed reporting delays
  y1<-sort(y1)
  dy2<-c(y1[1],diff(y1))
  m2=length(pi2ps)
  ss<-length(dy2)
  re_p02=p2ps/mu2ps+diag(m2)
  fks2=matrix(expAtv(t(p2ps),pi2ps,dy2[1])$eAtv,m2,ss)
  bks2=matrix(expAtv(p2ps,nu2ps,dy2[1])$eAtv,m2,ss)
  for(i in 2:ss){
    fks2[,i]=expAtv(t(p2ps),fks2[,i-1],dy2[i])$eAtv
    bks2[,i]=expAtv(p2ps,bks2[,i-1],dy2[i])$eAtv
  }
  cks2=matrix(pi2ps/sum(pi2ps*bks2[,ss]),ss,m2,byrow=T)
  for(i in (ss-1):1){
    cks2[i,]=expAtv(t(p2ps),cks2[i+1,],dy2[i+1])$eAtv+pi2ps/sum(pi2ps*bks2[,i])
  }
  bks2=cbind(nu2ps,bks2)
  Hs2=matrix(0,m2,m2)
  for(i in 1:ss){
    Hks2=matrix(0,m2,m2)
    alps2=matrix(bks2[,i],m2,Rmax+1)
    for(j in 2:(Rmax+1)) alps2[,j]=re_p02%*%alps2[,j-1]
    bets2=matrix(dpois(Rmax+1,mu2ps*dy2[i])*cks2[i,],Rmax+1,m2,byrow=T)
    for(j in Rmax:1){
      bets2[j,]=bets2[j+1,]%*%re_p02+dpois(j,mu2ps*dy2[i])*cks2[i,]
    }
    for(j in 1:(Rmax+1)){
      Hks2=Hks2+alps2[,j]%*%t(bets2[j,])
    }
    Hs2=Hs2+Hks2/mu2ps
  }
  Bs2=rep(0,m2)
  for(i in 1:ss){
    Bs2=Bs2+pi2ps*bks2[,i+1]/sum(pi2ps*bks2[,i+1])
  }
  Hs2<-t(Hs2)
  
  ###for missing data
  mm1<-m*m1+m2
  Hs3<-matrix(0,mm1,mm1)
  newCHks<-matrix(0,mm1,mm1)
  long_Hs3<-array(0,dim = c(length(censor),mm1,mm1))
  long_Hs3_transit<-array(0,dim = c(length(censor),mm1,mm1))
  Bs3<-rep(0,mm1)
  for(sub in 1:length(censor)){
    new_pips<-c(kronecker(pips,pi1ps),rep(0,m2))
    temp_nups<-rep(0,m*m1)
    for(i in 1:m){
      temp_nups[i*m1]<-temp_nups[i*m1]+exp(-sum(covariate_unreported[sub,]*est_alpha))*nu1ps[m1]
      # temp_nups[i*m1-1]<-temp_nups[i*m1-1]+exp(-sum(covariate_unreported[sub,]*est_alpha))*nu1ps[m1-1]
    }
    for(i in 1:m1){
      # temp_nups[(m-2)*m1+i]<-temp_nups[(m-2)*m1+i]+exp(-sum(covariate_unreported[sub,]*est_beta))*nups[m-1]
      temp_nups[(m-1)*m1+i]<-temp_nups[(m-1)*m1+i]+exp(-sum(covariate_unreported[sub,]*est_beta))*nups[m]
    }
    new_p0ps<-cbind((kronecker(exp(-sum(covariate_unreported[sub,]*est_beta))*p0ps,diag(rep(1,m1)))+kronecker(diag(rep(1,m)),exp(-sum(covariate_unreported[sub,]*est_alpha))*p1ps)),temp_nups%*%t(pi2ps))
    new_p0ps<-rbind(new_p0ps,cbind(matrix(0,m2,m*m1),p2ps))
    new_nups<-apply(new_p0ps,MARGIN = 1,sum)*(-1)
    
    hc=new_pips*expAtv(new_p0ps,rep(1,mm1),censor[sub])$eAtv
    sumhc<-sum(hc)
    Bs3<-Bs3+hc/sumhc
    mmups<-max(abs(diag(new_p0ps)))+0.1
    CHks=matrix(0,mm1,mm1)
    re_D0=new_p0ps/mmups+diag(mm1)
    es=rep(1,mm1)
    
    calps=matrix(es,mm1,Rmax+1)
    for(j in 2:(Rmax+1)) calps[,j]=re_D0%*%calps[,j-1]
    cbets=matrix(dpois(Rmax+1,mmups*censor[sub])*new_pips,Rmax+1,mm1,byrow=T)
    for(j in Rmax:1){
      cbets[j,]=cbets[j+1,]%*%re_D0+dpois(j,mmups*censor[sub])*new_pips
    }
    for(j in 1:(Rmax+1)){
      CHks=CHks+calps[,j]%*%t(cbets[j,])
    }
    CHks=t(CHks)/sumhc/mmups
    long_Hs3[sub,,]<-CHks
    Hs3<-Hs3+CHks
    newCHks<-newCHks+CHks*new_p0ps
    long_Hs3_transit[sub,,]<-CHks*new_p0ps
  }
  CHks<-Hs3
  
  
  
  
  ####################################M Step
  ###initial prob vector of lifetime T
  temp_vec1<-Bs
  for(i in 1:m){
    temp_vec1[i]<-temp_vec1[i]+sum(Bs3[((i-1)*m1+1):(i*m1)])
  }
  output_pips<-temp_vec1/sum(temp_vec1)
  ###initial prob vector of censoring time tau
  temp_vec2<-Bs1
  for(i in 1:m1){
    for(j in 1:m){
      temp_vec2[i]<-temp_vec2[i]+Bs3[(j-1)*m1+i]
    }
  }
  output_pi1ps<-temp_vec2/sum(temp_vec2)
  ###initial prob vector of reporting delay X
  # new_Hs1<-Hs1*p1ps
  new_Hs2<-Hs2*p2ps
  temp_vec3<-Bs2
  for(i in 1:m2){
    temp_vec3[i]<-temp_vec3[i]+sum(newCHks[1:(m*m1),(m*m1+i)])
  }
  output_pi2ps<-temp_vec3/sum(temp_vec3)
  
  current_beta1<-optim(par=c(est_beta,est_alpha),fn=beta_min,y=y,censor=censor,delta=delta,m=m,m1=m1,mm1=mm1,
                       long_Hs=long_Hs,long_Hs_transit=long_Hs_transit,
                       long_Hs1=long_Hs1,long_Hs1_transit=long_Hs1_transit,
                       long_Hs3=long_Hs3,long_Hs3_transit=long_Hs3_transit,
                       p0ps = p0ps, p1ps = p1ps, p2ps = p2ps,
                       covariate_reported = covariate_reported,
                       covariate_unreported = covariate_unreported)$par
  new_beta<-c(current_beta1[1:covariate_dim])
  new_alpha<-c(current_beta1[(1:covariate_dim)+covariate_dim])
  ####results
  list(output_pips,output_pi1ps,output_pi2ps,new_beta,new_alpha)
}


loglike <- function(beta,myalpha,pips,pi1ps,pi2ps,y,delta,y1,censor,p0ps,p1ps,p2ps,nups,nu1ps,nu2ps,covariate_reported,covariate_unreported){
  m=length(pips);m1=length(pi1ps);m2=length(pi2ps)
  mm1<-m*m1+m2
  lk<-sum(dphase(y1,ph=ph(alpha=pi2ps,Q=p2ps,xi=nu2ps),log=T))
  for(sub in 1:length(y)){
    if(delta[sub]==1){
      lk<-lk+dphase(y[sub],ph=ph(alpha = pips,Q=p0ps*exp(-sum(beta*covariate_reported[sub,])),xi=nups*exp(-sum(beta*covariate_reported[sub,]))),log = TRUE)
      lk<-lk+pphase(y[sub],ph=ph(alpha = pi1ps,Q=p1ps*exp(-sum(myalpha*covariate_reported[sub,])),xi=nu1ps*exp(-sum(myalpha*covariate_reported[sub,]))),lower.tail = FALSE,log.p=TRUE)
    }else{
      lk<-lk+pphase(y[sub],ph=ph(alpha = pips,Q=p0ps*exp(-sum(beta*covariate_reported[sub,])),xi=nups*exp(-sum(beta*covariate_reported[sub,]))),lower.tail = FALSE,log.p=TRUE)
      lk<-lk+dphase(y[sub],ph=ph(alpha = pi1ps,Q=p1ps*exp(-sum(myalpha*covariate_reported[sub,])),xi=nu1ps*exp(-sum(myalpha*covariate_reported[sub,]))),log = TRUE)
    }
  }
  for(sub in 1:length(censor)){
    new_pips<-c(kronecker(pips,pi1ps),rep(0,m2))
    temp_nups<-rep(0,m*m1)
    for(i in 1:m){
      temp_nups[i*m1]<-temp_nups[i*m1]+exp(-sum(covariate_unreported[sub,]*myalpha))*nu1ps[m1]
    }
    for(i in 1:m1){
      temp_nups[(m-1)*m1+i]<-temp_nups[(m-1)*m1+i]+exp(-sum(covariate_unreported[sub,]*beta))*nups[m]
    }
    new_p0ps<-cbind((kronecker(exp(-sum(covariate_unreported[sub,]*beta))*p0ps,diag(rep(1,m1)))+kronecker(diag(rep(1,m)),exp(-sum(covariate_unreported[sub,]*myalpha))*p1ps)),temp_nups%*%t(pi2ps))
    new_p0ps<-rbind(new_p0ps,cbind(matrix(0,m2,m*m1),p2ps))
    new_nups<-apply(new_p0ps,MARGIN = 1,sum)*(-1)
    lk<-lk+pphase(censor[sub],ph=ph(alpha=new_pips,Q=new_p0ps,xi=new_nups),lower.tail = FALSE,log.p=TRUE)
  }
  return(lk)
}


initialization <- function(s1, s2, s3, m1, m2, m3){
  s1<-s1*m1
  s2<-s2*m2
  s3<-s3*m3
  p0ps=diag(rep(-s1,m1))
  p1ps=diag(rep(-s2,m2))
  p2ps=diag(rep(-s3,m3))
  for(j in 1:(m1-1)){
    p0ps[j,j+1]=-p0ps[j,j]
  }
  nups<-apply(p0ps,MARGIN = 1,sum)*(-1)
  for(j in 1:(m2-1)){
    p1ps[j,j+1]=-p1ps[j,j]
  }
  nu1ps<-apply(p1ps,MARGIN = 1,sum)*(-1)
  for(j in 1:(m3-1)){
    p2ps[j,j+1]=-p2ps[j,j]
  }
  nu2ps<-apply(p2ps,MARGIN = 1,sum)*(-1)
  return(list(p0ps, p1ps, p2ps, nups, nu1ps, nu2ps))
}

cross_validation <- function(s1_candidate, s2_candidate, s3_candidate, fold, original_data_list, stopping_criteria = 1e-2){
  original_covariate <- original_data_list[[1]]
  original_ts <- original_data_list[[2]]
  original_tau <- original_data_list[[3]]
  original_lags <- original_data_list[[4]]
  original_censor <- original_data_list[[5]]
  my_split<-sample(1:sample_size,sample_size)
  if((ceiling(sample_size/fold)*fold-sample_size)>0){
    my_split<-c(my_split,rep(1e5, (ceiling(sample_size/fold)*fold-sample_size) ))
  }
  my_split<-matrix(my_split,fold,length(my_split)/fold)
  
  CVL_array<-array(0,dim = c(length(s1_candidate),length(s2_candidate),length(s3_candidate)))
  increment_array<-array(0,dim = c(length(s1_candidate),length(s2_candidate),length(s3_candidate)))
  
  for(gt1 in 1:length(s1_candidate)){
    for(gt2 in 1:length(s2_candidate)){
      for(gt3 in 1:length(s3_candidate)){
        s1<-s1_candidate[gt1]*m1
        s2<-s2_candidate[gt2]*m2
        s3<-s3_candidate[gt3]*m3
        temp_CVL<-rep(0,fold)
        temp_increment<-rep(0,fold)
        for(my_fold in 1:fold){
          temp_aaa<-c(my_split[-my_fold,])
          temp_aaa<-temp_aaa[which(temp_aaa<1e4)]
          covariate<-original_covariate[temp_aaa,]
          ts<-original_ts[temp_aaa]
          tau<-original_tau[temp_aaa]
          lags<-original_lags[temp_aaa]
          censor<-original_censor[temp_aaa]
          sample_size<-nrow(covariate)
          
          
          ids=which((pmin(ts,tau)+lags)<censor)
          lags=lags[ids]
          observed<-pmin(ts,tau)[ids]
          delta<-as.numeric(pmin(ts,tau)[ids]==ts[ids])
          covariate_reported<-covariate[ids,]
          covariate_unreported<-covariate[-ids,]
          censor<-censor[-ids]
          
          
          p0ps=diag(rep(-s1,m1))
          p1ps=diag(rep(-s2,m2))
          p2ps=diag(rep(-s3,m3))
          for(j in 1:(m1-1)){
            p0ps[j,j+1]=-p0ps[j,j]
          }
          nups<-apply(p0ps,MARGIN = 1,sum)*(-1)
          for(j in 1:(m2-1)){
            p1ps[j,j+1]=-p1ps[j,j]
          }
          nu1ps<-apply(p1ps,MARGIN = 1,sum)*(-1)
          for(j in 1:(m3-1)){
            p2ps[j,j+1]=-p2ps[j,j]
          }
          nu2ps<-apply(p2ps,MARGIN = 1,sum)*(-1)
          mu2ps=max(abs(diag(p2ps)))+0.1
          
          old_pips<-rep(1/m1,m1)
          old_pi1ps<-rep(1/m2,m2)
          old_pi2ps<-rep(1/m3,m3)
          
          old_beta<-rep(0,covariate_dim)#true_beta
          old_alpha<-rep(0,covariate_dim)#true_alpha#
          
          
          ############
          ites<-0
          Rmax=100
          ptm<-proc.time()
          increment<-Inf
          ###EM
          while((increment>stopping_criteria)&(ites<1000)){
            updates=emcore(observed,delta,lags,censor,old_pips,old_pi1ps,old_pi2ps,old_beta,old_alpha,Rmax,p0ps,p1ps,p2ps,nups,nu1ps,nu2ps,covariate_reported,covariate_unreported)
            ### save the updated parameter estimates
            new_pips<-updates[[1]]
            new_pi1ps<-updates[[2]]
            new_pi2ps<-updates[[3]]
            new_beta<-updates[[4]]
            new_alpha<-updates[[5]]
            ### compute the increment in observed-data log-likelihood
            increment<-loglike(new_beta,new_alpha,new_pips,new_pi1ps,new_pi2ps,observed,delta,lags,censor,p0ps,p1ps,p2ps,nups,nu1ps,nu2ps,covariate_reported,covariate_unreported)-loglike(old_beta,old_alpha,old_pips,old_pi1ps,old_pi2ps,observed,delta,lags,censor,p0ps,p1ps,p2ps,nups,nu1ps,nu2ps,covariate_reported,covariate_unreported)
            while((increment<0)&(Rmax<1500)){
              Rmax<-Rmax+500
              # print(Rmax)
              updates=emcore(observed,delta,lags,censor,old_pips,old_pi1ps,old_pi2ps,old_beta,old_alpha,Rmax,p0ps,p1ps,p2ps,nups,nu1ps,nu2ps,covariate_reported,covariate_unreported)
              ### save the updated parameter estimates
              new_pips<-updates[[1]]
              new_pi1ps<-updates[[2]]
              new_pi2ps<-updates[[3]]
              new_beta<-updates[[4]]
              new_alpha<-updates[[5]]
              ### compute the increment in observed-data log-likelihood
              increment<-loglike(new_beta,new_alpha,new_pips,new_pi1ps,new_pi2ps,observed,delta,lags,censor,p0ps,p1ps,p2ps,nups,nu1ps,nu2ps,covariate_reported,covariate_unreported)-loglike(old_beta,old_alpha,old_pips,old_pi1ps,old_pi2ps,observed,delta,lags,censor,p0ps,p1ps,p2ps,nups,nu1ps,nu2ps,covariate_reported,covariate_unreported)
            }
            Rmax=100
            ites=ites+1
            old_beta<-new_beta
            old_alpha<-new_alpha
            old_pips<-new_pips
            old_pi1ps<-new_pi1ps
            old_pi2ps<-new_pi2ps
            # print(c(ites,increment))
          }
          
          temp_increment[fold]<-increment
          
          #########################################validation
          temp_aaa<-c(my_split[my_fold,])
          temp_aaa<-temp_aaa[which(temp_aaa<1e4)]
          
          covariate<-original_covariate[temp_aaa,]
          ts<-original_ts[temp_aaa]
          tau<-original_tau[temp_aaa]
          lags<-original_lags[temp_aaa]
          censor<-original_censor[temp_aaa]
          sample_size<-nrow(covariate)
          
          
          ids=which((pmin(ts,tau)+lags)<censor)
          lags=lags[ids]
          observed<-pmin(ts,tau)[ids]
          delta<-as.numeric(pmin(ts,tau)[ids]==ts[ids])
          covariate_reported<-covariate[ids,]
          covariate_unreported<-covariate[-ids,]
          censor<-censor[-ids]
          temp_CVL[fold]<-loglike(new_beta,new_alpha,new_pips,new_pi1ps,new_pi2ps,observed,delta,lags,censor,p0ps,p1ps,p2ps,nups,nu1ps,nu2ps,covariate_reported,covariate_unreported)
        }
        CVL_array[gt1,gt2,gt3]<-sum(temp_CVL)
        increment_array[gt1,gt2,gt3]<-min(temp_increment)
      }
    }
  }
  
  temp<-which(CVL_array==max(CVL_array),arr.ind = TRUE)
  temp<-temp[1,]
  
  s1<-s1_candidate[temp[1]]
  s2<-s2_candidate[temp[2]]
  s3<-s3_candidate[temp[3]]
  
  return(c(s1, s2, s3))
}

EM_algo <- function(observed, delta, lags, censor, initial_pips, initial_pi1ps, initial_pi2ps, initial_beta, initial_alpha, stopping_criteria = 1e-2, max_iter = 1000){
  old_pips <- initial_pips
  old_pi1ps <- initial_pi1ps
  old_pi2ps <- initial_pi2ps
  old_beta <- initial_beta
  old_alpha <- initial_alpha
  ites<-0
  Rmax=100
  ptm<-proc.time()
  increment<-Inf
  ###EM
  while((increment>stopping_criteria)&(ites<max_iter)){
    updates=emcore(observed,delta,lags,censor,old_pips,old_pi1ps,old_pi2ps,old_beta,old_alpha,Rmax,p0ps,p1ps,p2ps,nups,nu1ps,nu2ps,covariate_reported,covariate_unreported)
    ### save the updated parameter estimates
    new_pips<-updates[[1]]
    new_pi1ps<-updates[[2]]
    new_pi2ps<-updates[[3]]
    new_beta<-updates[[4]]
    new_alpha<-updates[[5]]
    ### compute the increment in observed-data log-likelihood
    increment<-loglike(new_beta,new_alpha,new_pips,new_pi1ps,new_pi2ps,observed,delta,lags,censor,p0ps,p1ps,p2ps,nups,nu1ps,nu2ps,covariate_reported,covariate_unreported)-loglike(old_beta,old_alpha,old_pips,old_pi1ps,old_pi2ps,observed,delta,lags,censor,p0ps,p1ps,p2ps,nups,nu1ps,nu2ps,covariate_reported,covariate_unreported)
    while((increment<0)&(Rmax<1500)){
      Rmax<-Rmax+500
      # print(Rmax)
      updates=emcore(observed,delta,lags,censor,old_pips,old_pi1ps,old_pi2ps,old_beta,old_alpha,Rmax,p0ps,p1ps,p2ps,nups,nu1ps,nu2ps,covariate_reported,covariate_unreported)
      ### save the updated parameter estimates
      new_pips<-updates[[1]]
      new_pi1ps<-updates[[2]]
      new_pi2ps<-updates[[3]]
      new_beta<-updates[[4]]
      new_alpha<-updates[[5]]
      ### compute the increment in observed-data log-likelihood
      increment<-loglike(new_beta,new_alpha,new_pips,new_pi1ps,new_pi2ps,observed,delta,lags,censor,p0ps,p1ps,p2ps,nups,nu1ps,nu2ps,covariate_reported,covariate_unreported)-loglike(old_beta,old_alpha,old_pips,old_pi1ps,old_pi2ps,observed,delta,lags,censor,p0ps,p1ps,p2ps,nups,nu1ps,nu2ps,covariate_reported,covariate_unreported)
    }
    Rmax=100
    ites=ites+1
    old_beta<-new_beta
    old_alpha<-new_alpha
    old_pips<-new_pips
    old_pi1ps<-new_pi1ps
    old_pi2ps<-new_pi2ps
    cat(paste("EM iteration:", ites, "   increment in log-likelihood", increment, "\n"))
  }
  EM_time<-as.numeric((proc.time()-ptm)[3])
  EM_iteration<-ites
  ##########estimate
  return(list(new_beta, new_alpha, new_pips, new_pi1ps, new_pi2ps))
}


confidence_interval <- function(original_data_list, p0ps, p1ps, p2ps, nups, nu1ps, nu2ps, EM_res, stopping_criteria = 1e-2){
  covariate_dim <- 5
  covariate <- original_data_list[[1]]
  ts <- original_data_list[[2]]
  tau <- original_data_list[[3]]
  lags <- original_data_list[[4]]
  censor <- original_data_list[[5]]
  old_lags <- lags
  old_ts <- ts
  old_tau <- tau
  old_censor <- censor
  sample_size <- length(censor);#sample size
  
  status<-rep(0,sample_size)
  for(i in 1:sample_size){
    if(min(ts[i],tau[i])+lags[i]<censor[i]){
      if(ts[i]<tau[i]){
        status[i]<-1
      }
    }
  }
  
  ids=which((pmin(ts,tau)+lags)<censor)
  report_status <- as.numeric((pmin(ts,tau)+lags)<censor)
  lags=lags[ids]
  observed<-pmin(ts,tau)[ids]
  delta<-as.numeric(pmin(ts,tau)[ids]==ts[ids])
  covariate_reported<-covariate[ids,]
  covariate_unreported<-covariate[-ids,]
  ncensor=sample_size-length(lags);
  censor<-censor[-ids]
  
  m1<-ceiling((sample_size)^(1/4))
  m2<-ceiling((sample_size)^(1/4))
  m3<-ceiling((sample_size)^(1/4))
  
  s1 <- abs(diag(p0ps))[1]
  s2 <- abs(diag(p1ps))[1]
  s3 <- abs(diag(p2ps))[1]
  mu2ps=max(abs(diag(p2ps)))+0.1
  
  MLE_beta <- EM_res[[1]]
  MLE_alpha <- EM_res[[2]]
  pi_MLE1 <- EM_res[[3]]
  pi_MLE2 <- EM_res[[4]]
  pi_MLE3 <- EM_res[[5]]
  
  T_MLE1 <- p0ps
  T_MLE2 <- p1ps
  T_MLE3 <- p2ps
  
  xi_MLE1 <- apply(T_MLE1,MARGIN = 1,sum)*(-1)
  xi_MLE2 <- apply(T_MLE2,MARGIN = 1,sum)*(-1)
  xi_MLE3 <- apply(T_MLE3,MARGIN = 1,sum)*(-1)
  
  # this should correspond to the rdata on Nov 8; or change 2.5 to 2
  pertubation_constant<-rep(2.5,10)
  pertubation_constant[4]<-10
  
  emcore_fix_beta=function(y,delta,y1,censor,pips,pi1ps,pi2ps,est_beta,est_alpha){
    old_y<-y
    old_y1<-y1
    m=length(pips)
    Hs<-matrix(0,m,m)
    newHs<-matrix(0,m,m)
    long_Hs<-array(0,dim = c(length(y),m,m))
    long_Hs_transit<-array(0,dim = c(length(y),m,m))
    Bs<-rep(0,m)
    for(i in 1:length(y)){
      mups<-max(abs(diag(exp(-sum(covariate_reported[i,]*est_beta))*p0ps)))+0.1
      re_D0=exp(-sum(covariate_reported[i,]*est_beta))*p0ps/mups+diag(m)
      if(delta[i]==1){
        vec1<-nups*exp(-sum(covariate_reported[i,]*est_beta))
        hc<-pips*expAtv(A=exp(-sum(covariate_reported[i,]*est_beta))*p0ps,v=vec1,t=y[i])$eAtv
        sumhc<-sum(hc)
        Bs<-Bs+hc/sumhc
      }else{
        vec1<-rep(1,m)
        hc<-pips*expAtv(A=exp(-sum(covariate_reported[i,]*est_beta))*p0ps,v=vec1,t=y[i])$eAtv
        sumhc<-sum(hc)
        Bs<-Bs+hc/sumhc
      }
      CHks=matrix(0,m,m)
      calps=matrix(vec1,m,Rmax+1)
      for(j in 2:(Rmax+1)) calps[,j]=re_D0%*%calps[,j-1]
      cbets=matrix(dpois(Rmax+1,mups*y[i])*pips,Rmax+1,m,byrow=T)
      for(j in Rmax:1){
        cbets[j,]=cbets[j+1,]%*%re_D0+dpois(j,mups*y[i])*pips
      }
      for(j in 1:(Rmax+1)){
        CHks=CHks+calps[,j]%*%t(cbets[j,])
      }
      CHks<-t(CHks)/sumhc/mups
      long_Hs[i,,]<-CHks
      Hs<-Hs+CHks
      newHs<-newHs+CHks*exp(-sum(covariate_reported[i,]*est_beta))*p0ps
      long_Hs_transit[i,,]<-CHks*exp(-sum(covariate_reported[i,]*est_beta))*p0ps
    }
    
    ###consider the drop-out time
    m1=length(pi1ps)
    Hs1<-matrix(0,m1,m1)
    newHs1<-matrix(0,m1,m1)
    long_Hs1<-array(0,dim = c(length(y),m1,m1))
    long_Hs1_transit<-array(0,dim = c(length(y),m1,m1))
    Bs1<-rep(0,m1)
    for(i in 1:length(y)){
      mu1ps<-max(abs(diag(exp(-sum(covariate_reported[i,]*est_alpha))*p1ps)))+0.1
      re_D0=exp(-sum(covariate_reported[i,]*est_alpha))*p1ps/mu1ps+diag(m1)
      if(delta[i]==0){
        vec1<-nu1ps*exp(-sum(covariate_reported[i,]*est_alpha))
        hc<-pi1ps*expAtv(A=exp(-sum(covariate_reported[i,]*est_alpha))*p1ps,v=vec1,t=y[i])$eAtv
        sumhc<-sum(hc)
        Bs1<-Bs1+hc/sumhc
      }else{
        vec1<-rep(1,m1)
        hc<-pi1ps*expAtv(A=exp(-sum(covariate_reported[i,]*est_alpha))*p1ps,v=vec1,t=y[i])$eAtv
        sumhc<-sum(hc)
        Bs1<-Bs1+hc/sumhc
      }
      CHks1=matrix(0,m1,m1)
      calps=matrix(vec1,m1,Rmax+1)
      for(j in 2:(Rmax+1)) calps[,j]=re_D0%*%calps[,j-1]
      cbets=matrix(dpois(Rmax+1,mu1ps*y[i])*pi1ps,Rmax+1,m1,byrow=T)
      for(j in Rmax:1){
        cbets[j,]=cbets[j+1,]%*%re_D0+dpois(j,mu1ps*y[i])*pi1ps
      }
      for(j in 1:(Rmax+1)){
        CHks1=CHks1+calps[,j]%*%t(cbets[j,])
      }
      CHks1<-t(CHks1)/sumhc/mu1ps
      long_Hs1[i,,]<-CHks1
      Hs1<-Hs1+CHks1
      newHs1<-newHs1+CHks1*exp(-sum(covariate_reported[i,]*est_alpha))*p1ps
      long_Hs1_transit[i,,]<-CHks1*exp(-sum(covariate_reported[i,]*est_alpha))*p1ps
    }
    
    ###afterwards, we consider the observed reporting delays
    y1<-sort(y1)
    dy2<-c(y1[1],diff(y1))
    m2=length(pi2ps)
    ss<-length(dy2)
    re_p02=p2ps/mu2ps+diag(m2)
    fks2=matrix(expAtv(t(p2ps),pi2ps,dy2[1])$eAtv,m2,ss)
    bks2=matrix(expAtv(p2ps,nu2ps,dy2[1])$eAtv,m2,ss)
    for(i in 2:ss){
      fks2[,i]=expAtv(t(p2ps),fks2[,i-1],dy2[i])$eAtv
      bks2[,i]=expAtv(p2ps,bks2[,i-1],dy2[i])$eAtv
    }
    cks2=matrix(pi2ps/sum(pi2ps*bks2[,ss]),ss,m2,byrow=T)
    for(i in (ss-1):1){
      cks2[i,]=expAtv(t(p2ps),cks2[i+1,],dy2[i+1])$eAtv+pi2ps/sum(pi2ps*bks2[,i])
    }
    bks2=cbind(nu2ps,bks2)
    Hs2=matrix(0,m2,m2)
    for(i in 1:ss){
      Hks2=matrix(0,m2,m2)
      alps2=matrix(bks2[,i],m2,Rmax+1)
      for(j in 2:(Rmax+1)) alps2[,j]=re_p02%*%alps2[,j-1]
      bets2=matrix(dpois(Rmax+1,mu2ps*dy2[i])*cks2[i,],Rmax+1,m2,byrow=T)
      for(j in Rmax:1){
        bets2[j,]=bets2[j+1,]%*%re_p02+dpois(j,mu2ps*dy2[i])*cks2[i,]
      }
      for(j in 1:(Rmax+1)){
        Hks2=Hks2+alps2[,j]%*%t(bets2[j,])
      }
      Hs2=Hs2+Hks2/mu2ps
    }
    Bs2=rep(0,m2)
    for(i in 1:ss){
      Bs2=Bs2+pi2ps*bks2[,i+1]/sum(pi2ps*bks2[,i+1])
    }
    Hs2<-t(Hs2)
    
    ###for missing data
    mm1<-m*m1+m2
    Hs3<-matrix(0,mm1,mm1)
    newCHks<-matrix(0,mm1,mm1)
    long_Hs3<-array(0,dim = c(length(censor),mm1,mm1))
    long_Hs3_transit<-array(0,dim = c(length(censor),mm1,mm1))
    Bs3<-rep(0,mm1)
    for(sub in 1:length(censor)){
      new_pips<-c(kronecker(pips,pi1ps),rep(0,m2))
      temp_nups<-rep(0,m*m1)
      for(i in 1:m){
        temp_nups[i*m1]<-temp_nups[i*m1]+exp(-sum(covariate_unreported[sub,]*est_alpha))*nu1ps[m1]
        # temp_nups[i*m1-1]<-temp_nups[i*m1-1]+exp(-sum(covariate_unreported[sub,]*est_alpha))*nu1ps[m1-1]
      }
      for(i in 1:m1){
        # temp_nups[(m-2)*m1+i]<-temp_nups[(m-2)*m1+i]+exp(-sum(covariate_unreported[sub,]*est_beta))*nups[m-1]
        temp_nups[(m-1)*m1+i]<-temp_nups[(m-1)*m1+i]+exp(-sum(covariate_unreported[sub,]*est_beta))*nups[m]
      }
      new_p0ps<-cbind((kronecker(exp(-sum(covariate_unreported[sub,]*est_beta))*p0ps,diag(rep(1,m1)))+kronecker(diag(rep(1,m)),exp(-sum(covariate_unreported[sub,]*est_alpha))*p1ps)),temp_nups%*%t(pi2ps))
      new_p0ps<-rbind(new_p0ps,cbind(matrix(0,m2,m*m1),p2ps))
      new_nups<-apply(new_p0ps,MARGIN = 1,sum)*(-1)
      
      hc=new_pips*expAtv(new_p0ps,rep(1,mm1),censor[sub])$eAtv
      sumhc<-sum(hc)
      Bs3<-Bs3+hc/sumhc
      mmups<-max(abs(diag(new_p0ps)))+0.1
      CHks=matrix(0,mm1,mm1)
      re_D0=new_p0ps/mmups+diag(mm1)
      es=rep(1,mm1)
      
      calps=matrix(es,mm1,Rmax+1)
      for(j in 2:(Rmax+1)) calps[,j]=re_D0%*%calps[,j-1]
      cbets=matrix(dpois(Rmax+1,mmups*censor[sub])*new_pips,Rmax+1,mm1,byrow=T)
      for(j in Rmax:1){
        cbets[j,]=cbets[j+1,]%*%re_D0+dpois(j,mmups*censor[sub])*new_pips
      }
      for(j in 1:(Rmax+1)){
        CHks=CHks+calps[,j]%*%t(cbets[j,])
      }
      CHks=t(CHks)/sumhc/mmups
      long_Hs3[sub,,]<-CHks
      Hs3<-Hs3+CHks
      newCHks<-newCHks+CHks*new_p0ps
      long_Hs3_transit[sub,,]<-CHks*new_p0ps
    }
    CHks<-Hs3
    
    
    
    
    ####################################M Step
    ###initial prob vector of lifetime T
    temp_vec1<-Bs
    for(i in 1:m){
      temp_vec1[i]<-temp_vec1[i]+sum(Bs3[((i-1)*m1+1):(i*m1)])
    }
    output_pips<-temp_vec1/sum(temp_vec1)
    ###initial prob vector of censoring time tau
    temp_vec2<-Bs1
    for(i in 1:m1){
      for(j in 1:m){
        temp_vec2[i]<-temp_vec2[i]+Bs3[(j-1)*m1+i]
      }
    }
    output_pi1ps<-temp_vec2/sum(temp_vec2)
    ###initial prob vector of reporting delay X
    # new_Hs1<-Hs1*p1ps
    new_Hs2<-Hs2*p2ps
    temp_vec3<-Bs2
    for(i in 1:m2){
      temp_vec3[i]<-temp_vec3[i]+sum(newCHks[1:(m*m1),(m*m1+i)])
    }
    output_pi2ps<-temp_vec3/sum(temp_vec3)
    
    new_beta<-est_beta
    new_alpha<-est_alpha
    ####results
    list(output_pips,output_pi1ps,output_pi2ps)
  }
  
  loglike <- function(beta,myalpha,pips,pi1ps,pi2ps,y,delta,y1,censor){
    m=length(pips);m1=length(pi1ps);m2=length(pi2ps)
    mm1<-m*m1+m2
    lk<-sum(dphase(y1,ph=ph(alpha=pi2ps,Q=p2ps,xi=nu2ps),log=T))
    for(sub in 1:length(y)){
      if(delta[sub]==1){
        lk<-lk+dphase(y[sub],ph=ph(alpha = pips,Q=p0ps*exp(-sum(beta*covariate_reported[sub,])),xi=nups*exp(-sum(beta*covariate_reported[sub,]))),log = TRUE)
        lk<-lk+pphase(y[sub],ph=ph(alpha = pi1ps,Q=p1ps*exp(-sum(myalpha*covariate_reported[sub,])),xi=nu1ps*exp(-sum(myalpha*covariate_reported[sub,]))),lower.tail = FALSE,log.p=TRUE)
      }else{
        lk<-lk+pphase(y[sub],ph=ph(alpha = pips,Q=p0ps*exp(-sum(beta*covariate_reported[sub,])),xi=nups*exp(-sum(beta*covariate_reported[sub,]))),lower.tail = FALSE,log.p=TRUE)
        lk<-lk+dphase(y[sub],ph=ph(alpha = pi1ps,Q=p1ps*exp(-sum(myalpha*covariate_reported[sub,])),xi=nu1ps*exp(-sum(myalpha*covariate_reported[sub,]))),log = TRUE)
      }
    }
    for(sub in 1:length(censor)){
      new_pips<-c(kronecker(pips,pi1ps),rep(0,m2))
      temp_nups<-rep(0,m*m1)
      for(i in 1:m){
        temp_nups[i*m1]<-temp_nups[i*m1]+exp(-sum(covariate_unreported[sub,]*myalpha))*nu1ps[m1]
      }
      for(i in 1:m1){
        temp_nups[(m-1)*m1+i]<-temp_nups[(m-1)*m1+i]+exp(-sum(covariate_unreported[sub,]*beta))*nups[m]
      }
      new_p0ps<-cbind((kronecker(exp(-sum(covariate_unreported[sub,]*beta))*p0ps,diag(rep(1,m1)))+kronecker(diag(rep(1,m)),exp(-sum(covariate_unreported[sub,]*myalpha))*p1ps)),temp_nups%*%t(pi2ps))
      new_p0ps<-rbind(new_p0ps,cbind(matrix(0,m2,m*m1),p2ps))
      new_nups<-apply(new_p0ps,MARGIN = 1,sum)*(-1)
      lk<-lk+pphase(censor[sub],ph=ph(alpha=new_pips,Q=new_p0ps,xi=new_nups),lower.tail = FALSE,log.p=TRUE)
    }
    return(lk)
  }
  
  profile_nonparametric_part<-as.list(rep(0,2*covariate_dim))
  
  
  
  for(ppp in 1:(2*covariate_dim)){
    old_pips<-pi_MLE1
    old_pi1ps<-pi_MLE2
    old_pi2ps<-pi_MLE3
    
    temp_coef_MLE<-c(MLE_beta,MLE_alpha)
    temp_coef_MLE[ppp]<-temp_coef_MLE[ppp]+pertubation_constant[ppp]
    old_beta<-temp_coef_MLE[1:covariate_dim]
    old_alpha<-temp_coef_MLE[(covariate_dim+1):(2*covariate_dim)]
    
    
    ############
    ites<-0
    Rmax=100
    ptm<-proc.time()
    increment<-Inf
    ###EM
    while((increment>stopping_criteria)&(ites<1000)){
      updates=emcore_fix_beta(observed,delta,lags,censor,old_pips,old_pi1ps,old_pi2ps,old_beta,old_alpha)
      ### save the updated parameter estimates
      new_pips<-updates[[1]]
      new_pi1ps<-updates[[2]]
      new_pi2ps<-updates[[3]]
      ### compute the increment in observed-data log-likelihood
      increment<-loglike(old_beta,old_alpha,new_pips,new_pi1ps,new_pi2ps,observed,delta,lags,censor)-loglike(old_beta,old_alpha,old_pips,old_pi1ps,old_pi2ps,observed,delta,lags,censor)
      while((increment<0)&(Rmax<1500)){
        Rmax<-Rmax+500
        updates=emcore_fix_beta(observed,delta,lags,censor,old_pips,old_pi1ps,old_pi2ps,old_beta,old_alpha)
        ### save the updated parameter estimates
        new_pips<-updates[[1]]
        new_pi1ps<-updates[[2]]
        new_pi2ps<-updates[[3]]
        ### compute the increment in observed-data log-likelihood
        increment<-loglike(old_beta,old_alpha,new_pips,new_pi1ps,new_pi2ps,observed,delta,lags,censor)-loglike(old_beta,old_alpha,old_pips,old_pi1ps,old_pi2ps,observed,delta,lags,censor)
      }
      Rmax=100
      ites=ites+1
      old_pips<-new_pips
      old_pi1ps<-new_pi1ps
      old_pi2ps<-new_pi2ps
      # print(c(ites,increment))
    }
    ##########estimate
    pi_est1=new_pips;pi_est2=new_pi1ps;pi_est3=new_pi2ps
    xi_est1=nups;xi_est2=nu1ps;xi_est3=nu2ps
    T_est1=p0ps;T_est2=p1ps;T_est3=p2ps
    profile_nonparametric_part[[ppp]]<-list(old_beta,old_alpha,pi_est1,pi_est2,pi_est3)
  }
  
  individual_loglike<-function(beta,myalpha,pips,pi1ps,pi2ps,index){
    res<-0
    m=length(pips);m1=length(pi1ps);m2=length(pi2ps) 
    mm1<-m*m1+m2
    if(report_status[index]==1){
      if(status[index]==1){
        res<-res+dphase(old_ts[index],ph=ph(alpha = pips,Q=p0ps*exp(-sum(beta*covariate[index,])),xi=nups*exp(-sum(beta*covariate[index,]))),log = TRUE)
        res<-res+dphase(old_lags[index],ph=ph(alpha=pi2ps,Q=p2ps,xi=nu2ps),log=TRUE)
        res<-res+pphase(old_ts[index],ph=ph(alpha = pi1ps,Q=p1ps*exp(-sum(myalpha*covariate[index,])),xi=nu1ps*exp(-sum(myalpha*covariate[index,]))),lower.tail = FALSE,log.p = TRUE)
      }else{
        res<-res+pphase(old_tau[index],ph=ph(alpha = pips,Q=p0ps*exp(-sum(beta*covariate[index,])),xi=nups*exp(-sum(beta*covariate[index,]))),lower.tail = FALSE,log.p = TRUE)
        res<-res+dphase(old_lags[index],ph=ph(alpha=pi2ps,Q=p2ps,xi=nu2ps),log=TRUE)
        res<-res+dphase(old_tau[index],ph=ph(alpha = pi1ps,Q=p1ps*exp(-sum(myalpha*covariate[index,])),xi=nu1ps*exp(-sum(myalpha*covariate[index,]))),log = TRUE)
      }
    }else{
      new_pips<-c(kronecker(pips,pi1ps),rep(0,m2))
      temp_nups<-rep(0,m*m1)
      for(i in 1:m){
        temp_nups[i*m1]<-temp_nups[i*m1]+exp(-sum(covariate[index,]*myalpha))*nu1ps[m1]
      }
      for(i in 1:m1){
        temp_nups[(m-1)*m1+i]<-temp_nups[(m-1)*m1+i]+exp(-sum(covariate[index,]*beta))*nups[m]
      }
      new_p0ps<-cbind((kronecker(exp(-sum(covariate[index,]*beta))*p0ps,diag(rep(1,m1)))+kronecker(diag(rep(1,m)),exp(-sum(covariate[index,]*myalpha))*p1ps)),temp_nups%*%t(pi2ps))
      new_p0ps<-rbind(new_p0ps,cbind(matrix(0,m2,m*m1),p2ps))
      new_nups<-apply(new_p0ps,MARGIN = 1,sum)*(-1)
      res<-res+pphase(old_censor[index],ph=ph(alpha=new_pips,Q=new_p0ps,xi=new_nups),lower.tail = FALSE,log.p=TRUE)
    }
    return(res)
  }
  cov_mat_est<-matrix(0,(2*covariate_dim),(2*covariate_dim))
  for(i in 1:sample_size){
    term1<-individual_loglike(beta=MLE_beta,myalpha=MLE_alpha,pips=pi_MLE1,pi1ps=pi_MLE2,pi2ps=pi_MLE3,index=i)
    for(III in 1:(2*covariate_dim)){
      for(JJJ in III:(2*covariate_dim)){
        if(III==JJJ){
          term2<-individual_loglike(beta=profile_nonparametric_part[[III]][[1]],myalpha=profile_nonparametric_part[[III]][[2]],pips=profile_nonparametric_part[[III]][[3]],pi1ps=profile_nonparametric_part[[III]][[4]],pi2ps=profile_nonparametric_part[[III]][[5]],index=i)
          cov_mat_est[III,JJJ]<-cov_mat_est[III,JJJ]+((term2-term1)/(pertubation_constant[III]))^2
        }else{
          term2<-individual_loglike(beta=profile_nonparametric_part[[III]][[1]],myalpha=profile_nonparametric_part[[III]][[2]],pips=profile_nonparametric_part[[III]][[3]],pi1ps=profile_nonparametric_part[[III]][[4]],pi2ps=profile_nonparametric_part[[III]][[5]],index=i)
          term3<-individual_loglike(beta=profile_nonparametric_part[[JJJ]][[1]],myalpha=profile_nonparametric_part[[JJJ]][[2]],pips=profile_nonparametric_part[[JJJ]][[3]],pi1ps=profile_nonparametric_part[[JJJ]][[4]],pi2ps=profile_nonparametric_part[[JJJ]][[5]],index=i)
          cov_mat_est[III,JJJ]<-cov_mat_est[III,JJJ]+((term2-term1)/(pertubation_constant[III]))*((term3-term1)/(pertubation_constant[JJJ]))
          cov_mat_est[JJJ,III]<-cov_mat_est[JJJ,III]+((term2-term1)/(pertubation_constant[III]))*((term3-term1)/(pertubation_constant[JJJ]))
        }
      }
    }
  }
  cov_mat_est<-cov_mat_est/sample_size
  cov_mat_est<-solve(cov_mat_est)/sample_size 
  est_coef_MLE<-c(MLE_beta,MLE_alpha)
  CI_lower<-est_coef_MLE-sqrt(diag(cov_mat_est))*qnorm(0.975)
  CI_upper<-est_coef_MLE+sqrt(diag(cov_mat_est))*qnorm(0.975)
  return(list(CI_lower, CI_upper))
}