FITtestS<-function(GXp=pnorm(300),K=20) {
  #pnorm((Y-E(Y|X))/Var(Y|X)^0.5)
  # GXP fonction de repartition 
  # K= NB de classes
  N<-length(GXp)
  NR<-5000 
  #la taille d'ech pour Polya
  TTP<-rep(NA,6) #SORTIES
  ttpit<-cvm.test(GXp,"punif")
  TTP[1]<-ttpit[[1]]
  TTP[2]<-ttpit[[2]]
  Hpit<-hist(GXp,seq(0,1,1/K),plot=FALSE)
  chit<-sum(((Hpit[[2]]-N/K)^2)/(N/K))
  Pdiri<-bmixture::rdirichlet(NR, alpha =rep((N+1)/K,K) )
  Mnom<-matrix(NA,NR,K)
  for(j in 1:NR){
    Mnom[j,]<-as.vector(rmultinom(1,N,Pdiri[j,]))
  }
  Postpit<-apply(((Mnom-N/K)^2)/(N/K),1,"sum")
  TTP[3]<-chit
  TTP[4]<-length(which(Postpit>=chit))/NR
  liste = list(cvm=TTP[1], p.cvm= TTP[2],chi2= TTP[3], p.bayes.chi2=TTP[4],
               comptages= Hpit[[2]],matrix.simul.polya = Mnom )
 # return(list(Hpit[[2]],TTP,Mnom))
  return(liste)
}

#-------  fin de FITtestS et début de FITvalid -----------
FITvalid<-function(X=rnorm(300),mpred =0, sdpred=1,mcal=0,sdcal=1,K=20) {
  
  #X cible Y, mpred = E(Y|X), sdpred= Var(Y|X)^0.5
  # mcal moyenne de la loi normale de "calage" Y,
  # sdcal= ecart-type 
  
  N<-length(X)
  brk0<-qnorm((1:K-1)/K,mcal,sdcal)
  # classe égale proba sur calage
  brke<-c(max(min(X),brk0[1]),brk0[2:K],max(X))
  Gpred<-pnorm(X,mpred,sdpred)
  #F rep predictive
  Gcal<-pnorm(X,mcal,sdcal)
  #F rep calage
  ybr<-diff(pnorm(brke,mcal,sdcal))
  # limites de classes de calage
  MY<-matrix(NA,N,K)
  # MY variables individuelles classées (1 ds la classe ou la cible est apparue pour chaque jour)
  ITQ<-rep(NA,N)
  #ITQ intervalles interquartiles
  Ppred<-matrix(NA,N,K)
  # Ppred probabilité prédictive que la cible tombe dans une classe donnée
  for(j in 1:N){
    HG<-hist(Gpred[j],brke,plot=FALSE)
    Yj<-HG[[2]]
    aj<-diff(pnorm(brke,mpred[j],sdpred[j])) 
    Ppred[j,]<-aj
    MY[j,]<-Yj # vecteur de présence de la cible dans une classe
    ITQ[j]<-sdpred[j]*(qnorm(0.75,0,1)-qnorm(0.25,0,1)) # intervalle interquartile
  }
  BR<-mean(apply((Ppred-MY)^2,1,"sum"))
  ITQp<-mean(ITQ) # moyenne interquartiles
  ITQc<-sd(X)*(qnorm(0.75,0,1)-qnorm(0.25,0,1)) # moyenne interquartile du calage
  crpsX<-scoringRules::crps_norm(X,location = mpred, scale = sdpred)
   
  liste= list(ScoreCRPS=-mean(crpsX),
              CVScoreCRPS=sd(crpsX)/mean(crpsX),
              ScoreIgnorance= mean(scoringRules::logs_norm(X,location = mpred, scale = sdpred)),
              BrierSym=BR,
              InterQ=ITQp,
              CVInterQ=sd(ITQ)/mean(ITQ),
              InterQcal=ITQc,
              corr=cor(X,mpred))
  liste = map_dbl(liste,signif,digits=3)
  return(liste)
}

dessinePIT = function(Y=rnorm(100),mean=0, sd=1, K=20){
  n=length(Y)
  abet<-1/K + n/K
  bbet<- (K-1)/K + (K-1) * n/K
  lim_sup<-qbeta(0.95, abet,bbet)*n
  lim_inf<-qbeta(0.05, abet,bbet)*n
  comptages=graphics::hist((Y-mean)/sd,breaks = qnorm(seq(0,1,length.out=K+1)),
                           plot = FALSE)$counts
  p<-tibble(freq=comptages,x=(1:K)/K) %>% ggplot(aes(x=x,y=freq))+
    geom_bar(stat="identity", fill="blue")+
    geom_hline(aes(yintercept=lim_sup), col="red")+
    geom_hline(aes(yintercept=lim_inf), col="red")+
    geom_hline(aes(yintercept=n/K))
}
