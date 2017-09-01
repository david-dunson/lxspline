
lxfit  <- function(x,y,nsamps = 50000, nburn = 10000, a_list = c(1,1,1.1,1.2,1.3,1.7,2,2.6,3.5,5,9,12,24,30) ){
  #make sure the y variable is a matrix
  t_y = as.matrix(y,ncol=1)
  #center the x variable to be between 0 and 1
  t_x = as.matrix((x-min(x))/(max(x)-min(x)),ncol=1 )
  if (sum(dim(t_x)==dim(t_y))!=2){
    stop("X and Y vectors must be of the same dimension")
  }
  
  #error check the JK quantity making sure 1 is in the list and
  #it is sorted
  a_list = sort(a_list)
  if (a_list[1] != 1){
    stop("The minimum value in the anealling list must be 1")
  }
  
  if (length(a_list)==1){
    stop("Anealling list must have more than one element")
  }
  
  if (nsamps < nburn){
    stop("Number of samples must be greater than the number of burnin samples.")
  }
  
  if (nsamps < 0){
    stop("Number of samples must be greater than zero.")
  }
  
  if (nburn < 0){
    warning("Burn in samples less than zero. Setting burn in samples to zero.")
    nburn = 1; 
  }
  
  return (.LXsample(t_x,t_y,a_list,nsamps,nburn))
  
}


#.LXsample This function 
.LXsample <- function(t,y,JK,nsamps,nburn){
  #create the initial tree list
  a1 = 0; a2 = 0; 
  m.TREE = matrix(-1,nrow=6,3)
  m.TREE[1,1:3] = as.matrix(c(0,0.5,1))
  m.TREE[2,1:3] = as.matrix(c(0,0.25,0))
  m.TREE[3,1:3] = as.matrix(c(0,0.75,0))
  m.TREE[4,1:3] = as.matrix(c(0,0.25,0)) 
  m.TREE[5,1:3] = as.matrix(c(0, 1, 0))
  m.TREE[6,1:3] = as.matrix(c(0, 1, 0))
  CBX = shapespline2(m.TREE[1,]-.5, t-0.5,2)
  
  ONES = matrix(1,nrow=length(t),ncol=1)
  X = cbind(ONES,100*(CBX[,,3] - (a1+a2)*CBX[,,2]+ (a1*a2)*CBX[,,1]));
  BETAS = matrix(0,ncol(X),ncol=1)
   
  ptree = 0.5
  tau = 0.25 
  x.Pred <- matrix(0,nrow=length(y),ncol=nsamps) #mean prediction of the data points
  lamT = 1
  intLAM = 0.001
  p = 0.01
  
  BBETAS = vector("list",length(JK));
  BTAUS  = vector("list",length(JK)); 
  BA     = vector("list",length(JK)); 
  BX     = vector("list",length(JK));
  BCBX   = vector("list",length(JK)); 
  BTREE  = vector("list",length(JK)); 
  BCLAM  = vector("list",length(JK));
  BP     = vector("list",length(JK));
  
  for (i in 1:length(JK)){
    BBETAS[[i]] = BETAS; 
    BTAUS[[i]] = tau; 
    BCBX[[i]] = CBX; 
    if (i ==2) {a1 = -0.5; a2 = -0.5}
    if (i ==3) {a1 = 0.5; a2 = 0.5}
    BA[[i]] = c(a1,a2); 
    
    BX[[i]] = X; 
    BTREE[[i]] = m.TREE; 
    BCLAM[[i]] = lamT; 
    BP[[i]]    = p; 
    
  }
  
  MT  <- 1:(length(JK)-1)
  ADJ <- matrix(0,nrow=length(JK)-1,ncol=2)
  ADJ[,1] =  MT; ADJ[,2] =  MT+1; 
  ADJ = rbind(ADJ,c(1,3))
  ADJ = rbind(ADJ,c(1,4))
  ADJ = rbind(ADJ,c(2,4))
  ADJ = rbind(ADJ,c(2,5))
  
  
  rr = 1
  
  
  sa1      =  rep(0,nsamps)
  sa2      =  rep(0,nsamps)
  n_knots  =  rep(0,nsamps)
  h_lambda = rep(0,nsamps)
  h_p      = rep(0,nsamps)
  model_sample     = vector("list",nsamps);
  beta_sample      = vector("list",nsamps);
  knot_sample      = vector("list",nsamps); 
  tau_sample       = 1:nsamps
  
  
  pb <- txtProgressBar(min=0, max=nsamps, style=3)
  
  
  for (i in 1:nsamps){
    
    setTxtProgressBar(pb,i)    
  
    if (runif(1) < 0.5){ # simple update 
      
      jj <- sample(1:length(JK),1)
      KK = JK[jj];
      BETAS = BBETAS[[jj]]; 
      tau = BTAUS[[jj]];
      CBX = BCBX[[jj]];  
      A = BA[[jj]]; a1 = A[1]; a2 = A[2]; 
      X = BX[[jj]]; 
      m.TREE = BTREE[[jj]]; 
      lamT = BCLAM[[jj]]; 
      
      
      #################################
      ##HYPER PRIOR LAMBDA
      betaT = BETAS[2:length(BETAS)]
      betaT = betaT[betaT > 0]
      
      ta = length(betaT) + .2
      tb = sum(betaT) + 2
      tg_ub = pgamma(1e-5,ta,tb)
      lamT = qgamma(runif(1,tg_ub,1),ta,tb)
      
      #other hyper prior on probability of a zero
      betaT = BETAS[2:length(BETAS)]
      V = betaT == 0
      t.n = length(betaT)
      ta = 2 + sum(V)
      tb = 18 + t.n - sum(V)
      p = rbeta(1,ta,tb)
      
      #################################
      X = cbind(ONES,100*(CBX[,,3] - (a1+a2)*CBX[,,2]+(a1*a2)*CBX[,,1]));
      ###############################################################
      LAM = lamT
      BETAS = sampleBetas(y, X, BETAS,LAM, intLAM, p,tau/KK)
   
      ################################################################
      
      tB <- 0.5*t(X%*%BETAS-y)%*%(X%*%BETAS-y)/KK+1
      tA <- 0.5*length(y)/KK+ 1
      tau <- as.vector((1/tB)*rgamma(1,tA,1))
      
      ystar =  y - cbind(ONES,100*(CBX[,,3] - (a2)*CBX[,,2]))%*%BETAS;
      w = 100*(-CBX[,,2]+(a2)*CBX[,,1])%*%BETAS[-1]
      a1  =  sampA(ystar,w,tau/KK,0,1,-0.5,0.5)
      
      
      X = cbind(ONES,100*(CBX[,,3] - (a1+a2)*CBX[,,2]+(a1*a2)*CBX[,,1]));
      
      ystar =  y - cbind(ONES,100*(CBX[,,3] - (a1)*CBX[,,2]))%*%BETAS; 
      w = 100*(-CBX[,,2]+(a1)*CBX[,,1])%*%BETAS[-1]
      a2  =  sampA(ystar,w,tau/KK,0,1,-0.5,0.5)
      
      X = cbind(ONES,100*(CBX[,,3] - (a1+a2)*CBX[,,2]+(a1*a2)*CBX[,,1]));
      
      L = PROP.N.TREE(m.TREE);
      
      POS = L$pos;   
      
      
      if (L$INSERT==TRUE){
        
        NCBX = shapesplineInsert(as.matrix(L$m.TREE[1,]-0.5,nrow=1), as.matrix(t,ncol=1)-0.5, 0.1,2,CBX,L$pos)
       
        NX = cbind(ONES,100*(NCBX[,,3] - (a1+a2)*NCBX[,,2]+(a1*a2)*NCBX[,,1]));
        
        NBETAS = matrix(0,nrow=ncol(NX),ncol=1); 
        
        T_IDX <- (POS-1):(POS+2)+1;
        NBETAS[-T_IDX] = BETAS[-((POS-1):(POS+1)+1)]; 
        yStar = y - NX[,-T_IDX]%*%NBETAS[-T_IDX];
        LAM = 0*abs(L$m.TREE[1,]+10)*sqrt(((L$m.TREE[1,]-a1+10))^2*(L$m.TREE[1,]-a2+10)^2)+lamT; 
        
        N_R = sinsertBeta(as.matrix(yStar,ncol=1), NX[,T_IDX,drop=FALSE], tau/KK, p,as.matrix(LAM,nrow=1));
        
        T_IDX <- (POS-1):(POS+1)+1;
        yStar = y - X[,-T_IDX,drop=FALSE]%*%BETAS[-T_IDX,,drop=FALSE];
        LAM = 0*abs(m.TREE[1,]+10)*sqrt(((m.TREE[1,]-a1+10))^2*(m.TREE[1,]-a2+10)^2)+lamT;
        
        C_R =  sdeleteBeta(yStar,X[,T_IDX,drop=FALSE],tau/KK,p,as.matrix(LAM[T_IDX-1],ncol=1));
        
        lnLike = N_R$LPROB; 
        ldLike = C_R$LPROB;
        
        lnpMove = log(PROB.PROPOSAL(L$m.TREE,m.TREE))
        ldpMove = log(PROB.PROPOSAL(m.TREE,L$m.TREE))
        
        
        lnPrior = LOG.PROB.TREE(ptree,L$m.TREE)
        ldPrior = LOG.PROB.TREE(ptree, m.TREE)
        
        test = lnLike + lnpMove + lnPrior - (ldLike + ldpMove + ldPrior); 
        
        
      }
      if (L$INSERT==FALSE){
        
        NCBX = shapesplineDelete(as.matrix(L$m.TREE[1,]-0.5,ncol=1),as.matrix(t,ncol=1)-0.5,0.1,as.integer(2),CBX,as.integer(L$pos))
        
        NX = cbind(ONES,100*(NCBX[,,3] - (a1+a2)*NCBX[,,2]+(a1*a2)*NCBX[,,1]));
        
        
        T_IDX <- (POS-1):(POS+1)+1;
        NBETAS = BETAS[-(POS+1),,drop=F];
        
        yStar = y - NX[,-T_IDX,drop=FALSE]%*%NBETAS[-T_IDX,,drop=FALSE];
        LAM = lamT
        N_R = sdeleteBeta(yStar,NX[,T_IDX,drop=FALSE],tau/KK,p,as.matrix(c(LAM,LAM,LAM),ncol=1));
        
        T_IDX <- (POS-1):(POS+2)+1;
        yStar = y - X[,-T_IDX,drop=FALSE]%*%BETAS[-T_IDX,,drop=FALSE];
        LAM = lamT
        
        C_R =  sinsertBeta(yStar, X[,T_IDX,drop=FALSE], tau/KK, p,as.matrix(c(LAM,LAM,LAM,LAM),ncol=1));
        
        lnLike = N_R$LPROB; 
        ldLike = C_R$LPROB; 
        
        lnpMove = log(PROB.PROPOSAL(L$m.TREE,m.TREE))
        ldpMove = log(PROB.PROPOSAL(m.TREE,L$m.TREE))
        
        
        lnPrior = LOG.PROB.TREE(ptree,L$m.TREE)
        ldPrior = LOG.PROB.TREE(ptree, m.TREE)
        
        test = lnLike + lnpMove + lnPrior - (ldLike + ldpMove + ldPrior); 
        
      }
      
      
      if ( test > 0 || runif(1)<exp(test))
      {
        m.TREE = L$m.TREE
        CBX   = NCBX; 
        X = NX;   
        
        # now sample the new beta
        if (L$INSERT){ 
          T_IDX <- (POS-1):(POS+2)+1;
          cmat = matrix(c( 0,0,0,0,   
                           1,0,0,0,   
                           0,2,0,0, 
                           0,0,3,0,
                           0,0,0,4,
                           1,2,0,0, 
                           1,0,3,0, 
                           1,0,0,4,
                           0,2,3,0,
                           0,2,0,4,
                           0,0,3,4,
                           0,2,3,4,
                           1,0,3,4,
                           1,2,0,4,
                           1,2,3,0,
                           1,2,3,4),ncol=4, nrow=16,byrow=T);
          tBETAS = c(0,0,0,0);
          tBETAS[cmat[N_R$S_ELM,]]= rtmvn(as.matrix(N_R$MEAN),as.matrix(N_R$VAR));
          NBETAS[T_IDX] = tBETAS; 
          #  print(c(i,N_R$S_ELM)); 
        }else{
          T_IDX <- (POS-1):(POS+1)+1;  
          cmat = matrix(c( 0,0,0,   
                           1,0,0,   
                           0,2,0,
                           0,0,3,
                           1,2,0,
                           1,0,3,
                           0,2,3,
                           1,2,3) ,ncol=3, nrow=8,byrow=T);
          tBETAS = c(0,0,0); 
          tBETAS[cmat[N_R$S_ELM,]]= rtmvn(as.matrix(N_R$MEAN),as.matrix(N_R$VAR));
          NBETAS[T_IDX] = tBETAS;
          #  print(c(i,N_R$S_ELM)); 
        }          
        BETAS = NBETAS;     
      }
      BP[[jj]]  = p
      BBETAS[[jj]] = BETAS; 
      BTAUS[[jj]] = tau; 
      BCBX[[jj]] = CBX; 
      BA[[jj]] = c(a1,a2); 
      BX[[jj]] = X; 
      BTREE[[jj]] = m.TREE; 
      BCLAM[[jj]] = lamT; 
    } else{
      #PERFORM A SWAP MOVE
      kk <- sample(1:nrow(ADJ),1)
      L1 = BCLAM[[ADJ[kk,1]]]
      L2 = BCLAM[[ADJ[kk,2]]]
      B1 = BBETAS[[ADJ[kk,1]]]
      B2 = BBETAS[[ADJ[kk,2]]]
      X1 = BX[[ADJ[kk,1]]]
      X2 = BX[[ADJ[kk,2]]]
      T1 = BTAUS[[ADJ[kk,1]]]
      T2 = BTAUS[[ADJ[kk,2]]]
      K1 = JK[[ADJ[kk,1]]]
      K2 = JK[[ADJ[kk,2]]]
      
      SS1 = (y-X1%*%B1); SS1 = t(SS1)%*%SS1; 
      SS2 = (y-X2%*%B2); SS2 = t(SS2)%*%SS2; 
      #note all of the prior probabilities are independent of K1 or K2 and thus
      #cancel in the ratio
      DEN  = -0.5*T1/K1*SS1  - 0.5*T2/K2*SS2 + length(y)/(K1*2)*log(T1)+ length(y)/(K2*2)*log(T2); 
      NUM  = -0.5*T1/K2*SS1  - 0.5*T2/K1*SS2 + length(y)/(K2*2)*log(T1)+ length(y)/(K1*2)*log(T2);#swap 
      
      test = NUM - DEN; 
      
      if (test >  0 || runif(1) < exp(test)){
        l = ADJ[kk,1]
        m = ADJ[kk,2]
        TB = BBETAS[[m]]; 
        TT = BTAUS[[m]];
        TC = BCBX[[m]]; 
        TA = BA[[m]];
        TX = BX[[m]];
        TTR= BTREE[[m]];
        TCL= BCLAM[[m]]; 
        
        BBETAS[[m]]= BBETAS[[l]];
        BTAUS[[m]] = BTAUS[[l]];
        BCBX[[m]]  = BCBX[[l]]; 
        BA[[m]]    = BA[[l]];
        BX[[m]]    = BX[[l]];
        BTREE[[m]] = BTREE[[l]];
        BCLAM[[m]] = BCLAM[[l]];
        
        BBETAS[[l]] = TB; 
        BTAUS[[l]] = TT;
        BCBX[[l]] = TC; 
        BA[[l]] = TA;
        BX[[l]] = TX;
        BTREE[[l]] = TTR;
        BCLAM[[l]] = TCL; 
        
      }
      
      
    }  
    tau_sample[i]     = BTAUS[[rr]]; 
    model_sample[[i]] = BTREE[[rr]]
    BETAS = BBETAS[[rr]]
    beta_sample[[i]] = BETAS
    h_p[i] = BP[[rr]]
    h_lambda[i] = BCLAM[[rr]]
    X = BX[[rr]]
    n_knots[i] = ncol(X)
    TEMP = X%*%BETAS; 
    
    t.mtree = BTREE[[rr]]; 
    beta_sample[[i]] = BETAS; 
    knot_sample[[i]] = t.mtree[1,]
    
    r = qcopy(x.Pred,TEMP,as.integer(length(y)),as.integer(i))
    TA = BA[[rr]]; 
    sa1[i] = TA[1];
    sa2[i] = TA[2];
    
    
    
  }
  
  
  lmEST = rowMeans(x.Pred[,nburn:nsamps],na.rm=T)
  
  return(list(lmEST=lmEST,sa1=sa1,sa2=sa2,beta_sample=beta_sample,model_sample=model_sample,x.Pred=x.Pred))  
}


lxfitDich  <- function(x,y,nsamps = 50000, nburn = 10000, a_list = c(1,1,1.1,1.2,1.3,1.7,2,2.6,3.5,5,9,12,24,30) ){
  #make sure the y variable is a matrix
  t_y = as.matrix(y,ncol=1)
  #center the x variable to be between 0 and 1
  t_x = as.matrix((x-min(x))/(max(x)-min(x)),ncol=1 )
  if (sum(dim(t_x)==dim(t_y))!=2){
    stop("X and Y vectors must be of the same dimension")
  }
  
  #error check the JK quantity making sure 1 is in the list and
  #it is sorted
  a_list = sort(a_list)
  if (a_list[1] != 1){
    stop("The minimum value in the anealling list must be 1")
  }
  
  if (length(a_list)==1){
    stop("Anealling list must have more than one element")
  }
  
  if (nsamps < nburn){
    stop("Number of samples must be greater than the number of burnin samples.")
  }
  
  if (nsamps < 0){
    stop("Number of samples must be greater than zero.")
  }
  
  if (nburn < 0){
    warning("Burn in samples less than zero. Setting burn in samples to zero.")
    nburn = 1; 
  }
  
  return (.LXsampleDich(t_x,t_y,a_list,nsamps,nburn))
  
}



#.LXsample This function 
.LXsampleDich <- function(t,xy,JK,nsamps,nburn){
  #create the initial tree list
  a1 = 0; a2 = 0; 
  m.TREE = matrix(-1,nrow=6,3)
  m.TREE[1,1:3] = as.matrix(c(0,0.5,1))
  m.TREE[2,1:3] = as.matrix(c(0,0.25,0))
  m.TREE[3,1:3] = as.matrix(c(0,0.75,0))
  m.TREE[4,1:3] = as.matrix(c(0,0.25,0)) 
  m.TREE[5,1:3] = as.matrix(c(0, 1, 0))
  m.TREE[6,1:3] = as.matrix(c(0, 1, 0))
  CBX = shapespline2(m.TREE[1,]-.5, t-0.5,2)
  
  ONES = matrix(1,nrow=length(t),ncol=1)
  X = cbind(ONES,100*(CBX[,,3] - (a1+a2)*CBX[,,2]+ (a1*a2)*CBX[,,1]));
  BETAS = matrix(0,ncol(X),ncol=1)
   
  ptree = 0.5
  tau = 0.25 
  x.Pred <- matrix(0,nrow=length(xy),ncol=nsamps) #mean prediction of the data points
  lamT = 1
  intLAM = 0.001
  p = 0.01
  
  BBETAS = vector("list",length(JK));
  BTAUS  = vector("list",length(JK)); 
  BA     = vector("list",length(JK)); 
  BX     = vector("list",length(JK));
  BCBX   = vector("list",length(JK)); 
  BTREE  = vector("list",length(JK)); 
  BCLAM  = vector("list",length(JK));
  BP     = vector("list",length(JK));
  
  for (i in 1:length(JK)){
    BBETAS[[i]] = BETAS; 
    BTAUS[[i]] = tau; 
    BCBX[[i]] = CBX; 
    if (i ==2) {a1 = -0.5; a2 = -0.5}
    if (i ==3) {a1 = 0.5; a2 = 0.5}
    BA[[i]] = c(a1,a2); 
    
    BX[[i]] = X; 
    BTREE[[i]] = m.TREE; 
    BCLAM[[i]] = lamT; 
    BP[[i]]    = p; 
    
  }
  
  MT  <- 1:(length(JK)-1)
  ADJ <- matrix(0,nrow=length(JK)-1,ncol=2)
  ADJ[,1] =  MT; ADJ[,2] =  MT+1; 
  ADJ = rbind(ADJ,c(1,3))
  ADJ = rbind(ADJ,c(1,4))
  ADJ = rbind(ADJ,c(2,4))
  ADJ = rbind(ADJ,c(2,5))
  
  
  rr = 1
  
  
  sa1      =  rep(0,nsamps)
  sa2      =  rep(0,nsamps)
  n_knots  =  rep(0,nsamps)
  h_lambda = rep(0,nsamps)
  h_p      = rep(0,nsamps)
  model_sample     = vector("list",nsamps);
  beta_sample      = vector("list",nsamps);
  knot_sample      = vector("list",nsamps); 
  tau_sample       = 1:nsamps
  
  
  pb <- txtProgressBar(min=0, max=nsamps, style=3)
  
  
  for (i in 1:nsamps){
    
    setTxtProgressBar(pb,i)    
  
    if (runif(1) < 0.5){ # simple update 
 
    
      jj <- sample(1:length(JK),1)
      KK = JK[jj];
      BETAS = BBETAS[[jj]]; 
      tau = BTAUS[[jj]];
      CBX = BCBX[[jj]];  
      A = BA[[jj]]; a1 = A[1]; a2 = A[2]; 
      X = BX[[jj]]; 
      m.TREE = BTREE[[jj]]; 
      lamT = BCLAM[[jj]]; 
     
      #################################################
      #################################################
      # Albert and Chib Probit data Augmentation
      #################################################
      X = cbind(ONES,100*(CBX[,,3] - (a1+a2)*CBX[,,2]+(a1*a2)*CBX[,,1]));
      ac.mean <- X%*%BETAS
      y <- TnormV(ac.mean*(1-2*(xy == 0)),rep(1,length(xy)))*(1-2*(xy == 0))  #flip the mean if we observed a 0
      							 #this way we can simulat all positive truncated
      							 #values and flip them back when we need to									

      
      # invisible(readline(prompt="Press [enter] to continue"))
      #################################
      ##HYPER PRIOR LAMBDA
      betaT = BETAS[2:length(BETAS)]
      betaT = betaT[betaT > 0]
      
      ta = length(betaT) + .2
      tb = sum(betaT) + 2
      tg_ub = pgamma(1e-5,ta,tb)
      lamT = qgamma(runif(1,tg_ub,1),ta,tb)
      
      #other hyper prior on probability of a zero
      betaT = BETAS[2:length(BETAS)]
      V = betaT == 0
      t.n = length(betaT)
      ta = 2 + sum(V)
      tb = 18 + t.n - sum(V)
      p = rbeta(1,ta,tb)
      
      #################################
      X = cbind(ONES,100*(CBX[,,3] - (a1+a2)*CBX[,,2]+(a1*a2)*CBX[,,1]));
      ###############################################################
      LAM = lamT
      BETAS = sampleBetas(y, X, BETAS,LAM, intLAM, p,tau/KK)
   
      ################################################################
      
      tB <- 0.5*t(X%*%BETAS-y)%*%(X%*%BETAS-y)/KK+1
      tA <- 0.5*length(y)/KK+ 1
      tau <- as.vector((1/tB)*rgamma(1,tA,1))
      
      ystar =  y - cbind(ONES,100*(CBX[,,3] - (a2)*CBX[,,2]))%*%BETAS;
      w = 100*(-CBX[,,2]+(a2)*CBX[,,1])%*%BETAS[-1]
      a1  =  sampA(ystar,w,tau/KK,0,1,-0.5,0.5)
      
      
      X = cbind(ONES,100*(CBX[,,3] - (a1+a2)*CBX[,,2]+(a1*a2)*CBX[,,1]));
      
      ystar =  y - cbind(ONES,100*(CBX[,,3] - (a1)*CBX[,,2]))%*%BETAS; 
      w = 100*(-CBX[,,2]+(a1)*CBX[,,1])%*%BETAS[-1]
      a2  =  sampA(ystar,w,tau/KK,0,1,-0.5,0.5)
      
      X = cbind(ONES,100*(CBX[,,3] - (a1+a2)*CBX[,,2]+(a1*a2)*CBX[,,1]));
      
      L = PROP.N.TREE(m.TREE);
      
      POS = L$pos;   
      
      
      if (L$INSERT==TRUE){
        
        NCBX = shapesplineInsert(as.matrix(L$m.TREE[1,]-0.5,nrow=1), as.matrix(t,ncol=1)-0.5, 0.1,2,CBX,L$pos)
       
        NX = cbind(ONES,100*(NCBX[,,3] - (a1+a2)*NCBX[,,2]+(a1*a2)*NCBX[,,1]));
        
        NBETAS = matrix(0,nrow=ncol(NX),ncol=1); 
        
        T_IDX <- (POS-1):(POS+2)+1;
        NBETAS[-T_IDX] = BETAS[-((POS-1):(POS+1)+1)]; 
        yStar = y - NX[,-T_IDX]%*%NBETAS[-T_IDX];
        LAM = 0*abs(L$m.TREE[1,]+10)*sqrt(((L$m.TREE[1,]-a1+10))^2*(L$m.TREE[1,]-a2+10)^2)+lamT; 
        
        N_R = sinsertBeta(as.matrix(yStar,ncol=1), NX[,T_IDX,drop=FALSE], tau/KK, p,as.matrix(LAM,nrow=1));
        
        T_IDX <- (POS-1):(POS+1)+1;
        yStar = y - X[,-T_IDX,drop=FALSE]%*%BETAS[-T_IDX,,drop=FALSE];
        LAM = 0*abs(m.TREE[1,]+10)*sqrt(((m.TREE[1,]-a1+10))^2*(m.TREE[1,]-a2+10)^2)+lamT;
        
        C_R =  sdeleteBeta(yStar,X[,T_IDX,drop=FALSE],tau/KK,p,as.matrix(LAM[T_IDX-1],ncol=1));
        
        lnLike = N_R$LPROB; 
        ldLike = C_R$LPROB;
        
        lnpMove = log(PROB.PROPOSAL(L$m.TREE,m.TREE))
        ldpMove = log(PROB.PROPOSAL(m.TREE,L$m.TREE))
        
        
        lnPrior = LOG.PROB.TREE(ptree,L$m.TREE)
        ldPrior = LOG.PROB.TREE(ptree, m.TREE)
        
        test = lnLike + lnpMove + lnPrior - (ldLike + ldpMove + ldPrior); 
        
        
      }
      if (L$INSERT==FALSE){
        
        NCBX = shapesplineDelete(as.matrix(L$m.TREE[1,]-0.5,ncol=1),as.matrix(t,ncol=1)-0.5,0.1,as.integer(2),CBX,as.integer(L$pos))
        
        NX = cbind(ONES,100*(NCBX[,,3] - (a1+a2)*NCBX[,,2]+(a1*a2)*NCBX[,,1]));
        
        
        T_IDX <- (POS-1):(POS+1)+1;
        NBETAS = BETAS[-(POS+1),,drop=F];
        
        yStar = y - NX[,-T_IDX,drop=FALSE]%*%NBETAS[-T_IDX,,drop=FALSE];
        LAM = lamT
        N_R = sdeleteBeta(yStar,NX[,T_IDX,drop=FALSE],tau/KK,p,as.matrix(c(LAM,LAM,LAM),ncol=1));
        
        T_IDX <- (POS-1):(POS+2)+1;
        yStar = y - X[,-T_IDX,drop=FALSE]%*%BETAS[-T_IDX,,drop=FALSE];
        LAM = lamT
        
        C_R =  sinsertBeta(yStar, X[,T_IDX,drop=FALSE], tau/KK, p,as.matrix(c(LAM,LAM,LAM,LAM),ncol=1));
        
        lnLike = N_R$LPROB; 
        ldLike = C_R$LPROB; 
        
        lnpMove = log(PROB.PROPOSAL(L$m.TREE,m.TREE))
        ldpMove = log(PROB.PROPOSAL(m.TREE,L$m.TREE))
        
        
        lnPrior = LOG.PROB.TREE(ptree,L$m.TREE)
        ldPrior = LOG.PROB.TREE(ptree, m.TREE)
        
        test = lnLike + lnpMove + lnPrior - (ldLike + ldpMove + ldPrior); 
        
      }
      
      
      if ( test > 0 || runif(1)<exp(test))
      {
        m.TREE = L$m.TREE
        CBX   = NCBX; 
        X = NX;   
        
        # now sample the new beta
        if (L$INSERT){ 
          T_IDX <- (POS-1):(POS+2)+1;
          cmat = matrix(c( 0,0,0,0,   
                           1,0,0,0,   
                           0,2,0,0, 
                           0,0,3,0,
                           0,0,0,4,
                           1,2,0,0, 
                           1,0,3,0, 
                           1,0,0,4,
                           0,2,3,0,
                           0,2,0,4,
                           0,0,3,4,
                           0,2,3,4,
                           1,0,3,4,
                           1,2,0,4,
                           1,2,3,0,
                           1,2,3,4),ncol=4, nrow=16,byrow=T);
          tBETAS = c(0,0,0,0);
          tBETAS[cmat[N_R$S_ELM,]]= rtmvn(as.matrix(N_R$MEAN),as.matrix(N_R$VAR));
          NBETAS[T_IDX] = tBETAS; 
          #  print(c(i,N_R$S_ELM)); 
        }else{
          T_IDX <- (POS-1):(POS+1)+1;  
          cmat = matrix(c( 0,0,0,   
                           1,0,0,   
                           0,2,0,
                           0,0,3,
                           1,2,0,
                           1,0,3,
                           0,2,3,
                           1,2,3) ,ncol=3, nrow=8,byrow=T);
          tBETAS = c(0,0,0); 
          tBETAS[cmat[N_R$S_ELM,]]= rtmvn(as.matrix(N_R$MEAN),as.matrix(N_R$VAR));
          NBETAS[T_IDX] = tBETAS;
          #  print(c(i,N_R$S_ELM)); 
        }          
        BETAS = NBETAS;     
      }
      BP[[jj]]  = p
      BBETAS[[jj]] = BETAS; 
      BTAUS[[jj]] = tau; 
      BCBX[[jj]] = CBX; 
      BA[[jj]] = c(a1,a2); 
      BX[[jj]] = X; 
      BTREE[[jj]] = m.TREE; 
      BCLAM[[jj]] = lamT; 
    } else{
      #PERFORM A SWAP MOVE
      kk <- sample(1:nrow(ADJ),1)
      L1 = BCLAM[[ADJ[kk,1]]]
      L2 = BCLAM[[ADJ[kk,2]]]
      B1 = BBETAS[[ADJ[kk,1]]]
      B2 = BBETAS[[ADJ[kk,2]]]
      X1 = BX[[ADJ[kk,1]]]
      X2 = BX[[ADJ[kk,2]]]
      T1 = BTAUS[[ADJ[kk,1]]]
      T2 = BTAUS[[ADJ[kk,2]]]
      K1 = JK[[ADJ[kk,1]]]
      K2 = JK[[ADJ[kk,2]]]
      
      ###########################################################
      #ALBERT and CHIB AUGMENTATION step
      #
      ###########################################################
      ac.mean <- X1%*%B1
      y1 <- TnormV(ac.mean*(1-2*(xy == 0)),rep(1,length(xy)))*(1-2*(xy == 0))  #flip the mean if we observed a 0
      							 #this way we can simulat all positive truncated
      							 #values and flip them back when we need to									
      
      ac.mean <- X2%*%B2
      y2 <- TnormV(ac.mean*(1-2*(xy == 0)),rep(1,length(xy)))*(1-2*(xy == 0)) #flip the mean if we observed a 0
      								 #this way we can simulat all positive truncated
      								 #values and flip them back when we need to									
 
      ##################################################################
      SS1 = (y1-X1%*%B1); SS1 = t(SS1)%*%SS1; 
      SS2 = (y2-X2%*%B2); SS2 = t(SS2)%*%SS2; 
      #note all of the prior probabilities are independent of K1 or K2 and thus
      #cancel in the ratio
      DEN  = -0.5*T1/K1*SS1  - 0.5*T2/K2*SS2 + length(xy)/(K1*2)*log(T1)+ length(xy)/(K2*2)*log(T2); 
      NUM  = -0.5*T1/K2*SS1  - 0.5*T2/K1*SS2 + length(xy)/(K2*2)*log(T1)+ length(xy)/(K1*2)*log(T2);#swap 
      
      test = NUM - DEN; 
      
      if (test >  0 || runif(1) < exp(test)){
        l = ADJ[kk,1]
        m = ADJ[kk,2]
        TB = BBETAS[[m]]; 
        TT = BTAUS[[m]];
        TC = BCBX[[m]]; 
        TA = BA[[m]];
        TX = BX[[m]];
        TTR= BTREE[[m]];
        TCL= BCLAM[[m]]; 
        
        BBETAS[[m]]= BBETAS[[l]];
        BTAUS[[m]] = BTAUS[[l]];
        BCBX[[m]]  = BCBX[[l]]; 
        BA[[m]]    = BA[[l]];
        BX[[m]]    = BX[[l]];
        BTREE[[m]] = BTREE[[l]];
        BCLAM[[m]] = BCLAM[[l]];
        
        BBETAS[[l]] = TB; 
        BTAUS[[l]] = TT;
        BCBX[[l]] = TC; 
        BA[[l]] = TA;
        BX[[l]] = TX;
        BTREE[[l]] = TTR;
        BCLAM[[l]] = TCL; 
        
      }
      
      
    }  
    tau_sample[i]     = BTAUS[[rr]]; 
    model_sample[[i]] = BTREE[[rr]]
    BETAS = BBETAS[[rr]]
    beta_sample[[i]] = BETAS
    h_p[i] = BP[[rr]]
    h_lambda[i] = BCLAM[[rr]]
    X = BX[[rr]]
    n_knots[i] = ncol(X)
    TEMP = X%*%BETAS; 
    
    t.mtree = BTREE[[rr]]; 
    beta_sample[[i]] = BETAS; 
    knot_sample[[i]] = t.mtree[1,]
    
    r = qcopy(x.Pred,TEMP,as.integer(length(xy)),as.integer(i))
    TA = BA[[rr]]; 
    sa1[i] = TA[1];
    sa2[i] = TA[2];
    
    
    
  }
  
  
  lmEST = rowMeans(x.Pred[,nburn:nsamps],na.rm=T)
  
  return(list(lmEST=lmEST,sa1=sa1,sa2=sa2,beta_sample=beta_sample,model_sample=model_sample,x.Pred=x.Pred))  
}

#########################################
# Calculate Posterior model probabilities 
# nburn - is the number of burnin samples to ignore
# nsamps -  is the maximum sample - needs to be less than 
# the number of samples used in the fit
#########################################

calculateProbs <- function(fit,nburn,nsamps){
  sa1 = fit$sa1
  sa2 = fit$sa2
  if (nsamps > length(sa2)){
    stop("Number of samples (nsamps) must be less than or equal to the total number of 
          samples")
    
  }
  SHAPES =  rep(0,nsamps)
  for (ii in nburn:nsamps) {
    a.s <- c(sa1[ii],sa2[ii])
    betas = fit$beta_sample[[ii]]
    knots = c(fit$model_sample[[ii]][1,])
    betas = betas[2:length(betas)]
    flat  = rep(0,length(knots)-1)
    n = length(betas)
    #decide 'base shape'
    shape = 1; #'n' shape
    shape =  ((sa1[ii] > -.5)*(sa2[ii]> -.5)*(sa1[ii] < 0.5)*(sa2[ii]< 0.5))*4 + shape; #s shaped
    shape = ((sa1[ii] <= -0.5)*(sa2[ii] >= 0.5)+(sa1[ii] >= 0.5)*(sa2[ii] <= -0.5))*2 + shape# Monotone Decreasing
    shape = ((sa1[ii] <= -0.5)*(sa2[ii] <= -0.5)+(sa1[ii] >= 0.5)*(sa2[ii] >= 0.5))*3 + shape # Monotone Increasing
    shape = ((sa1[ii] > -.5)*(sa1[ii]< 0.5)*(sa2[ii] <= -0.5)+(sa2[ii] > -.5)*(sa2[ii]< 0.5)*(sa1[ii] <= -0.5))*1 +shape#J-shaped
    
    flat[1] = (betas[1] == 0) + (betas[2] ==  0)
    flat[2] = (betas[2] == 0)
    flat[length(flat)] = (betas[n] == 0) + (betas[n-1] == 0)
    flat[length(flat)-1] = (betas[n-1] == 0)
    if (length(flat) > 2){
      for (jj in 3:(n-2)){
        flat[jj-2] = flat[jj-2] + (betas[jj] == 0)
        flat[jj-1] = flat[jj-1] + (betas[jj] == 0)
        flat[jj] = flat[jj] + (betas[jj] == 0)
      }
    }
    intLoc1 = sum(sa1[ii] > knots-0.5)
    intLoc2 = sum(sa2[ii] > knots-0.5)
    temp = sort(c(intLoc1,intLoc2)+1)
    temp = c(1,temp,length(flat)+2)
    flat = c(1,( flat ==  3),1) #if there are three flat betas the region is flat
    
    #different 'checks' depending on the starting shape of the curve  
    if (shape == 5){
      left = (prod(flat[temp[1]:temp[2]]) ==1)
      center = (prod(flat[temp[2]:temp[3]]) ==1)
      right = (prod(flat[temp[3]:temp[4]]) ==1)
      if  (center){
        shape = 4
      }else if (left){
        if  (right){shape = 2}else{shape=3}
      }else  if(right){
        shape = 1    
      }
    }else if(shape == 1 ){
      temp  = unique(temp)
      left = (prod(flat[temp[1]:temp[2]]) ==1)
      right = (prod(flat[temp[2]:temp[3]]) ==1)
      if (left){
        shape  = 3
      }else  if(right){
        shape  = 4    
      }
    }else if(shape == 2){
      temp  = unique(temp)
      left = (prod(flat[temp[1]:temp[2]]) ==1)
      right = (prod(flat[temp[2]:temp[3]]) ==1)
      if (left){
        shape=3
      }else  if(right){
        shape = 2    
      }
    }
    SHAPES[ii] = shape  
  }
  
  pr1 = mean(SHAPES == 5) #S SHAPED
  pr2 = mean(SHAPES == 3) #Monotone Decreasing
  pr3 = mean(SHAPES == 4) #Monotone Increasing
  pr4 = mean(SHAPES == 2) #J shaped
  pr5 = 1-pr1-pr2-pr3-pr4 # n- shaped
  return (list(pr1=pr1,pr2=pr2,pr3=pr3,pr4=pr4,pr5=pr5,SHAPES = SHAPES))
}

