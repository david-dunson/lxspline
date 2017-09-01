#include <RcppArmadillo.h>
#include <iostream>

using namespace std; 
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

/*     N      INTEGER, the number of variables.
*     LOWER  REAL, array of lower integration limits.
*     UPPER  REAL, array of upper integration limits.
*     INFIN  INTEGER, array of integration limits flags:
*            if INFIN(I) < 0, Ith limits are (-infinity, infinity);
*            if INFIN(I) = 0, Ith limits are (-infinity, UPPER(I)];
*            if INFIN(I) = 1, Ith limits are [LOWER(I), infinity);
*            if INFIN(I) = 2, Ith limits are [LOWER(I), UPPER(I)].
*     CORREL REAL, array of correlation coefficients; the correlation
*            coefficient in row I column J of the correlation matrix
*            should be stored in CORREL( J + ((I-2)*(I-1))/2 ), for J < I.
*     MAXPTS INTEGER, maximum number of function values allowed. This 
*            parameter can be used to limit the time taken. A 
*            sensible strategy is to start with MAXPTS = 1000*N, and then
*            increase MAXPTS if ERROR is too large.
*     ABSEPS REAL absolute error tolerance.
*     RELEPS REAL relative error tolerance.
*     ERROR  REAL estimated absolute error, with 99% confidence level.
*     VALUE  REAL estimated value for the integral
*     INFORM INTEGER, termination status parameter:
*            if INFORM = 0, normal completion with ERROR < EPS;
*            if INFORM = 1, completion with ERROR > EPS and MAXPTS 
*                           function vaules used; increase MAXPTS to 
*                           decrease ERROR;
*            if INFORM = 2, N > 20 or N < 1.
*/
extern "C" {
  double  bvnd_(double * DH, double *DK, double *R);
  double  tvtl_(int * NU, double *HK, double *R, double *EPSI);
  double  sadmvn_(int * N, double * LOWER, double* UPPER, 
                  int * INFIN, double* CORREL, int *MAXPTS, 
                  double * ABSEPS,double * RELEPS, double * E, double * VALUE, int *INFORM );
}


// [[Rcpp::export]]
double SADMVN(arma::mat M, arma::mat C){
  int    N = 4; 
  int    MAXPTS = 2000; 
  double RELEPS = 0; 
  double ABSEPS = 1e-6; 
  double E_E     = 0.0; 
  double returnV = 0.0; 
  double UPPER[4]   = {std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity()}; //POPULATED NEVER USED
  int INFIN[4]   = {1,1,1,1}; 
  int INFORM; 
  
  
  // Set up the correlation matrix
  int gk = 0; 
  double B[6]; 
  for (int i = 1; i < 4; i++){
    M(i-1,0) = M(i-1,0)/sqrt(C(i-1,i-1)); 
    for (int j = 0; j < i; j++){
      B[gk] = C(i,j)/sqrt(C(j,j)*C(i,i));
      gk++; 
    }
  }
  M(3,0) = M(3,0)/sqrt(C(3,3)); 
  
  M = -M; 
  sadmvn_(&N,M.memptr(),UPPER,INFIN,B,&MAXPTS,&ABSEPS,&RELEPS,&E_E,&returnV,&INFORM); 
  
  return returnV; 
  
  
}


// [[Rcpp::export]]
double BivNProb(NumericVector mean, NumericVector cv){
  // Function BivNProb:
  //Wrapper function that computes the probability  that X > 0
  // by calling the fortran function bvnd_
  // assumes everythign is great no error checking
  double rv;
  double a =  -1.0*mean[0]/sqrt(cv[0]);
  double b =  -1.0*mean[1]/sqrt(cv[3]);
  double rho = cv[1]/sqrt(cv[0]*cv[3]);
  rv = bvnd_( &a,  &b,  &rho);
  return rv;
}

// [[Rcpp::export]]
double TriNProb(NumericVector mean, NumericVector cv){
  // Function BivNProb:
  //Wrapper function that computes the probability  that X > 0
  // by calling the fortran function bvnd_
  // assumes everythign is great no error checking
  double m[3];
  double R[3];
  double RV;  
  int df = 0;
  double EPS = 2e-16;
  
  //  compute the means
  m[0] = double(mean[0])/sqrt(cv[0]);
  m[1] = double(mean[1])/sqrt(cv[4]);
  m[2] = double(mean[2])/sqrt(cv[8]);
  
  // compute the correlations
  R[0]  = cv[1]/sqrt(cv[0]*cv[4]);
  R[1]  = cv[2]/sqrt(cv[0]*cv[8]);
  R[2]  = cv[5]/sqrt(cv[4]*cv[8]);
  
  RV = tvtl_( &df,  m,  R, &EPS);
  return RV;
}


// Generate a truncated normal truncated at Zero
// no error checking ... built for speed

int genTruncNormZ(double* mean, double* sd, double* rV){
  
  if (*mean < 0){
    // do everything on the log scale  - more numerically stable
    // when mean/sd < -8 
    double ucutoff = -1.0*R::pnorm(0.0,-1.0*(*mean),*sd,true, true);
    GetRNGstate();
    double cv = R::rexp(1) + ucutoff;
    PutRNGstate();
    *rV = -1.0*R::qnorm(-1.0*cv,-1.0*(*mean),*sd,true,true);
    if (*rV < 0){*rV = 1e-5;}//very rare numerical problem
  }else{
    // no real problems when 0 < mean
    double lcutoff = R::pnorm(0.0,*mean,*sd,true, false);
    GetRNGstate();
    double tunif = R::runif(lcutoff,1.0);
    PutRNGstate();
    *rV = R::qnorm(tunif,*mean,*sd,true,false);
  }
  
  
  
  
  return  0;
}
// [[Rcpp::export]]
NumericVector Tnorm(NumericVector mean, NumericVector sd){
  double rval; 
  double m = mean[0]; double sd1 = sd[0]; 
  genTruncNormZ( &m, &sd1, &rval); 
  return wrap(rval);
  
}

// [[Rcpp::export]]
NumericVector TnormV(NumericVector mean, NumericVector sd){
  NumericVector returnY = mean; 
  double rval;	
  double m; double sd1; 
  for (int i = 0; i < returnY.length() ; i++){
	  m = mean[i];
	  sd1 = sd[i];
	  genTruncNormZ( &m, &sd1, &rval);
	  returnY[i] = rval; 
  } 
  return wrap(returnY);
  
}
// [[Rcpp::export]]
NumericVector sampleBetas(NumericVector ttY,  NumericVector ttX, NumericVector tbetas,
                          NumericVector LAM,  NumericVector intLAM, NumericVector p,
                          NumericVector tau){
  
  
  arma::mat Y       = Rcpp::as<arma::mat>(ttY); 
  arma::mat X       = Rcpp::as<arma::mat>(ttX); 
  arma::mat betas   = Rcpp::as<arma::mat>(tbetas); 
  
  
  arma::mat tX = X.col(0); // first column
  Y = Y -  X.submat(0,1,X.n_rows-1,X.n_cols-1) * betas.submat(1,0,X.n_cols-1,0);
  
  arma::mat tV = 1/((tau[0])*(tX.t()*tX)+intLAM[0]);
  arma::mat tM = ((tau[0])*(tX.t()*Y)+intLAM[0]*0.0)*tV;
  GetRNGstate();
  betas(0,0) = R::rnorm(tM(0,0),sqrt(tV(0,0)));
  PutRNGstate();
  Y = Y - tX*betas(0,0);
  double PL = LAM[0];
  //  sample  each  constrained beta one at a time
  for (int i = 1; i < X.n_cols; i++){
    
    tX = X.submat(0,i,X.n_rows-1,i);
    Y  = Y + tX*betas(i,0); // the  previous iteration this value was r
    // removed from  the  y vector  'add it back'
    
    tV = ((tau[0])*(tX.t()*tX)).i();
    tM = ((tau[0])*(tX.t()*Y)-PL)*tV;
    
    
    double A = log(p[0])+R::dnorm(0.0,tM(0,0),sqrt(tV(0,0)),true);
    double B = log((1.0-p[0])*(PL))+R::pnorm(0.0,tM(0,0),sqrt(tV(0,0)),false,true);
    
    if (A > B){
      A = A-A;
      B = B-A;
    }else{
      A = A-B;
      B = B-B;
    }
    
    double PZERO = exp(A)/(exp(A)+exp(B));
    //  Rprintf("A: %f B:%f m: %f  tV %f\n",A,B,tM(0,0),tV(0,0));
    
    if (R::runif(0.0,1.0) < PZERO){
      betas(i,0)  = 0.0;
    }else{
      double rV = 0.0;
      double mean = tM(0,0); double sd = sqrt(tV(0,0));
      
      genTruncNormZ( &mean, &sd, &rV);
      if (rV < 0.0){rV = 1e-16;}
      betas(i,0) = rV;
      
    }
    // Rprintf("A: %f B:%f m: %f  tV %f betas:%f \n",A,B,tM(0,0),tV(0,0),betas(i,0));
    Y  = Y - tX*betas(i,0); // add the new *regressor* to the residual
    
  }
  
  
  return wrap(betas); 
}


// [[Rcpp::export]]
List sinsertBeta(NumericVector tY, NumericVector Xidx, NumericVector ctau,
                 NumericVector tp, NumericVector lam){
  
  // initialize the proper values
  double p = tp[0];
  arma::mat Y = Rcpp::as<arma::mat>(tY); 
  arma::mat X = Rcpp::as<arma::mat>(Xidx);
  double tau = ctau[0];
  
  arma::mat tLAM = Rcpp::as<arma::mat>(lam); 
  arma::mat LAM(tLAM.n_rows,tLAM.n_rows);  LAM.diag()  = tLAM;
  arma::mat INTC(16,2);
  arma::mat PMAT(16,1); 
  
  INTC.zeros();
  arma::mat B = (tau)*(X.t()*X);      arma::mat U = (tau)*(X.t()*Y)-lam[0];
  
  List uniMeans(4); List uniVars(4);
  List bivMeans(6); List bivVars(6); 
  List TriMeans(4); List TriVars(4); 
  
  List returnV(4); 
  StringVector names(4);
  
  names(0) = "VAR";    names(1) = "MEAN"; 
  names(2) = "LPROB";  names(3) = "S_ELM"; 
  
  returnV.attr("names") = names;
  
  double cnum = arma::cond(B); 
  
  if ( cnum > 1e5){
    returnV[0] = wrap(0);
    returnV[1] = wrap(0); 
    returnV[2] = wrap(log(0)); 
    returnV[3] = wrap(1); 
    return returnV; 
  }
  
  PMAT(0,0) = (p)*(p)*(p)*(p); 
  
  arma::mat tV;            arma::mat tM;
  arma::mat TEMP; double pr = 0.0;
  // The  univariate values
  for (int i  = 0; i < 4;i++ ){
    PMAT(i+1,0) = p*p*p*(1-p);
    tV = (B(i,i));
    tV = tV.i();
    tM = tV*(U(i,0));
    TEMP = -.5*tM.t()*solve(tV,tM);
    INTC(i+1,0) = TEMP(0,0);
    //Total positive area for the univariate normal distribution
    pr = 1.0-R::pnorm(0,tM(0,0),sqrt(tV(0,0)),true,false);
    
    if (pr <= 0 || isnan(pr)){pr = 0;}// rare numerical error when probability is essentially zero
    TEMP = log(pr*sqrt(tV)*lam[0])+0.5*log(2.0*M_PI);
    INTC(i+1,1) = TEMP(0,0);
    uniMeans[i] = wrap(tM(0,0)); uniVars[i] = wrap(tV(0,0));
  }
  //bivariate normal case
  arma::mat tU; 
  arma::mat tB; 
  int counter = 0; 
  for (int i  = 0; i < 3;i++ ){
    for (int j = i+1; j < 4; j++){
      PMAT(counter + 5, 0) = p*p*(1-p)*(1-p);
      tU = join_cols(U.row(i),U.row(j)); 
      tB = join_rows(B.col(i),B.col(j));
      for ( int k = 3; k >= 0; k--){ 
        if (k != i && k != j)
          tB.shed_row(k);
      }
      tV = (tB);     tV = tV.i();
      tM = tV*(tU);
      
      bivVars[counter] = wrap(tV);
      bivMeans[counter] = wrap(tM);
      TEMP = -.5*tM.t()*solve(tV,tM);
      INTC(counter+5,0) = TEMP(0,0);
      pr = BivNProb(wrap(tM),wrap(tV));  if (pr <= 0 || isnan(pr)){pr = 0;}// numerical instability for small probabilities
      TEMP = log(pr*(sqrt(abs(det(tV)))*lam[0]*lam[0]))+log(2.0*M_PI);
      
      INTC(counter+5,1) = TEMP(0,0);
      counter++; 
    }
  }
  
  //trivariate case 
  for (int i =0; i<4; i++){
    PMAT(i+11,0) = p*(1-p)*(1-p)*(1-p);
    tU = U; tU.shed_row(i); 
    tB = B; tB.shed_col(i); tB.shed_row(i);
    tV = (tB);    tV = tV.i();
    tM = tV*(tU);
    TriVars[i] = wrap(tV);
    TriMeans[i] = wrap(tM);
    TEMP = -.5*tM.t()*solve(tV,tM);
    INTC(i+11,0) = TEMP(0,0);
    pr = TriNProb(wrap(tM), wrap(tV));   if (pr <= 0 || isnan(pr)){pr = 0;}// numerical instability for small probabilities
    
    TEMP = log(pr*(sqrt(abs(det(tV)))*lam[0]*lam[0]*lam[0]))+1.5*log(2.0*M_PI);
    INTC(i+11,1) = TEMP(0,0);
  }
  
  // 4-variate normal
  PMAT(15,0) = (1-p)*(1-p)*(1-p)*(1-p);
  tU = U;    tB = B; 
  tV = (tB);    tV = tV.i();
  tM = tV*(tU);
  TEMP = -.5*tM.t()*solve(tV,tM);
  INTC(15,0) = TEMP(0,0);
  pr = SADMVN(tM, tV);   if (pr <= 2e-16 || isnan(pr)){pr = 0;} // rare numerical error caused by low probabilities and a large 
  // unstable determinant 
  
  TEMP = log(pr*(sqrt(abs(det(tV)))*lam[0]*lam[0]*lam[0]*lam[0]))+2*log(2.0*M_PI);
  if (TEMP(0,0) > 10){TEMP(0,0)=log(0);} // another rare numerical problem
  INTC(15,1) = TEMP(0,0);
  
  arma::mat TT = INTC.col(1)-INTC.col(0);
  
  //cout << TT << endl; 
  arma::mat MT = max(TT);
  arma::mat SMULT = PMAT%exp(TT-repmat(MT,16,1));
  arma::mat log_prob = -0.5*(tau)*(Y.t()*Y) + log(sum(SMULT)) + MT;
  arma::mat POST_PROBS = exp(log(SMULT)-log(repmat(sum(SMULT),16,1))); //posterior sampling prob
  GetRNGstate();
  // get the cumulative  density
  double rvalue = R::runif(0,1); int rval = -1;
  
  for (int i = 1; i < 16; i++){  POST_PROBS(i,0) +=  POST_PROBS(i-1,0); }
  for (int i = 0; i < 16 && rval ==-1; i++){if (rvalue <= POST_PROBS(i,0)){ rval = i;}}
  
  
  rval = rval+1;
  PutRNGstate();
  NumericVector rVar, rMeans;
  
  switch(rval){
  case 16:
    rVar = wrap(tV); rMeans = wrap(tM);
    break;
  case 1:
    rVar  = 0.0;    rMeans = 0.0;
    break;
  case 2: case 3: case 4: case 5: // univariate chosen
    rVar   = uniVars[rval-2];
    rMeans = uniMeans[rval-2];
    break;
  case 12: case 13: case 14: case 15:
    rVar   = TriVars[rval-12]; 
    rMeans = TriMeans[rval-12];
    break; 
  default:
    rVar = bivVars[rval-6];
  rMeans = bivMeans[rval-6];
  break;
  }
  
 
  returnV[0] = wrap(rVar);
  returnV[1] = wrap(rMeans); 
  returnV[2] = wrap(log_prob(0,0)); 
  returnV[3] = wrap(rval);
  
  
 
  return returnV; 
}
// [[Rcpp::export]]
List sdeleteBeta(NumericVector tY, NumericVector Xidx, NumericVector ctau,
                 NumericVector tp, NumericVector lam){
  
  // initialize the proper values
  arma::mat Y = Rcpp::as<arma::mat>(tY); 
  arma::mat X = Rcpp::as<arma::mat>(Xidx);
  double tau = ctau[0];
  
  arma::mat tLAM = Rcpp::as<arma::mat>(lam); 
  arma::mat LAM(tLAM.n_rows,tLAM.n_rows);  LAM.diag()  = tLAM;
  arma::mat INTC(8,2);
  
  INTC.zeros();
  arma::mat B = (tau)*(X.t()*X);      arma::mat U = (tau)*(X.t()*Y);
  double p = tp[0];
  
  arma::mat tV;            arma::mat tM;
  arma::mat uniMeans(3,1); arma::mat uniVars(3,1);
  arma::mat bivMeans(2,3); arma::mat bivVars(2,6);
  
  arma::mat TEMP; double pr = 0.0;
  arma::mat PMAT(8,1); PMAT(0,0) = (p)*(p)*(p);
  
  // start with the univariate values
  for (int i  = 1; i <4;i++ ){
    PMAT(i,0) = p*p*(1-p);
    tV = (B(i-1,i-1));
    tV = tV.i();
    tM = tV*(U(i-1,0)-lam[0]);
    TEMP = -.5*tM.t()*solve(tV,tM);
    INTC(i,0) = TEMP(0,0);
    //Total positive area for the univariate normal distribution
    pr = 1-R::pnorm(0,tM(0,0),sqrt(tV(0,0)),true,false);
    TEMP = log(pr*sqrt(tV)*lam[0])+0.5*log(2.0*M_PI);
    INTC(i,1) = TEMP(0,0);
    uniMeans(i-1,0) = tM(0,0); uniVars(i-1,0) = tV(0,0);
    
  }
  //  bivariate values don't worry about looping as it just adds extra
  //  First  Condition
  
  
  PMAT(4,0) = p*(1-p)*(1-p);
  arma::mat tU = U.submat(0,0,1,0);
  arma::mat tB = B.submat(0,0,1,1);
  tLAM = LAM.submat(0,0,1,1);
  tV = (tB);
  tV = tV.i();
  tM = tV*(tU-lam[0]);
  bivVars.submat(0,0,1,1) = tV;
  bivMeans.submat(0,0,1,0) =  tM;
  TEMP = -.5*tM.t()*solve(tV,tM);
  INTC(4,0) = TEMP(0,0);
  pr = BivNProb(wrap(tM),wrap(tV));
  if (pr < 0 || isnan(pr)){pr = 0.0;} // rare numerical error
  TEMP = log(pr*(sqrt(abs(det(tV)))*lam[0]*lam[0]))+log(2.0*M_PI);
  INTC(4,1) = TEMP(0,0);
  
  // Second Condition
  
  PMAT(5,0) = p*(1-p)*(1-p);
  tB(0,0) = B(0,0);tB(0,1) = B(0,2);
  tB(1,0) = B(2,0);tB(1,1) = B(2,2);
  tLAM(0,0) = lam[0]; tLAM(1,1) = lam[2];
  tU(0,0) = U(0,0); tU(1,0) = U(2,0);
  tV = (tB);
  tV = tV.i();
  tM = tV*(tU-lam[0]);
  bivVars.submat(0,2,1,3) = tV;
  bivMeans.submat(0,1,1,1) = tM;
  TEMP = -.5*tM.t()*solve(tV,tM);
  INTC(5,0) = TEMP(0,0);
  
  pr = BivNProb(wrap(tM), wrap(tV));
  if (pr < 0 || isnan(pr)){pr = 0.0;} // rare numerical error
  TEMP = log(pr*(sqrt(abs(det(tV)))*lam[0]*lam[0]))+log(2.0*M_PI);
  INTC(5,1) = TEMP(0,0);
  
  // Third Condition
  PMAT(6,0) = p*(1-p)*(1-p);
  tU = U.submat(1,0,2,0);
  tB = B.submat(1,1,2,2); tLAM = LAM.submat(1,1,2,2);
  tV = (tB);
  tV = tV.i();
  tM = tV*(tU-lam[0]);
  bivVars.submat(0,4,1,5) = tV;
  bivMeans.submat(0,2,1,2) = tM;
  TEMP = -.5*tM.t()*solve(tV,tM);
  INTC(6,0) = TEMP(0,0);
  pr = BivNProb(wrap(tM), wrap(tV));
  if (pr < 0 || isnan(pr)){pr = 0.0;} // rare numerical error
  TEMP = log(pr*(sqrt(abs(det(tV)))*lam[0]*lam[0]))+log(2.0*M_PI);
  INTC(6,1) = TEMP(0,0);
  
  // TRIVARIATE CASE
  PMAT(7,0) = (1-p)*(1-p)*(1-p);
  tV = (B);
  tV = tV.i();
  tM = tV*(U-lam[0]);
  TEMP = -.5*tM.t()*solve(tV,tM);
  INTC(7,0) = TEMP(0,0);
  pr=   TriNProb(wrap(tM), wrap(tV));
  if (pr < 0 || isnan(pr)){pr = 0.0;} // rare numerical error
  TEMP = log(pr*sqrt(abs(det(tV)))*lam[0]*lam[0]*lam[0])+1.5*log(2.0*M_PI);
  INTC(7,1) = TEMP(0,0);
  
  
  
  // figure out the log Pbility  of the move
  
  
  arma::mat TT = INTC.col(1)-INTC.col(0);
  //cout << TT << endl; 
  arma::mat MT = max(TT);
  arma::mat SMULT = PMAT%exp(TT-repmat(MT,8,1));
  
  arma::mat log_prob = -0.5*(tau)*(Y.t()*Y) + log(sum(SMULT)) + MT;
  
  arma::mat POST_PROBS = exp(log(SMULT)-log(repmat(sum(SMULT),8,1))); //posterior sampling prob
  GetRNGstate();
  // get the cumulative  density
  double rvalue = R::runif(0,1); int rval = -1;
  
  for (int i = 1; i < 8; i++){
    POST_PROBS(i,0) +=  POST_PROBS(i-1,0);
    
  }
  
  for (int i = 0; i < 8 && rval ==-1; i++){
    if (  rvalue <= POST_PROBS(i,0)){
      rval = i;
    }
  }
  
  rval = rval+1;
  
  PutRNGstate();
  arma::mat rVar, rMeans;
  
  switch(rval){
  case 8:
    
    rVar = tV;
    rMeans = tM;
    break;
  case 1:
    
    rVar  = 0.0;
    rMeans = 0.0;
    break;
  case 2:
  case 3:
  case 4:
    
    
    rVar = uniVars(rval-2,0);
    rMeans = uniMeans(rval-2,0);
    break;
  default:
    rVar = bivVars.submat(0,(rval-5)*2,1,(rval-5)*2+1);
  rMeans = bivMeans.submat(0,rval-5,1,rval-5);
  break;
  }
  List returnV(4); 
  returnV[0] = wrap(rVar);
  returnV[1] = wrap(rMeans); 
  returnV[2] = wrap(log_prob(0,0)); 
  returnV[3] = wrap(rval);
  
  StringVector names(4);
  
  names(0) = "VAR"; 
  names(1) = "MEAN"; 
  names(2) = "LPROB";
  names(3) = "S_ELM"; 
  
  returnV.attr("names") = names;
  
  return returnV; 
}

// [[Rcpp::export]]
NumericVector shapesplineInsert(NumericVector k, NumericVector tp,NumericVector txi,
                                NumericVector tdeg, NumericVector tCBX,NumericVector tpos){
  
  arma::mat knots = Rcpp::as<arma::mat>(k); 
  arma::mat t     = Rcpp::as<arma::mat>(tp);
  
  
  int deg = tdeg[0];
  int pos = tpos[0]-1; // adjust  the position for C++ and extra knots
  knots.insert_rows(0,1); knots.insert_rows(0,1); 
  knots.insert_rows(knots.n_rows,1); knots.insert_rows(knots.n_rows,1);
  knots(0,0) = knots(2,0); knots(1,0) = knots(2,0);
  knots(knots.n_rows-1,0) = knots(knots.n_rows-3,0);
  knots(knots.n_rows-2,0) = knots(knots.n_rows-3,0); 
  
  IntegerVector arrayDims = tCBX.attr("dim");
  arma::cube CBX(tCBX.begin(),arrayDims[0],arrayDims[1],arrayDims[2]);
  arma::cube RVAL = arma::cube(arrayDims[0],arrayDims[1]+1,arrayDims[2]);
  
  RVAL.zeros();
  //      Rprintf("%d %d %d %d",pos,TV.n_rows,TV.n_cols,TV.n_slices);
  RVAL.tube(0,0,RVAL.n_rows-1,pos-1) = CBX.tube(0,0,CBX.n_rows-1,pos-1);
  RVAL.tube(0,pos+1,RVAL.n_rows-1,RVAL.n_cols-1) = CBX.tube(0,pos,CBX.n_rows-1,CBX.n_cols-1);
  
  for (int j = 0; j <= deg ; j++){
    for (int i = pos-1;  i  <=  pos+2; i++){
      double h1 = knots(i,0); double h2 = knots(i+1,0); 
      double h3 = knots(i+2,0); double h4 = knots(i+3,0);
      
      double    t1=(((h2-h1)*(h3-h1))); 
      double    t2=(((h3-h1)*(h3-h2)));      
      double    t3=(((h4-h2)*(h3-h2))); 
      double    t4=(((h4-h3)*(h4-h2)));      
      double    z1,z2,z3,z4; 
      if (t1 == 0){ z1 = 0;} else{z1 = 1/t1;}
      if (t2 == 0){ z2 = 0;} else{z2 = 1/t2;}
      if (t3 == 0){ z3 = 0;} else{z3 = 1/t3;}
      if (t4 == 0){ z4 = 0;} else{z4 = 1/t4;}
      
      
      arma::umat LTH1 = (t <= h1);
      arma::umat LTH2 = (t <= h2);
      arma::umat LTH3 = (t <= h3);
      arma::umat LTH4 = (t <= h4); 
      
      
      arma::mat S1 = 1.0/double(3+j)*pow(t,3+j)-2.0*h1/(2+j)*pow(t,2+j)+pow(t,1+j)*pow(h1,2)/double(1+j); S1 = z1*S1; 
      arma::mat BS1= S1*0 + ( 1.0/double(3+j)*pow(h1,3+j)-2.0*h1/double(2+j)*pow(h1,2+j)+(pow(h1,1+j))*pow(h1,2)/double(1+j))*z1;  
      arma::mat AS1 = S1*0 + ( 1.0/double(3+j)*pow(h2,3+j)-2.0*h1/double(2+j)*pow(h2,2+j)+(pow(h2,1+j))*pow(h1,2)/double(1+j))*z1;
      
      arma::mat S2 = h3/double(2+j)*pow(t,2+j) - 1.0/double(3+j)*pow(t,3+j)-1.0/double(1+j)*pow(t,j+1)*h1*h3+h1/double(2+j)*pow(t,2+j); S2 = S2*z2; 
      arma::mat BS2 = S2*0 + (h3/double(2+j)*pow(h2,2+j) - 1.0/double(3+j)*pow(h2,3+j)-1.0/double(1+j)*pow(h2,j+1)*h1*h3+h1/double(2+j)*pow(h2,2+j))*z2;
      arma::mat AS2 = S2*0 + (h3/double(2+j)*pow(h3,2+j) - 1.0/double(3+j)*pow(h3,3+j)-1.0/double(1+j)*pow(h3,j+1)*h1*h3+h1/double(2+j)*pow(h3,2+j))*z2;
      
      arma::mat S3 = h4/double(2+j)*pow(t,2+j) - 1.0/double(3+j)*pow(t,3+j) - h2*h4/double(1+j)*pow(t,1+j) + h2/double(2+j)*pow(t,2+j); S3 = S3*z3;
      arma::mat BS3 = S3*0 + (h4/double(2+j)*pow(h2,2+j) - 1.0/double(3+j)*pow(h2,3+j) - h2*h4/double(1+j)*pow(h2,1+j) + h2/double(2+j)*pow(h2,2+j) )*z3;
      AS2 = AS2 + (h4/double(2+j)*pow(h3,2+j) - 1.0/double(3+j)*pow(h3,3+j) - h2*h4/double(1+j)*pow(h3,1+j) + h2/double(2+j)*pow(h3,2+j) )*z3;
      
      arma::mat S4 = 1.0/double(1+j)*pow(h4,2)*pow(t,1+j) - 2.0*h4/double(2+j)*pow(t,2+j) + 1.0/double(3+j)*pow(t,3+j); S4 = S4*z4; 
      arma::mat BS4 = S4*0 + (1.0/double(1+j)*pow(h4,2)*pow(h3,1+j) - 2.0*h4/double(2+j)*pow(h3,2+j) + 1.0/double(3+j)*pow(h3,3+j))*z4; 
      arma::mat AS3 = S4*0 + (1.0/double(1+j)*pow(h4,2)*pow(h4,1+j) - 2.0*h4/double(2+j)*pow(h4,2+j) + 1.0/double(3+j)*pow(h4,3+j))*z4; 
      
      arma::mat A  = S1*0;
      A  = A + ((S1-BS1)%LTH2 + (AS1-BS1)%(1-LTH2))%(1-LTH1);
      A  = A + ((S2+S3-(BS2+BS3))%LTH3 + (AS2-(BS2+BS3))%(1-LTH3))%(1-LTH2);
      A  = A + ((S4-BS4)%LTH4 + (AS3-BS4)%(1-LTH4))%(1-LTH3) ;
      RVAL.slice(j).submat(arma::span(),arma::span(i)) = A; 
    }
  }
  
  
  return wrap(RVAL); 
  
  // return R_NilValue;
  
}

// [[Rcpp::export]]
NumericVector shapesplineDelete(NumericVector k, NumericVector tp,NumericVector txi,
                                NumericVector tdeg, NumericVector tCBX,NumericVector tpos){
  arma::mat knots = Rcpp::as<arma::mat>(k); 
  arma::mat t     = Rcpp::as<arma::mat>(tp);
  
  
  
  
  int deg = tdeg[0];
  int pos = tpos[0] - 1; // adjust  the position for C++ indexing
  
  knots.insert_rows(0,1); knots.insert_rows(0,1); 
  knots.insert_rows(knots.n_rows,1); knots.insert_rows(knots.n_rows,1);
  knots(0,0) = knots(2,0); knots(1,0) = knots(2,0);
  knots(knots.n_rows-1,0) = knots(knots.n_rows-3,0);
  knots(knots.n_rows-2,0) = knots(knots.n_rows-3,0); 
  
  IntegerVector arrayDims = tCBX.attr("dim");
  arma::cube CBX(tCBX.begin(),arrayDims[0],arrayDims[1],arrayDims[2]);
  arma::cube RVAL = arma::cube(arrayDims[0],arrayDims[1]-1,arrayDims[2]);
  
  RVAL.zeros();
  //      Rprintf("%d %d %d %d",pos,TV.n_rows,TV.n_cols,TV.n_slices);
  RVAL.tube(0,0,RVAL.n_rows-1,pos-1) = CBX.tube(0,0,CBX.n_rows-1,pos-1);
  RVAL.tube(0,pos,RVAL.n_rows-1,RVAL.n_cols-1) = CBX.tube(0,pos+1,CBX.n_rows-1,CBX.n_cols-1);
  
  for (int j = 0; j <= deg ; j++){
    for (int i = pos-1;  i  <=  pos+1; i++){
      double h1 = knots(i,0); double h2 = knots(i+1,0); 
      double h3 = knots(i+2,0); double h4 = knots(i+3,0);
      
      double    t1=(((h2-h1)*(h3-h1))); 
      double    t2=(((h3-h1)*(h3-h2)));      
      double    t3=(((h4-h2)*(h3-h2))); 
      double    t4=(((h4-h3)*(h4-h2)));      
      double    z1,z2,z3,z4; 
      if (t1 == 0){ z1 = 0;} else{z1 = 1/t1;}
      if (t2 == 0){ z2 = 0;} else{z2 = 1/t2;}
      if (t3 == 0){ z3 = 0;} else{z3 = 1/t3;}
      if (t4 == 0){ z4 = 0;} else{z4 = 1/t4;}
      
      
      arma::umat LTH1 = (t <= h1);
      arma::umat LTH2 = (t <= h2);
      arma::umat LTH3 = (t <= h3);
      arma::umat LTH4 = (t <= h4); 
      
      
      arma::mat S1 = 1.0/double(3+j)*pow(t,3+j)-2.0*h1/(2+j)*pow(t,2+j)+pow(t,1+j)*pow(h1,2)/double(1+j); S1 = z1*S1; 
      arma::mat BS1= S1*0 + ( 1.0/double(3+j)*pow(h1,3+j)-2.0*h1/double(2+j)*pow(h1,2+j)+(pow(h1,1+j))*pow(h1,2)/double(1+j))*z1;  
      arma::mat AS1 = S1*0 + ( 1.0/double(3+j)*pow(h2,3+j)-2.0*h1/double(2+j)*pow(h2,2+j)+(pow(h2,1+j))*pow(h1,2)/double(1+j))*z1;
      
      arma::mat S2 = h3/double(2+j)*pow(t,2+j) - 1.0/double(3+j)*pow(t,3+j)-1.0/double(1+j)*pow(t,j+1)*h1*h3+h1/double(2+j)*pow(t,2+j); S2 = S2*z2; 
      arma::mat BS2 = S2*0 + (h3/double(2+j)*pow(h2,2+j) - 1.0/double(3+j)*pow(h2,3+j)-1.0/double(1+j)*pow(h2,j+1)*h1*h3+h1/double(2+j)*pow(h2,2+j))*z2;
      arma::mat AS2 = S2*0 + (h3/double(2+j)*pow(h3,2+j) - 1.0/double(3+j)*pow(h3,3+j)-1.0/double(1+j)*pow(h3,j+1)*h1*h3+h1/double(2+j)*pow(h3,2+j))*z2;
      
      arma::mat S3 = h4/double(2+j)*pow(t,2+j) - 1.0/double(3+j)*pow(t,3+j) - h2*h4/double(1+j)*pow(t,1+j) + h2/double(2+j)*pow(t,2+j); S3 = S3*z3;
      arma::mat BS3 = S3*0 + (h4/double(2+j)*pow(h2,2+j) - 1.0/double(3+j)*pow(h2,3+j) - h2*h4/double(1+j)*pow(h2,1+j) + h2/double(2+j)*pow(h2,2+j) )*z3;
      AS2 = AS2 + (h4/double(2+j)*pow(h3,2+j) - 1.0/double(3+j)*pow(h3,3+j) - h2*h4/double(1+j)*pow(h3,1+j) + h2/double(2+j)*pow(h3,2+j) )*z3;
      
      arma::mat S4 = 1.0/double(1+j)*pow(h4,2)*pow(t,1+j) - 2.0*h4/double(2+j)*pow(t,2+j) + 1.0/double(3+j)*pow(t,3+j); S4 = S4*z4; 
      arma::mat BS4 = S4*0 + (1.0/double(1+j)*pow(h4,2)*pow(h3,1+j) - 2.0*h4/double(2+j)*pow(h3,2+j) + 1.0/double(3+j)*pow(h3,3+j))*z4; 
      arma::mat AS3 = S4*0 + (1.0/double(1+j)*pow(h4,2)*pow(h4,1+j) - 2.0*h4/double(2+j)*pow(h4,2+j) + 1.0/double(3+j)*pow(h4,3+j))*z4; 
      
      arma::mat A  = S1*0;
      A  = A + ((S1-BS1)%LTH2 + (AS1-BS1)%(1-LTH2))%(1-LTH1);
      A  = A + ((S2+S3-(BS2+BS3))%LTH3 + (AS2-(BS2+BS3))%(1-LTH3))%(1-LTH2);
      A  = A + ((S4-BS4)%LTH4 + (AS3-BS4)%(1-LTH4))%(1-LTH3) ;
      RVAL.slice(j).submat(arma::span(),arma::span(i)) = A; 
    }
  }
  
  return wrap(RVAL);
  //return R_NilValue;
  
}


double rtn(double mean, double sd){
  double rV = 0.0;
  
  if (mean < 0){
    // do everything on the log scale  - more numerically stable
    double ucutoff = -1.0*R::pnorm(0.0,-1.0*(mean),sd,true, true);
    GetRNGstate();
    double cv = R::rexp(1) + ucutoff;
    PutRNGstate();
    rV = -1.0*R::qnorm(-1.0*cv,-1.0*(mean),sd,true,true);
  }else{
    double lcutoff = R::pnorm(0.0,mean,sd,true, false);
    GetRNGstate();
    double tunif = R::runif(lcutoff,1.0);
    PutRNGstate();
    rV = R::qnorm(tunif,mean,sd,true,false);
  }
  
  return  rV;
  
}
// [[Rcpp::export]]
NumericVector rtmvn(NumericVector tMean, NumericVector tVar){
  arma::mat COV   = Rcpp::as<arma::mat>(tVar);
  arma::mat MEAN  = Rcpp::as<arma::mat>(tMean);
  arma::mat tV12_22inv_tM2(1,MEAN.n_rows);
  arma::mat tV(1,MEAN.n_rows);
  arma::mat tV12_22inv(1,MEAN.n_rows*(MEAN.n_rows-1));
  arma::mat RV = MEAN;
  
  
  
  if (MEAN.n_rows == 1){
    RV(0,0) =  rtn(MEAN(0,0),sqrt(COV(0,0)));
  }else{
    for (int i = 0; i < int(MEAN.n_rows);i++){
      arma::mat TEMP = COV;
      arma::mat TEMP2 = COV.cols(i,i);
      arma::mat TEMP3 = MEAN;
      
      TEMP3.shed_row(i);
      TEMP.shed_row(i);
      TEMP.shed_col(i);
      
      TEMP2.shed_row(i);
      tV12_22inv.submat(0,i*(MEAN.n_rows-1),0,(i+1)*(MEAN.n_rows-1)-1) = TEMP2.t()*TEMP.i();
      tV.submat(0,i,0,i) =  sqrt(COV.submat(i,i,i,i) - tV12_22inv.submat(0,i*(MEAN.n_rows-1),0,(i+1)*(MEAN.n_rows-1)-1)*TEMP2); // variances conditional on the mean
      tV12_22inv_tM2.submat(0,i,0,i)  =  MEAN.submat(i,0,i,0) - TEMP2.t()*solve(TEMP,TEMP3);
    }
    // set the starting value
    
    for (int i = 0 ; i< int(MEAN.n_rows);i++){
      if (RV(i,0) < 0.0){
        RV(i,0) = 0.01;
      }
    }
    
    arma::mat TempV = RV;
    
    
    for (int i = 0; i < 20; i++){
      // gibbs sampler
      for (int j = 0; j < int(MEAN.n_rows); j++){
        TempV = RV;
        TempV.shed_row(j);
        TempV = tV12_22inv_tM2.submat(0,j,0,j)+ tV12_22inv.submat(0,j*(MEAN.n_rows-1),0,(j+1)*(MEAN.n_rows-1)-1)*TempV;
        RV(j,0) = rtn(TempV(0,0),tV(0,j));
      }
    }
    
  }
  
  
  // return the R value
  return wrap(RV);
}

// [[Rcpp::export]]
SEXP qcopy(SEXP a,SEXP b, SEXP c,SEXP d){
  double *x =  REAL(a);
  double *y =  REAL(b);
  int n =  INTEGER(c)[0];
  int i =  INTEGER(d)[0];
  memmove(&x[n*(i-1)],y,n*sizeof(double));
  return R_NilValue;
}

NumericVector Tnorm(NumericVector mean, NumericVector sd); 
RcppExport SEXP sourceCpp_0_Tnorm(SEXP meanSEXP, SEXP sdSEXP) {
  BEGIN_RCPP
  Rcpp::RObject __result;
  Rcpp::RNGScope __rngScope;
  Rcpp::traits::input_parameter< NumericVector >::type mean(meanSEXP);
  Rcpp::traits::input_parameter< NumericVector >::type sd(sdSEXP);
  __result = Rcpp::wrap(Tnorm(mean, sd));
  return __result;
  END_RCPP
}

SEXP qcopy(SEXP a, SEXP b, SEXP c, SEXP d);
RcppExport SEXP sourceCpp_0_qcopy(SEXP aSEXP, SEXP bSEXP, SEXP cSEXP, SEXP dSEXP) {
  BEGIN_RCPP
  Rcpp::RObject __result;
  Rcpp::RNGScope __rngScope;
  Rcpp::traits::input_parameter< SEXP >::type a(aSEXP);
  Rcpp::traits::input_parameter< SEXP >::type b(bSEXP);
  Rcpp::traits::input_parameter< SEXP >::type c(cSEXP);
  Rcpp::traits::input_parameter< SEXP >::type d(dSEXP);
  __result = Rcpp::wrap(qcopy(a, b, c, d));
  return __result;
  END_RCPP
}

double BivNProb(NumericVector mean, NumericVector cv);
RcppExport SEXP sourceCpp_0_BivNProb(SEXP meanSEXP, SEXP cvSEXP) {
  BEGIN_RCPP
  Rcpp::RObject __result;
  Rcpp::RNGScope __rngScope;
  Rcpp::traits::input_parameter< NumericVector >::type mean(meanSEXP);
  Rcpp::traits::input_parameter< NumericVector >::type cv(cvSEXP);
  __result = Rcpp::wrap(BivNProb(mean, cv));
  return __result;
  END_RCPP
}
// TriNProb
double TriNProb(NumericVector mean, NumericVector cv);
RcppExport SEXP sourceCpp_0_TriNProb(SEXP meanSEXP, SEXP cvSEXP) {
  BEGIN_RCPP
  Rcpp::RObject __result;
  Rcpp::RNGScope __rngScope;
  Rcpp::traits::input_parameter< NumericVector >::type mean(meanSEXP);
  Rcpp::traits::input_parameter< NumericVector >::type cv(cvSEXP);
  __result = Rcpp::wrap(TriNProb(mean, cv));
  return __result;
  END_RCPP
}
// sampleBetas
NumericVector sampleBetas(NumericVector ttY, NumericVector ttX, NumericVector tbetas, NumericVector LAM, NumericVector intLAM, NumericVector p, NumericVector tau);
RcppExport SEXP sourceCpp_0_sampleBetas(SEXP ttYSEXP, SEXP ttXSEXP, SEXP tbetasSEXP, SEXP LAMSEXP, SEXP intLAMSEXP, SEXP pSEXP, SEXP tauSEXP) {
  BEGIN_RCPP
  Rcpp::RObject __result;
  Rcpp::RNGScope __rngScope;
  Rcpp::traits::input_parameter< NumericVector >::type ttY(ttYSEXP);
  Rcpp::traits::input_parameter< NumericVector >::type ttX(ttXSEXP);
  Rcpp::traits::input_parameter< NumericVector >::type tbetas(tbetasSEXP);
  Rcpp::traits::input_parameter< NumericVector >::type LAM(LAMSEXP);
  Rcpp::traits::input_parameter< NumericVector >::type intLAM(intLAMSEXP);
  Rcpp::traits::input_parameter< NumericVector >::type p(pSEXP);
  Rcpp::traits::input_parameter< NumericVector >::type tau(tauSEXP);
  __result = Rcpp::wrap(sampleBetas(ttY, ttX, tbetas, LAM, intLAM, p, tau));
  return __result;
  END_RCPP
}
// sinsertBeta
List sinsertBeta(NumericVector tY, NumericVector Xidx, NumericVector ctau, NumericVector tp, NumericVector lam);
RcppExport SEXP sourceCpp_0_sinsertBeta(SEXP tYSEXP, SEXP XidxSEXP, SEXP ctauSEXP, SEXP tpSEXP, SEXP lamSEXP) {
  BEGIN_RCPP
  Rcpp::RObject __result;
  Rcpp::RNGScope __rngScope;
  Rcpp::traits::input_parameter< NumericVector >::type tY(tYSEXP);
  Rcpp::traits::input_parameter< NumericVector >::type Xidx(XidxSEXP);
  Rcpp::traits::input_parameter< NumericVector >::type ctau(ctauSEXP);
  Rcpp::traits::input_parameter< NumericVector >::type tp(tpSEXP);
  Rcpp::traits::input_parameter< NumericVector >::type lam(lamSEXP);
  __result = Rcpp::wrap(sinsertBeta(tY, Xidx, ctau, tp, lam));
  return __result;
  END_RCPP
}
// sdeleteBeta
List sdeleteBeta(NumericVector tY, NumericVector Xidx, NumericVector ctau, NumericVector tp, NumericVector lam);
RcppExport SEXP sourceCpp_0_sdeleteBeta(SEXP tYSEXP, SEXP XidxSEXP, SEXP ctauSEXP, SEXP tpSEXP, SEXP lamSEXP) {
  BEGIN_RCPP
  Rcpp::RObject __result;
  Rcpp::RNGScope __rngScope;
  Rcpp::traits::input_parameter< NumericVector >::type tY(tYSEXP);
  Rcpp::traits::input_parameter< NumericVector >::type Xidx(XidxSEXP);
  Rcpp::traits::input_parameter< NumericVector >::type ctau(ctauSEXP);
  Rcpp::traits::input_parameter< NumericVector >::type tp(tpSEXP);
  Rcpp::traits::input_parameter< NumericVector >::type lam(lamSEXP);
  __result = Rcpp::wrap(sdeleteBeta(tY, Xidx, ctau, tp, lam));
  return __result;
  END_RCPP
}

// shapesplineInsert2
NumericVector shapesplineInsert(NumericVector k, NumericVector tp, NumericVector txi, NumericVector tdeg, NumericVector tCBX, NumericVector tpos);
RcppExport SEXP sourceCpp_0_shapesplineInsert(SEXP kSEXP, SEXP tpSEXP, SEXP txiSEXP, SEXP tdegSEXP, SEXP tCBXSEXP, SEXP tposSEXP) {
  BEGIN_RCPP
  Rcpp::RObject __result;
  Rcpp::RNGScope __rngScope;
  Rcpp::traits::input_parameter< NumericVector >::type k(kSEXP);
  Rcpp::traits::input_parameter< NumericVector >::type tp(tpSEXP);
  Rcpp::traits::input_parameter< NumericVector >::type txi(txiSEXP);
  Rcpp::traits::input_parameter< NumericVector >::type tdeg(tdegSEXP);
  Rcpp::traits::input_parameter< NumericVector >::type tCBX(tCBXSEXP);
  Rcpp::traits::input_parameter< NumericVector >::type tpos(tposSEXP);
  __result = Rcpp::wrap(shapesplineInsert(k, tp, txi, tdeg, tCBX, tpos));
  return __result;
  END_RCPP
}
// shapesplineDelete2
NumericVector shapesplineDelete(NumericVector k, NumericVector tp, NumericVector txi, NumericVector tdeg, NumericVector tCBX, NumericVector tpos);
RcppExport SEXP sourceCpp_0_shapesplineDelete(SEXP kSEXP, SEXP tpSEXP, SEXP txiSEXP, SEXP tdegSEXP, SEXP tCBXSEXP, SEXP tposSEXP) {
  BEGIN_RCPP
  Rcpp::RObject __result;
  Rcpp::RNGScope __rngScope;
  Rcpp::traits::input_parameter< NumericVector >::type k(kSEXP);
  Rcpp::traits::input_parameter< NumericVector >::type tp(tpSEXP);
  Rcpp::traits::input_parameter< NumericVector >::type txi(txiSEXP);
  Rcpp::traits::input_parameter< NumericVector >::type tdeg(tdegSEXP);
  Rcpp::traits::input_parameter< NumericVector >::type tCBX(tCBXSEXP);
  Rcpp::traits::input_parameter< NumericVector >::type tpos(tposSEXP);
  __result = Rcpp::wrap(shapesplineDelete(k, tp, txi, tdeg, tCBX, tpos));
  return __result;
  END_RCPP
}

// rtmvn
NumericVector rtmvn(NumericVector tMean, NumericVector tVar);
RcppExport SEXP sourceCpp_0_rtmvn(SEXP tMeanSEXP, SEXP tVarSEXP) {
  BEGIN_RCPP
  Rcpp::RObject __result;
  Rcpp::RNGScope __rngScope;
  Rcpp::traits::input_parameter< NumericVector >::type tMean(tMeanSEXP);
  Rcpp::traits::input_parameter< NumericVector >::type tVar(tVarSEXP);
  __result = Rcpp::wrap(rtmvn(tMean, tVar));
  return __result;
  END_RCPP
}

double SADMVN(arma::mat M, arma::mat C);
RcppExport SEXP sourceCpp_0_SADMVN(SEXP MSEXP, SEXP CSEXP) {
  BEGIN_RCPP
  Rcpp::RObject __result;
  Rcpp::RNGScope __rngScope;
  Rcpp::traits::input_parameter< arma::mat >::type M(MSEXP);
  Rcpp::traits::input_parameter< arma::mat >::type C(CSEXP);
  __result = Rcpp::wrap(SADMVN(M, C));
  return __result;
  END_RCPP
}
