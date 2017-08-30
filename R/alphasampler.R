sampA<-function(Y,w,tau,E,V,amin,amax){
#Sample the A matrix given a Y and a precomputed W matrix
#which is based upon the values 

  EV = 1/(tau*(t(w)%*%w + 1));   
  ES = EV*(tau*t(w)%*%Y ); 
  UB = pnorm(1,ES,sqrt(EV)); 
  LB = pnorm(-1,ES,sqrt(EV)); 
  alpha = qnorm(runif(1,LB,UB),ES,sqrt(EV));           
 
  return(alpha);                                             
}

