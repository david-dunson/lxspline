#library('mnormt')
#library('tmvtnorm')

## 
#Function: shapeSpline: computes a various factors from a shape spline
#matrix
#Output : X  a length(t) x length (knots) x deg+1 Array
# containing the necessary information to calculate a shapespline matrix
# given alpha
##
shapespline <- function(knots, t,xi,deg){
  knots = sort(knots); 
  ti = length(knots); 
  knots = c( knots[1] - (knots[2]-knots[1]), knots,knots[ti]+(knots[ti]-knots[ti-1]) ); 

  X = array(0,dim=c(length(t),length(knots)-2,deg+1));
  
  for (j in 0:deg){
    for ( i in 2:(length(knots)-1)){
      h1 = knots[i-1]
      h2 = knots[i]
      h3 = knots[i+1];  
      z1=(2/((h3-h1)*(h2-h1))); 
      z2=(2/((h3-h1)*(h3-h2)));      
      
      #vector if t is less than the 
      #given knot h1 through h3
      LTH1 = (t <= h1);
      LTH2 = (t <= h2); 
      LTH3 = (t <= h3); 
      
      # perform the necessary integrations for the triangle
      # distribution
      S1 = z1*(1/(j+2)*h2^(j+2)- h1/(j+1)*h2^(j+1)+ (1/(j+1)-1/(j+2))*h1^(j+2)); 
      S2 = z2*(h3/(j+1)*h3^(j+1)-1/(j+2)*h3^(j+2) + (h2^(j+2)/(j+2)-h3*h2^(j+1)/(j+1)));
      
      A1 = z1*(1/(j+2)*t^(j+2)- h1/(j+1)*t^(j+1)+ (1/(j+1)-1/(j+2))*h1^(j+2)); 
      A2 = S1 + z2*(h3/(j+1)*t^(j+1)-1/(j+2)*t^(j+2) + (h2^(j+2)/(j+2)-h3*h2^(j+1)/(j+1))); 
      
      X[,i-1,j+1] = ((1-LTH1)*LTH2*A1 + (1-LTH1)*(1-LTH2)*LTH3*A2 + (1-LTH1)*(1-LTH2)*(1-LTH3)*(S1+S2)) #*(1/h3)^deg; 
     
    }
  } 
  return(X)
}

shapespline2 <- function(knots,t,deg){
  X = array(0,dim=c(length(t),length(knots)+1,deg+1));
  knots = c(min(knots),min(knots), knots, max(knots),max(knots)); 
  
  
  for (j in 0:deg){
    for ( i in 1:(length(knots)-3)){
      h1 = knots[i]
      h2 = knots[i+1]
      h3 = knots[i+2] 
      h4 = knots[i+3]
      
      t1=(((h2-h1)*(h3-h1))); 
      t2=(((h3-h1)*(h3-h2)));      
      t3=(((h4-h2)*(h3-h2))); 
      t4=(((h4-h3)*(h4-h2)));      
      if (t1 == 0){ z1 = 0} else{z1 = 1/t1}
      if (t2 == 0){ z2 = 0} else{z2 = 1/t2}
      if (t3 == 0){ z3 = 0} else{z3 = 1/t3}
      if (t4 == 0){ z4 = 0} else{z4 = 1/t4}
      
      #vector if t is less than the 
      #given knot h1 through h4
      LTH1 = (t <= h1);
      LTH2 = (t <= h2); 
      LTH3 = (t <= h3); 
      LTH4 = (t <= h4); 
      
      # perform the necessary integrations
      
      S1 = 1/(3+j)*t^(3+j)-2*h1/(2+j)*t^(2+j)+t^(1+j)*h1^2/(1+j); S1 = z1*S1; 
      BS1= S1*0 + ( 1/(3+j)*h1^(3+j)-2*h1/(2+j)*h1^(2+j)+(h1^(1+j))*h1^2/(1+j))*z1  
      AS1 = S1*0 + ( 1/(3+j)*h2^(3+j)-2*h1/(2+j)*h2^(2+j)+(h2^(1+j))*h1^2/(1+j))*z1
      
      S2 = h3/(2+j)*t^(2+j) - 1/(3+j)*t^(3+j)-1/(1+j)*t^(j+1)*h1*h3+h1/(2+j)*t^(2+j); S2 = S2*z2; 
      BS2 = S2*0 + (h3/(2+j)*h2^(2+j) - 1/(3+j)*h2^(3+j)-1/(1+j)*h2^(j+1)*h1*h3+h1/(2+j)*h2^(2+j))*z2;
      AS2 = S2*0 + (h3/(2+j)*h3^(2+j) - 1/(3+j)*h3^(3+j)-1/(1+j)*h3^(j+1)*h1*h3+h1/(2+j)*h3^(2+j))*z2;
      
      S3 = h4/(2+j)*t^(2+j) - 1/(3+j)*t^(3+j) - h2*h4/(1+j)*t^(1+j) + h2/(2+j)*t^(2+j); S3 = S3*z3;
      BS3 = S3*0 + (h4/(2+j)*h2^(2+j) - 1/(3+j)*h2^(3+j) - h2*h4/(1+j)*h2^(1+j) + h2/(2+j)*h2^(2+j) )*z3;
      AS2 = AS2 + (h4/(2+j)*h3^(2+j) - 1/(3+j)*h3^(3+j) - h2*h4/(1+j)*h3^(1+j) + h2/(2+j)*h3^(2+j) )*z3;
      
      S4 = 1/(1+j)*h4^2*t^(1+j) - 2*h4/(2+j)*t^(2+j) + 1/(3+j)*t^(3+j); S4 = S4*z4; 
      BS4 = S4*0 + (1/(1+j)*h4^2*h3^(1+j) - 2*h4/(2+j)*h3^(2+j) + 1/(3+j)*h3^(3+j))*z4; 
      AS3 = S4*0 + (1/(1+j)*h4^2*h4^(1+j) - 2*h4/(2+j)*h4^(2+j) + 1/(3+j)*h4^(3+j))*z4; 
      
      ##Now you do the fundamental theorem on each
      ##portion to integrate it all up
      A  = S1*0
      A  = A + ((S1-BS1)*LTH2 + (AS1-BS1)*(1-LTH2))*(1-LTH1)
      A  = A + ((S2+S3-(BS2+BS3))*LTH3 + (AS2-(BS2+BS3))*(1-LTH3))*(1-LTH2)
      A  = A + ((S4-BS4)*LTH4 + (AS3-BS4)*(1-LTH4))*(1-LTH3) 
      
      
      X[,i,j+1] = A ; 
      
    }
  } 
  return(X)
}
