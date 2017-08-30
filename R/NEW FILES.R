#generate a tree
m.TREE = matrix(-1,nrow=6,3)


m.TREE[1,1:3] = as.matrix(c(0,0.5,1))
m.TREE[2,1:3] = as.matrix(c(0,0.25,0))
m.TREE[3,1:3] = as.matrix(c(0,0.75,0))
m.TREE[4,1:3] = as.matrix(c(0,0.25,0)) 
m.TREE[5,1:3] = as.matrix(c(0, 1, 0))
m.TREE[6,1:3] = as.matrix(c(0, 1, 0))


PROB.PROPOSAL <- function(m.TREE,n.TREE){
  new.C = ncol(n.TREE)
  old.C = ncol(m.TREE)
  
  prob = .5; 
  if (new.C > old.C){
        if (old.C == 3){
          return(prob)
        }
        return(.5*(1/sum(m.TREE[5:6,])));  
  }else{
        return(.5*(1/sum(m.TREE[5,]*m.TREE[6,])));    
  }
}


PROP.N.TREE <- function(m.TREE){
  INSERT = T
  LEFT = F; 
  if (ncol(m.TREE)==3){ #must insert a node at 0.5
    if (runif(1)<0.5){
       LEFT = T; 
       m.TREE = TREE.INSERT(2,T,m.TREE)
       pos = 2; 
    }else{
       m.TREE = TREE.INSERT(2,F,m.TREE)     
       pos = 3; 
    }
    
  }else{
      temp = 1:ncol(m.TREE); 
      if (runif(1) < 0.5){ 
        lm = sum(m.TREE[5,])
        rm = sum(m.TREE[6,])
        n = lm + rm; 
        ts = sample(1:n,1); 
        if (ts <= lm){ # left insert
            LEFT = T; 
            pos = temp*m.TREE[5,]        
            pos = pos[pos > 0]
            m.TREE = TREE.INSERT(pos[ts],T,m.TREE)
            pos = pos[ts]
        }else{ # right insert 
            pos = temp*m.TREE[6,]        
            pos = pos[pos > 0]
            m.TREE = TREE.INSERT(pos[ts-lm],F,m.TREE)
            pos = pos[ts - lm] +1 
        }
      }else{
        INSERT = F; 
        available = m.TREE[5,]*m.TREE[6,]; 
        n = sum(available)
        ts = sample(1:n,1); 
        pos = available*temp; 
        pos = pos[pos > 0]
        m.TREE = TREE.DELETE(pos[ts],m.TREE)
        pos = pos[ts]
      }
  }
    
  return(list(INSERT=INSERT,m.TREE=m.TREE,pos=pos,LEFT=LEFT))
}

TREE.INSERT <- function(pos,left,m.TREE){
    #basic error checking 
    if (pos == 1 || pos == ncol(m.TREE)){
      return( m.TREE); 
    }
      
    new.col = matrix(0,nrow=6,1)
    new.col[5] = 1; new.col[6] = 1; 
    new.col[4] = m.TREE[4,pos]/2; 
    
    if (left == T){
      if (m.TREE[5,pos] == 0){ #cant add the tree
        return(m.TREE); 
      } 
      m.TREE[5,pos] = 0; 
      new.col[1] = m.TREE[1,pos]-m.TREE[4,pos]
      new.col[2] = new.col[1] - new.col[4]
      new.col[3] = new.col[1] + new.col[4]
      m.TREE = cbind(m.TREE[,1:(pos-1)],new.col,m.TREE[,pos:ncol(m.TREE)])
    }else{
      if (m.TREE[6,pos] == 0){ #cant add the tree
        return(m.TREE); 
      }
      m.TREE[6,pos] = 0;       
      new.col[1] = m.TREE[1,pos]+m.TREE[4,pos]
      new.col[2] = new.col[1] - new.col[4]
      new.col[3] = new.col[1] + new.col[4]
      m.TREE = cbind(m.TREE[,1:(pos)],new.col,m.TREE[,(pos+1):ncol(m.TREE)])
    }
    
    return(m.TREE)
}

TREE.DELETE <- function(pos,m.TREE){
  #basic error checking 
  if (pos == 1 || pos == ncol(m.TREE)){
    return( m.TREE); 
  }
  if (m.TREE[5,pos] == 0 || m.TREE[6,pos] == 0){
    return(m.TREE); 
  }
  
  if (m.TREE[1,pos] == m.TREE[3,pos-1]){
    m.TREE[6,pos-1] = 1; #right node now open
  }else{
    m.TREE[5,pos+1] = 1; #left node now open
  }
  
  m.TREE = cbind(m.TREE[,1:(pos-1)], m.TREE[,(pos+1):ncol(m.TREE)])
  return(m.TREE); 
}

LOG.PROB.TREE <-function (p,m.TREE){

  p.temp = p^(log(1/m.TREE[4,2:(ncol(m.TREE)-1)]))
  succ <- colSums(1-m.TREE[5:6,2:(ncol(m.TREE)-1),drop=F])
  fail <- colSums(m.TREE[5:6,2:(ncol(m.TREE)-1),drop=F]) 
  return(sum(log(p.temp^(succ)*(1-p.temp)^(fail))))
}

