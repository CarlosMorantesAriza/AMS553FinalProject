henry_model=function(size_a=1, size_b=1, size_c=1, size_d=1, size_e=1, b_factor=seq(0.1, 0.5, 0.1), d_prop=0.1,
                     runtime=1000, dim=100){
  set.seed(1)
  ticker=0
  sizes = c(size_a, size_b,size_c,size_d, size_e)
  b = b_factor / sizes #birth chance
  d = b_factor*d_prop / sizes #death chance
  
  death <- function(spec) {
    return(rbernoulli(1, p = 1 - d[spec]) * spec)
  }
  
  birth <- function(neighbors) {
    neighbors = if(length(which(neighbors==0)!=0)) neighbors[-which(neighbors==0)] # https://stackoverflow.com/questions/32581950/removing-zeros-in-a-list-of-vectors
    if (length(neighbors) == 0) {
      new = 0
    }
    else {
      chosen = sample(neighbors, 1)
      new = rbernoulli(1, p = b[chosen]) * chosen
    }
    return(new)
  }
  
  
  

  #landscape_dt = rep(0 , dim^2)
  landscape_dt = sample(1:5, dim^2, replace=T)
  landscape = matrix(landscape_dt, nrow = dim,ncol = dim, byrow=TRUE)
  
  
  
  maxtime = runtime
  
  diversity = rep(0 , maxtime)
  
  harvest_count = rep(0, 5)
  
  for (t in 1:maxtime) {
    
    #natural death
    for (i in 1:dim) { 
      for (j in 1:dim) {
        if (landscape[i,j] != 0) {
          spec = landscape[i,j]
          landscape[i,j] = death(spec)
        }
      }
    }
    
    #harvesting: based on relative abundances
    counts = rep(0,6)
    for (i in 1:dim) {
      for (j in 1:dim) {
        spec = landscape[i,j]
        counts[spec] = counts[spec] + 1
      }
    }
    counts = counts[-1]
    counts = counts / sum(counts)
    harvest_probs = counts
    #harvest_probs = rep(.2, 5)
    for (i in 1:dim) {
      for (j in 1:dim) {
        if (landscape[i,j] != 0) {
          spec = landscape[i,j]
          trial = rbernoulli(1, p = 1 - harvest_probs[spec])
          if (trial == 0){
            harvest_count[spec] = harvest_count[spec] + 1
          }
          landscape[i,j] = trial * spec
        }
      }
    }
    
    # points = sort(ceiling(runif(4) * dim))
    # 
    # for (i in points[1]: points[2]) {
    #   for (j in  points[3]: points[4]) {
    #     if (landscape[i,j] != 0) {
    #       spec = landscape[i,j]
    #       trial = rbernoulli(1, p = .5)
    #       if (trial == 0){
    #         harvest_count[spec] = harvest_count[spec] + 1
    #       }
    #       landscape[i,j] = trial * spec
    #     }
    #   }
    # }
    
    
    prev = landscape
    
    #corner births
    if (prev[1,1] == 0) {
      landscape[1,1] = birth(c(prev[1,2], prev[2,1], prev[2,2]))
    }
    if (prev[dim,1] == 0) {
      landscape[dim,1] = birth(c(prev[dim-1,1], prev[dim,2], prev[dim-1,2]))
    }
    if (prev[1,dim] == 0) {
      landscape[1,dim] = birth(c(prev[1,dim-1], prev[2,dim], prev[2,dim-1]))
    }
    if (prev[dim,dim] == 0) {
      landscape[dim,dim] = birth(c(prev[dim-1,dim], prev[dim,dim-1], prev[dim-1,dim-1]))
    }
    
    #non corner edges
    for (i in 1:1) { 
      for (j in 2:dim-1) {
        if (prev[i,j] == 0) {
          neighbors = c(prev[i+1,j+1], prev[i+1,j-1], prev[i+1,j], prev[i,j+1], prev[i,j-1])
          landscape[i,j] = birth(neighbors)
        }
      }
    }
    
    for (i in dim:dim) { 
      for (j in 2:dim-1) {
        if (prev[i,j] == 0) {
          neighbors = c(prev[i-1,j-1], prev[i-1,j+1], prev[i,j+1], prev[i,j-1], prev[i-1,j])
          landscape[i,j] = birth(neighbors)
        }
      }
    }
    
    for (i in 2:dim-1) { 
      for (j in 1:1) {
        if (prev[i,j] == 0) {
          neighbors = c(prev[i+1,j+1], prev[i-1,j+1], prev[i+1,j], prev[i,j+1], prev[i-1,j])
          landscape[i,j] = birth(neighbors)
        }
      }
    }
    
    for (i in 2:dim-1) { 
      for (j in dim:dim) {
        if (prev[i,j] == 0) {
          neighbors = c(prev[i-1,j-1], prev[i+1,j-1], prev[i+1,j], prev[i,j-1], prev[i-1,j])
          landscape[i,j] = birth(neighbors)
        }
      }
    }
    
    #interior points
    for (i in 2:dim-1) { 
      for (j in 2:dim-1) {
        if (prev[i,j] == 0) {
          neighbors = c(prev[i+1,j+1], prev[i-1,j-1], prev[i+1,j-1], prev[i-1,j+1], prev[i+1,j], prev[i,j+1], prev[i,j-1], prev[i-1,j])
          landscape[i,j] = birth(neighbors)
        }
      }
    }
    counts = table(landscape)
    counts = counts[-1]
    counts = counts / sum(counts)
    diversity[t] = length(counts[counts > .1])
  }
  params=c(size_a, size_b, size_c, size_d, size_e, b_factor, d_prop,
           runtime, dim)
  names(params)<-c("size_a", "size_b", "size_c", "size_d", 
                   "size_e", "b_factor_1", "b_factor_2", "b_factor_3", "b_factor_4"
                   , "b_factor_5", "d_prop",
                   "runtime", "dim")
  return(list(outcome = diversity, parameters = params
  ))
}

