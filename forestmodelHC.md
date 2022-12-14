---
title: "forest stuff"
author: "h r c"
date: "2022-11-23"
output: html_document
---

```{r}
library(purrr)
```


```{r cars}
set.seed(1)

sizes = c(10, 10, 10, 10, 10)
b = 2 / sizes #birth chance
d = 1 / sizes #death chance
```

```{r}
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

harvest <- function(landscape) {
   return()
}
```


```{r pressure, echo=FALSE}
dim = 100
#landscape_dt = rep(0 , dim^2)
landscape_dt = sample(1:5, dim^2, replace=T)
landscape = matrix(landscape_dt, nrow = dim,ncol = dim, byrow=TRUE)
```

```{r}
maxtime = 2000

diversity = rep(0 , maxtime)

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
  counts = table(landscape)
  counts = counts[-1]
  counts = counts / sum(counts)
  harvest_probs = counts
  for (i in 1:dim) {
    for (j in 1:dim) {
      if (landscape[i,j] != 0) {
        spec = landscape[i,j]
        landscape[i,j] = rbernoulli(1, p = 1 - harvest_probs[spec]) * spec
      }
    }
  }
  
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
```

```{r}
plot(diversity)
```

```{r}
counts = table(landscape)
counts[-1]
```
