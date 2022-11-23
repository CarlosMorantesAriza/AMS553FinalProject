---
title: "forest stuff"
author: "h r c"
date: "2022-11-23"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

```{r}
library(purrr)
```


```{r cars}
sizes = c(2, 7, 9, 10, 11)
b = 1 / sizes
d = 2 / sizes
```


```{r pressure, echo=FALSE}
dim = 100
landscape_dt <- rep(0 , dim^2)
landscape <- matrix(landscape_dt, nrow = dim,ncol = dim, byrow=TRUE)
#death
for (i in 1:dim) { 
  for (j in 1:dim) {
    if (landscape[i,j] != 0) {
      spec = landscape[i,j]
      landscape[i,j] = rbernoulli(1, p = d[spec]) * 1
    }
  }
}
prev = landscape
#birth
for (i in 2:dim-1) { 
  for (j in 2:dim-1) {
    if (prev[i,j] == 0) {
      neighbors = c(prev[i,j], prev[i+1,j+1], prev[i-1,j-1], prev[i+1,j-1], prev[i-1,j+1], prev[i+1,j], prev[i,j+1], prev,prev[i,j-1])
    }
  }
}
```
