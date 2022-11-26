library(purrr)
library(tidyverse)
library(MASS)
library(abind)
library(future)
library(furrr)
library(parallel)
library(tictoc)
library(future.apply)
trabajadores=detectCores()-8
source('henryscode.R')
parameter_tibble=
  expand.grid(
    size_a=c(1,5,10,20), 
    size_b=c(2),
    size_c=c(3),
    size_d=c(4),
    size_e=c(5),
    d_prop=seq(0.1,0.5,0.1),
    runtime=100,
    dim=100
    )


tic('I finished running')
plan(multisession, workers = trabajadores)
results = 
  future_pmap(
    .l = parameter_tibble,
    .f = henry_model,
    .options=furrr_options(seed = TRUE)
  )



toc()
df=as_tibble(bind_rows(lapply(results, results_sorting)))
