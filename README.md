# AMS553FinalProject

Forests are an important source of biodiversity and economic resources globally. An important question is to determine how to harvest trees sustainably while also conserving biodiversity. Forest systems are difficult to study empirically because they operate on large spatial and temporal scales. Therefore, to address this question, simulation modeling is an ideal tool. For this project, we propose two simulate a simplified forest with 5 species without any harvesting, and then calculate biodiversity metrics. Then, we propose to simulate the same forest testing different proportions and sizes of trees harvested, and calculate the same biodiversity metrics. This will allow us to make sustainable logging recommendations. 

## The model

$$N_{TOT_{t+1}}= \sum_{i} r_{i} N_{i_{t}} \left(1-\frac{\sum_{i} N_{i_{t}}}{K}\right) - \sum_{i}h_{i}N_{i_{t}}$$


### Parameter

| Symbol | Name  | Notes       |
| ------|:------ | ----------- |
|$r_{i_{t}}$| Growth rate | One per species per time |
| $h_{i_{t}}$      | Harvest rate       |   One per species per time          |
| $K$  | Total forest size       |  Single value |


## Growth rates

$$r_{i_{t}}$$


**Notes on the order of the order of the simulation**

## Initial conditions:

### Initial conditions options.
* All trees start with the same abundance
* Tree species start with random abundances.
* Tree species start with abundances inversely proportional to their growth rate means.


# Non-harvesting model

1. Draw number of individuals per species that die in that cycle.
   * **Model options**
     * Draw 'die/not-die' state for each individual ( $d \sim Bernoulli(p_{i})$ where $p_{i}$ is the probability of an individual of species $i$ dying from causes different than harvesting). Computationally inefficient. Boring :) .
     * Draw 'die/not-die' state for the current population of species $i$ ( $d_{t} \sim Binom(N_{i_{t}},p_{i})$ ). Faster, little less boring.
2. Generate the number of individuals per species that enter the population.
   * **Model options**
     * Two steps approach:
      * Step 1: Generate number of saplings per tree per species.
      * Step 2: Draw 'recruit/not-recruit' state for each sapling (I assume this should be something like $Binom(\frac{\sum_i N_i}{K})$).
     

## In each iteration

1) Generate random number of individuals per species that die in the iteration.
2) Generate the number of individuals per species that enter the population.










