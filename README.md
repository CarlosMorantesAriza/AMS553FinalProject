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
| $S$ | Size | Arbitrary selected sizes. Large size = low death rate, low growth rate.  |


## Growth rates

$$r_{i_{t}}$$


**Notes on the order of the order of the simulation**

## Initial conditions:

### Initial conditions options.
* All trees start with the same abundance
* Tree species start with random abundances.
* Tree species start with abundances inversely proportional to their growth rate means.


# Non-harvesting model

## In each iteration

1. Draw number of individuals per species that die in that cycle.
   * **Model options**
     * Draw 'die/not-die' state for each individual ( $d \sim Bernoulli(p_{i})$ where $p_{i}$ is the probability of an individual of species $i$ dying from causes different than harvesting). Computationally inefficient. Boring :) .
     * Draw 'die/not-die' state for the current population of species $i$ ( $d_{t} \sim Binom(N_{i_{t}},p_{i})$ ). Faster, little less boring.
2. Generate the number of individuals per species that enter the population.
   * **Model options**
     * Two steps approach:
      * Step 1: Generate number of saplings produced per tree per species $Pois(\frac{size_{i}}{size_{max}})$ .
      * Step 2: Draw 'recruit/not-recruit' state for each sapling (maybe this should be something like $Bernoulli(\frac{\sum_i N_i}{K})$ ). Or draw 'recruit/not-recruit' state for all saplings from a single species (maybe this should be something like $Binom(saplings_{i}, \frac{saplings_{i}}{\sum_{i} saplings_{i}})$ )
     *  Single step approach: Generate number of recruited trees per species (maybe something like $Pois(\frac{\sum_i N_i}{K})$ or $Pois(\frac{K- \sum_i N_i}{K} \frac{N_i}{\sum_i N_i})$ )   
     
3) Measure the number of individuals per species.

# Harvesting model


1. Draw number of individuals per species that die in that cycle.
   * **Model options**
     * Draw 'die/not-die' state for each individual ( $d \sim Bernoulli(p_{i})$ where $p_{i}$ is the probability of an individual of species $i$ dying from causes different than harvesting). Computationally inefficient. Boring :) .
     * Draw 'die/not-die' state for the current population of species $i$ ( $d_{t} \sim Binom(N_{i_{t}},p_{i})$ ). Faster, little less boring.
2. Generate the number of individuals per species that enter the population.
   * **Model options**
     * Two steps approach:
      * Step 1: Generate number of saplings produced per tree per species.
      * Step 2: Draw 'recruit/not-recruit' state for each sapling (maybe this should be something like $Bernoulli(\frac{\sum_i N_i}{K})$ ). Or draw 'recruit/not-recruit' state for all saplings from a single species (maybe this should be something like $Binom(saplings_{i}, \frac{saplings_{i}}{\sum_{i} saplings_{i}})$ )
     *  Single step approach: Generate number of recruited trees per species (maybe something like $Pois(\frac{\sum_i N_i}{K})$ or $Pois(\frac{K- \sum_i N_i}{K} \frac{N_i}{\sum_i N_i})$ )   
     
3. Select trees to harvest. 
Each species has its own harvesting rate (long lived low, short-lived high)
Harvest by proportion:
  * What is the effect of harvesting y percent of each species every z years.


| Species | Proportion  | Time step |
| ------|:------ | ----------- |
| A | $j \in [0.01, 0.05, 0.1 ] $ | $t \in [ 1,5,10 ] $ |
| B | $j \in [ 0.01, 0.05, 0.1 ] $ | $t \in [ 1,5,10 ] $ |
| C | $j \in [ 0.01, 0.05, 0.1 ] $ | $t \in [ 1,5,10 ] $ |
| D | $j \in [ 0.01, 0.05, 0.1 ] $ | $t \in [ 1,5,10 ] $ |
| E | $j \in [ 0.01, 0.05, 0.1 ] $ | $t \in [ 1,5,10 ] $ |



  * **Model options**
    * Spatially explicit model: 




Scenario Spreadsheet: https://docs.google.com/spreadsheets/d/19jHNMtIbE_P-EOxpflBBQBJnVxJdUdxf1SgPhaa7HUs/edit?usp=sharing



# Relative abundances

$$r_i= \sum_{i}^{N} frac{n_i}{N}$$
