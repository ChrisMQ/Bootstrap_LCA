##-------------------------------------------------------------------------------------#
## What Can Bootstrap Do in Latent Class Analysis? A Whole-Process Bootstrap Framework
## Author: Meng Qiu
## Last update: February 2024
##-------------------------------------------------------------------------------------#
## Empirical example: Stouffer-Toby Data
##-------------------------------------------------------------------------------------#

## Load libraries
library(poLCA) 
library(dplyr)
library(descr)
library(combinat)

## Load functions
source("functions.R")

## Load Stouffer-Toby data
data(values)

## Part 1: Model selection & bootstrap selection rates
# Arguments:
# data: Observed data
# nrep: Number of times to estimate the model for the original data
# nrep.bt: Number of times to estimate the model for a bootstrapped data
# maxiter: The maximum number of iterations through which the estimation algorithm will cycle
# K: The maximum number of classes in the class enumeration process
# B: Number of bootstrap samples
# conf.level: Confidence level
# verbose: Logical, indicating whether LCA.bt.NB should show a progress bar
set.seed(123)
part1 <- LCA.bt.NB(data=values, nrep=20, nrep.bt=1, maxiter=10000, 
                   K=5, B=500, conf.level=0.95, verbose=T)
# Table 4
part1$select.single # number of classes chosen by each index
part1$indx.value    # values of indices for k varying from 1 to K 
# Table 5
part1$select.rate

## Part 2: Detect local dependence
# Arguments:
# data: Observed data
# nclass: The selected number of classes
part2 <- bvr.pVal(data=values, nclass=2)
# Table 6
part2

## Part 3: Parameter estimation & classification uncertainty
# Arguments:
# MS.out: The output from the LCA.bt.NB function
# nclass: The selected number of classes
# conf.level: Confidence level
part3 <- LCA.bt.EST(MS.out=part1, nclass=2, conf.level=.95)
# Table 7
part3$Estimates
# Table 8
part3$CI
# Table 9
part3$Entropy
