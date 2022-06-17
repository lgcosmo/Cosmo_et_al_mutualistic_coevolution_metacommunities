using DrWatson
using CSV
using Distributions
using DataFrames
using Statistics
using LinearAlgebra
using Random
using RCall
using Plots
using StatsPlots

#Loading functions

include(srcdir("main_functions.jl"))
include(srcdir("aux_functions.jl"))

#Creating metacommunity

G=moore_neighborhood(n_patches=100, n=10, periodic=true) #Creating transition matrix for dispersal

#Simulation

#The function coevo_metacom take as input several parameters and return a dataframe containing the variables used in the main text
#The model parameters are:
#n_sp - Initial species richness of the metacommunity
#climchange  - Amount of directional climate change at each time step
#mi - Strength of mutualistic interactions
#α - Sensitivity to trait matching of mutualistic coevolution
#ρ - Sensitivity of species adaptive landscape/suitability to changes in trait values
#σ - Additive genetic variance
#flow - Total fraction of gene flow that each patch can receive from neighbours
#tmax - Maximum simulation time
#sim - Variable to identify the simulation

r=coevo_metacom(n_sp=16, G=G, climchange=0.0, mi=0.5, α=0.1, ρ=0.1, σ=1.0, flow=0.05, tmax=500, sim=1)

#Plotting occupancy over time for each species

@df r plot(:time, :occupancy, group=:species)