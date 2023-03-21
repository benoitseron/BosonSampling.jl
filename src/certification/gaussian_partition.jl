begin
    using Revise

    using BosonSampling
    using Plots
    using ProgressMeter
    using Distributions
    using Random
    using Test
    using ArgCheck
    using StatsBase
    using ColorSchemes
    using Interpolations
    using Dierckx
    using LinearAlgebra
    using PrettyTables
    using LaTeXStrings
    using JLD
    using AutoHashEquals
    using LinearRegression
    using DataStructures
    using Parameters
    using UnPack
    using Optim


    using Permanents
    using Plots
    using Test
    using Combinatorics
    using Random
    using IterTools
    using Statistics
    using LinearAlgebra #removed so as to be able to use generic types such as BigFloats, can be put back if needed
    using PolynomialRoots
    using StatsBase
    #using JLD
    using CSV
    using DataFrames
    using Tables
    using Plots#; plotly() #plotly allows to make "dynamic" plots where you can point the mouse and see the values, they can also be saved like that but note that this makes them heavy (and html instead of png)
    using PrettyTables #to display large matrices
    using Roots
    using BenchmarkTools
    using Optim
    using ProgressMeter
    using ProgressBars
    using Parameters
    using ArgCheck
    using Distributions
    using Luxor
    using AutoHashEquals
    using LinearRegression
    using HypothesisTests
    using SimpleTraits
    using Parameters
    using UnPack
    using Dates
    using JLD
    using DelimitedFiles
    using BosonSampling
    using Plots
    using ProgressMeter
    using Distributions
    using Random
    using Test
    using ArgCheck
    using StatsBase
    using ColorSchemes
    using Interpolations
    using Dierckx
    using LinearAlgebra
    using PrettyTables
    using LaTeXStrings
    using JLD
    using AutoHashEquals
    using LinearRegression
    using DataStructures
    using Parameters
    using UnPack
    using Combinatorics
    using BosonSampling
    using Plots
    using ProgressMeter
    using Distributions
    using Random
    using Test
    using ArgCheck
    using ColorSchemes
    using Interpolations
    using Dierckx
    using LinearAlgebra
    using PrettyTables
    using LaTeXStrings
    using JLD
    using AutoHashEquals
    using LinearRegression
    using DataStructures
    using Parameters
    using UnPack
    using Optim
    using DelimitedFiles
    using Dates
    using DelimitedFiles    

    using ITensors
end


# r - squeezing
# x - generatinc funtion
# η - vector of fourier coefficients
# m - modes
# Q - 
# Λ -
# δ -
# λ - thermal parameters

# the lambdas are created page 24
# Q page 23
# C just above

m = 4
r = 1.5 * ones(m)
λ = 1 * ones(m)

C_array = [0.5 - 1/(1+ λ[j] * exp(2*r[j])) for j in 1:m]

C = diagm(C_array)