using Plots
using Test
using Combinatorics
using Random
using IterTools
using Statistics
using Distributions
using LinearAlgebra #removed so as to be able to use generic types such as BigFloats, can be put back if needed
using GenericLinearAlgebra
using PolynomialRoots
using FFTW
using StatsBase
using JLD
using CSV
using DataFrames
using Tables
using Plots#; plotly() #plotly allows to make "dynamic" plots where you can point the mouse and see the values, they can also be saved like that but note that this makes them heavy (and html instead of png)
using PrettyTables #to display large matrices
using Roots
using BenchmarkTools
using Optim
using ProgressMeter

const ATOL = 1e-10

include("src/constructors.jl")
include("src/hpsd.jl")
include("src/io.jl")
include("src/types.jl")
include("src/tools.jl")
include("src/permanents.jl")
include("src/proba_tools.jl")
include("src/scattering.jl")
include("src/special_matrices.jl")
include("src/tools.jl")
include("src/type_function.jl")
include("src/distributions/theoretical_distribution.jl")
include("src/distributions/noisy_distribution.jl")
include("src/boson_samplers/classical_sampler.jl")
include("src/boson_samplers/metropolis_sampler.jl")
include("src/boson_samplers/noisy_sampler.jl")
include("src/boson_samplers/cliffords_sampler.jl")
include("src/methods.jl")
