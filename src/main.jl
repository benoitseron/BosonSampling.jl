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

include("constructors.jl")
include("hpsd.jl")
include("io.jl")
include("matrix_tests.jl")
include("types.jl")
include("tools.jl")
include("output_statistics.jl")
include("permanents.jl")
include("proba_tools.jl")
include("scattering.jl")
include("sampling_func.jl")
include("special_matrices.jl")
include("tools.jl")
include("type_function.jl")
include("methods.jl")
