module BosonSampling

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

# using Distributed
# using SharedArrays

@consts begin
    ATOL = 1e-10
    SAFETY_FACTOR_FULL_BUNCHING = 10
end

include("special_matrices.jl")
include("matrix_tests.jl")
include("proba_tools.jl")
include("circuits/circuit_elements.jl")
include("types/type_functions.jl")
include("types/types.jl")
include("scattering.jl")

include("bunching/bunching.jl")
include("partitions/legacy.jl")
include("partitions/partition_expectation_values.jl")
include("partitions/partitions.jl")


include("boson_samplers/tools.jl")
include("boson_samplers/classical_sampler.jl")
include("boson_samplers/cliffords_sampler.jl")
include("boson_samplers/methods.jl")
include("boson_samplers/metropolis_sampler.jl")
include("boson_samplers/noisy_sampler.jl")
include("boson_samplers/sample.jl")
include("boson_samplers/gaussian_sampler.jl")

include("distributions/noisy_distribution.jl")
include("distributions/theoretical_distribution.jl")
include("distributions/full_distribution.jl")

include("permanent_conjectures/bapat_sunder.jl")
include("permanent_conjectures/counter_example_functions.jl")
include("permanent_conjectures/counter_example_numerical_search.jl")
include("permanent_conjectures/permanent_on_top.jl")

include("certification/experimental_data_generation.jl")
include("certification/bayesian.jl")
include("certification/correlators.jl")

include("visual.jl")

include("loop/loop_functions.jl")

include("experiments/data_conversion.jl")
include("experiments/coincidences.jl")

include("distinguishability/parametrizations.jl")



permanent = ryser

for n in names(@__MODULE__; all=true)
    if Base.isidentifier(n) && n ∉ (Symbol(@__MODULE__), :eval, :include)
        @eval export $n
    end
end

end
