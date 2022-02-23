module BosonSampling

using Permanents
using Plots
using Test
using Combinatorics
using Random
using IterTools
using Statistics
using LinearAlgebra
using PolynomialRoots
using StatsBase
using JLD
using CSV
using DataFrames
using Tables
using PrettyTables
using Roots
using BenchmarkTools
using Optim
using ProgressMeter
using LaTeXStrings

const ATOL = 1e-10

include("special_matrices.jl")
include("matrix_tests.jl")
include("proba_tools.jl")
include("tools.jl")
include("circuits/circuit_elements.jl")

include("type_functions.jl")
include("types.jl")
include("scattering.jl")
include("bunching/bunching.jl")
include("partitions/legacy.jl")

include("permanent_conjectures/bapat_sunder.jl")
include("permanent_conjectures/counter_example_functions.jl")
include("permanent_conjectures/counter_example_numerical_search.jl")
include("permanent_conjectures/permanent_on_top.jl")

include("boson_samplers/classical_sampler.jl")
include("boson_samplers/cliffords_sampler.jl")
include("boson_samplers/metropolis_sampler.jl")
include("boson_samplers/noisy_sampler.jl")

include("distributions/noisy_distribution.jl")
include("distributions/theoretical_distribution.jl")

permanent = ryser

for n in names(@__MODULE__; all=true)
    if Base.isidentifier(n) && n ∉ (Symbol(@__MODULE__), :eval, :include)
        @eval export $n
    end
end

end
