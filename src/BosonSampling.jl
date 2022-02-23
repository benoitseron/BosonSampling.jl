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
using Parameters
using ArgCheck
using Distributions

const ATOL = 1e-10

include("special_matrices.jl")
include("matrix_tests.jl")
include("proba_tools.jl")
include("type_functions.jl")
include("types.jl")
include("scattering.jl")

include("circuits/circuit_elements.jl")

include("bunching/bunching.jl")
include("partitions/legacy.jl")

include("permanent_conjectures/bapat_sunder.jl")
include("permanent_conjectures/counter_example_functions.jl")
include("permanent_conjectures/counter_example_numerical_search.jl")
include("permanent_conjectures/permanent_on_top.jl")

permanent = ryser

### temporary permanent before update of Julia rep ###

function is_a_square_three_tensor(W::Array)

	length(size(W)) == 3 && all(size(W) .== size(W)[1])

end


function naive_tensor(W::Array)

    if is_a_square_three_tensor(W)
        n = size(W)[1]

        sum_tensor_diag(sigma,rho, W) = prod(W[i,sigma[i],rho[i]] for i = 1:n)
        return sum(sum_tensor_diag(sigma,rho, W) for sigma in permutations(collect(1:n)) for rho in permutations(collect(1:n)))
    else
        error("tensor permanent implemented only for square 3-indices tensors")
    end

end

function ryser_tensor(W::Array)

	"""implements the tensor permanent according to Ryser's factorization
	without using Gray factoring

	corresponds to Eq. 21 of https://arxiv.org/abs/1410.7687"""

	@argcheck is_a_square_three_tensor(W) "tensor permanent implemented only for square 3-indices tensors"

	n = size(U)[1]
    nstring = collect(1:n)
    sub_nstring = collect(powerset(nstring))
    sub_nstring = sub_nstring[2:length(sub_nstring)]
    res = 0

    function delta_set(S1, S2)
        if S1 == S2
            return 1
        else
            return 0
        end
    end

    for r = 1:length(sub_nstring)
        for s = r:length(sub_nstring)
            R = sub_nstring[r]
            S = sub_nstring[s]

            t = prod(
                sum(
					W[rr,ss,j] for rr in R for ss in S
                ) for j = 1:n
            )
            res +=
                (2 - delta_set(S, R)) * (-1)^(length(S) + length(R)) * real(t)
        end
    end
    return res

end

for n in names(@__MODULE__; all=true)
    if Base.isidentifier(n) && n ∉ (Symbol(@__MODULE__), :eval, :include)
        @eval export $n
    end
end

end
