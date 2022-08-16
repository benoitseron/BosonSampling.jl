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


function foo(n, r_s, r_z)

    U = copy(rand_haar(n))

    S = rand_gram_matrix_rank(n,r_s)

    @argcheck rank(S) == r_s

    phases = zeros(Complex, n)

    phases[1:r_z] = exp.(1im * rand(r_z)) .- 1

    diagm(phases)

    pr = U' * diagm(phases) * U

    pr = convert(Matrix{ComplexF64}, pr)

    @show rank(pr .* S)

end


for i in 1:1 
    foo(10,2,2)
end
