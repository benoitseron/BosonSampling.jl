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

    n = 5
    r_s = 4
    r_z = 3

    U = copy(rand_haar(n))

    S = rand_gram_matrix_rank(n,r_s)

    @argcheck rank(S) == r_s

    phases = zeros(Complex, n)

    phases[1:r_z] = exp.(1im * rand(r_z)) .- 1

    diagm(phases)

    pr = U' * diagm(phases) * U

    pr = convert(Matrix{ComplexF64}, pr)

    @argcheck rank(pr) == r_z "doesn't hold"

end


for i in 1:10000
    foo(10,4,2)
end
