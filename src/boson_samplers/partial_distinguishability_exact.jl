# obtain the S matrix
# make the cholesky decomposition
# use it to separate modes into indistinguishable parts between photons
# interfere each indistinguishable part in its own interferometer (so up to n interferometers)
# recombine 

using Revise
using BosonSampling
using LinearAlgebra

n = 3
m = n
r = 2

interf = RandHaar(m)
TIn = GramMatrix


input_state = Input{TIn}(first_modes(n,m))
i = input_state


S = rand_gram_matrix(n)
GramMatrix{ComplexF64}(n, )


    function Input{T}(r::ModeOccupation, n::Int, m::Int, G::GramMatrix, distinguishability_param::Union{Real,Nothing}) where {T<:InputType}
        new{T}(r,n,m,G, distinguishability_param)
    end


S = rand_gram_matrix_from_orthonormal_basis(n,r)
