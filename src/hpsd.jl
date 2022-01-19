# hermitian positive semi definite matrices
# follows https://arxiv.org/abs/1609.02416

include("special_matrices.jl")
include("matrix_tests.jl")
include("scattering.jl")



function sample_H_matrix(n::Int,m::Int)

    """samples a unitary of size m according to the Haar measure and returns the
    H matrix corresponding to n photons sent going from the first n top input
    modes to the top n output modes"""

    U = rand_haar(m)

    input_state = [i <= n ? 1 : 0 for i in 1:m]
    output_state = input_state # variations of this are implicit in the Haar random sampling

    H_matrix(U, input_state, output_state)

end

### example ###

@test is_positive_semidefinite(sample_H_matrix(5,4))
