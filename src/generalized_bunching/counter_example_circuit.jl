include("counter_example_functions.jl")

omega_5 = exp(2/5 * 1im * pi)
M_dagger_unnormalized = 1/sqrt(2) .* [[sqrt(2) 0]; [0 sqrt(2)]; [1 1]; [1 omega_5]; [1 omega_5^2]; [1 omega_5^3]; [1 omega_5^4]]

H_circuit_unnormalized = M_dagger_unnormalized * M_dagger_unnormalized' # unnormalized meaning that it is a gram matrix though the normalization of H must come from U
S_circuit = copy(conj(H_circuit_unnormalized'))

### starting from the M expression ###

M_dagger = sqrt(2/7) .* M_dagger_unnormalized

M_dagger' * M_dagger

H = M_dagger * M_dagger'
H_circuit = H

S = copy(conj((M_dagger_unnormalized * M_dagger_unnormalized')))

@test violates_bapat_sunder(H,S)

U = add_columns_to_make_square_unitary(M_dagger)

V = U'

# the scattering matrix is manifestly U' and not U
