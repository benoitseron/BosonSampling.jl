using ProgressMeter
include("counter_example_circuit.jl")

function schur_matrix(H)

    """computes the Schur matrix as defined in Eq. 1 of Linear Algebra and its Applications 490 (2016) 196–201"""

    @test is_positive_semidefinite(H)

    n = size(H,1)
    schur_mat = Matrix{eltype(H)}(undef,(factorial(n), factorial(n)))

    @showprogress for (i,sigma) in enumerate(permutations(collect(1:n)))
        for (j,rho) in enumerate(permutations(collect(1:n)))
            schur_mat[i,j] = prod([H[sigma[k], rho[k]] for k in 1:n])
        end
    end

    schur_mat

end

schur_drury = schur_matrix(H)

# spectrum_schur_drury = eigvals(schur_drury)
#
# maximum(spectrum_schur_drury) > abs(permanent_ryser(H))
#
# maximum(spectrum_schur_drury) / abs(permanent_ryser(H))

eigen_dec = eigen(schur_drury)

theta = eigen_dec.vectors[:, findmax(eigen_dec.values)[2]]
# eigenvector corresponding to the greatest eigenvalue

@test theta' * schur_drury * theta ≈ findmax(eigen_dec.values)[1]
# health check


function find_permutation_index(this_perm, permutation_array)
    findfirst(x->x == this_perm, permutation_array)
end

function multiply_permutations(a,b)

    a[b[:]]
end

function J_array(theta, n)

    """
    returns the J as defined in in Eq.10 of https://arxiv.org/abs/1509.01561, with the permutations coming in the order given by permutations(collect(1:n))
    """

    J = zeros(eltype(theta), factorial(n))

    permutation_array = collect(permutations(collect(1:n)))

    @showprogress for (i,sigma) in enumerate(permutations(collect(1:n)))
        for (j,tau) in enumerate(permutations(collect(1:n)))
            J[i] += conj(theta[j]) * theta[find_permutation_index(multiply_permutations(tau, sigma), permutation_array)]
        end
    end

    J

end

J_drury = J_array(theta, 7)
@test J_drury[1] ≈ 1

function density_matrix_from_J(J,n)

    """density matrix associated to a J function as defined in https://arxiv.org/abs/1509.01561, computed through Eq. 46"""

    density_matrix = zeros(eltype(J), factorial(n), factorial(n))
    permutation_array = collect(permutations(collect(1:n)))

    @showprogress for (i,sigma) in enumerate(permutations(collect(1:n)))
        for (j,tau) in enumerate(permutations(collect(1:n)))
            density_matrix[j,find_permutation_index(multiply_permutations(sigma, tau), permutation_array)] += J[i]
        end
    end

    density_matrix/factorial(n)

end

rho_drury = density_matrix_from_J(J_drury, 7)

@test abs(tr(rho_drury)) ≈ 1. atol = 1e-10







permutation_array = collect(permutations(collect(1:7)))

CSV.write("counter_examples/mixed_state_input_drury/permutation_array.csv",  DataFrame(transpose(permutation_array),:auto), header=false)
CSV.write("counter_examples/mixed_state_input_drury/schur_drury.csv",  DataFrame(schur_drury,:auto), header=false)
CSV.write("counter_examples/mixed_state_input_drury/J_drury.csv",  DataFrame(transpose(J_drury),:auto), header=false)
CSV.write("counter_examples/mixed_state_input_drury/rho_drury.csv",  DataFrame(rho_drury,:auto), header=false)

J_drury
