# functions relating to more general input to the generalized bunching conjecture
# instead of pure state inputs

"""

    schur_matrix(H)

computes the Schur matrix as defined in Eq. 1 of Linear Algebra and its Applications 490 (2016) 196–201
"""
function schur_matrix(H)

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

function find_permutation_index(this_perm, permutation_array)
    findfirst(x->x == this_perm, permutation_array)
end

function multiply_permutations(a,b)
    a[b[:]]
end

"""

    J_array(theta, n)

returns the `J` as defined in in Eq.10 of [Universality of Generalized Bunching and Efficient Assessment of Boson Sampling](https://arxiv.org/abs/1509.01561), with the permutations coming in the order given by `permutations(collect(1:n))`
"""
function J_array(theta, n)



    J = zeros(eltype(theta), factorial(n))

    permutation_array = collect(permutations(collect(1:n)))

    @showprogress for (i,sigma) in enumerate(permutations(collect(1:n)))
        for (j,tau) in enumerate(permutations(collect(1:n)))
            J[i] += conj(theta[j]) * theta[find_permutation_index(multiply_permutations(tau, sigma), permutation_array)]
        end
    end

    J

end


"""
    density_matrix_from_J(J,n)

density matrix associated to a J function as defined in [Universality of Generalized Bunching and Efficient Assessment of Boson Sampling](https://arxiv.org/abs/1509.01561), computed through Eq. 46.
"""
function density_matrix_from_J(J,n)

    density_matrix = zeros(eltype(J), factorial(n), factorial(n))
    permutation_array = collect(permutations(collect(1:n)))

    @showprogress for (i,sigma) in enumerate(permutations(collect(1:n)))
        for (j,tau) in enumerate(permutations(collect(1:n)))
            density_matrix[j,find_permutation_index(multiply_permutations(sigma, tau), permutation_array)] += J[i]
        end
    end

    density_matrix/factorial(n)

end
