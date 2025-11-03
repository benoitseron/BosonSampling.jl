"""
    duplicate_row_col(A::Array, lengths::Vector{Int})

Returns a matrix where the column and row ``i`` of ``A`` are repeated ``lengths_i`` times.
"""
function duplicate_row_col(A::Array, lengths::Vector{Int})

    x = axes(A,1)
    @argcheck length(x) == length(lengths)

    res = similar(x, sum(lengths))
    i = 1

    for idx in 1:length(x)
        tmp = x[idx]
        for kdx in 1:lengths[idx]
            res[i] = tmp
            i += 1
        end
    end

    if length(size(A)) == 2
        return A[res,res]
    else
        return A[res]
    end

end

"""
    williamson_decomp(V::AbstractMatrix)

Find the Williamson decomposition of the positive semidefinite matrix ``V``.
Returns ``D`` and ``S`` such that ``V = S D S^T
!!! note "Reference"
    [The Walrus documentation](https://the-walrus.readthedocs.io/en/latest/code/decompositions.html?highlight=williamson#thewalrus.decompositions.williamson)
"""
function williamson_decomp(V::AbstractMatrix)
    @argcheck isposdef(V) && size(V)[1] == size(V)[2] && size(V)[1] % 2 == 0

    n = div(LinearAlgebra.checksquare(V), 2)
    J = symplectic_mat(n)
    v = real.(sqrt(inv(V)))
    inter = v * J * v

    F = schur(inter)
    T = F.T
    Z = F.Z

    perm = vcat([i for i in 1:2:2n-1], [i+1 for i in 1:2:2n-1])
    x = [0 1; 1 0]
    id = Matrix{eltype(V)}(I, 2, 2)
    seq = []
    for i in 1:2:2n-1
        T[i,i+1] > 0 ? push!(seq, id) : push!(seq, x)
    end

    p = seq[1]
    for i in 2:length(seq)
        p = direct_sum(p, seq[i])
    end

    Zp = Z * p
    Zp = Zp[:, perm]
    Tp = p * T * p
    d = [1/Tp[i,i+1] for i in 1:2:2n-1]
    D = diagm(vcat(d,d))

    S = transpose(inv(v * Zp * sqrt(D)))
    return D, S
end

function LaplaceExpansion(perm_mat, full_mat)

    s_full = size(full_mat)[1]
    s_perm = size(perm_mat)[2]
    res = Vector{Float64}(undef, s_full)
    global v_perms = Vector{ComplexF64}(undef, s_perm)

    Threads.@threads for i in 1:s_perm
        v_perms[i] = ryser(perm_mat[:,1:end .!= i])
    end

    @simd for i in 1:s_full
        rowA = full_mat[i,:]
        res[i] = abs.(dot(rowA, v_perms')).^2
    end

    return res

end

function total_variation_distance(p, q)

    if length(p) > length(q)
        vcat(q, zeros(length(p)-length(q)))
    elseif length(q) > length(p)
        vcat(p, zeros(length(q)-length(p)))
    end

    return sum(abs(p[i]-q[i]) for i = 1:length(p))

end

function collect_sub_mat_perm(U)
    sub_mat = [remove_row_col(U, [], [i]) for i in 1:size(U)[2]]
    return [fast_glynn_perm(m) for m in sub_mat]
end

"""
    gram_to_coefficients(S::Matrix; atol=1e-10)

Extract internal degrees of freedom coefficients from a Gram matrix.

Given a Gram matrix S where S[i,j] = ⟨photon_i, photon_j⟩, returns a matrix V
where V[i,k] is the coefficient of photon i in the k-th internal degree of freedom.

The Gram matrix can be reconstructed as: S = V * V'

# Arguments
- `S::Matrix`: n×n Gram matrix (Hermitian, positive semi-definite, diagonal = 1)
- `atol::Float64`: Absolute tolerance for rank detection (default: 1e-10)

# Returns
- `V::Matrix`: n×r matrix where r is the detected rank, V[i,k] = coefficient of photon i in internal mode k

# Algorithm
1. Verifies that S is a valid Gram matrix
2. Attempts Cholesky decomposition: S = L * L'
3. If Cholesky fails (rank-deficient), falls back to eigenvalue decomposition
4. Automatically detects effective rank by filtering near-zero columns/eigenvalues
5. Returns V matrix representing the internal degrees of freedom

# Example
```julia
S = rand_gram_matrix_from_orthonormal_basis(3, 2)
V = gram_to_coefficients(S)
S_reconstructed = V * V'  # Should equal S
```
"""
function gram_to_coefficients(S::Matrix; atol=1e-10)
    # Check that S is a valid Gram matrix
    check_is_gram_matrix(S, atol)

    n = size(S, 1)

    # Attempt Cholesky decomposition: S = L * L'
    # L is lower triangular matrix
    chol = cholesky(S, check=false)

    # Check if Cholesky decomposition succeeded
    if chol.info != 0
        # Cholesky failed (likely rank-deficient or numerical issues)
        # Fall back to eigenvalue decomposition

        eig = eigen(Hermitian(S))
        eigenvals = real.(eig.values)
        eigenvecs = eig.vectors

        # Keep only positive eigenvalues above threshold
        positive_idx = eigenvals .> atol

        # V = eigenvectors * sqrt(eigenvalues)
        V = eigenvecs[:, positive_idx] * Diagonal(sqrt.(eigenvals[positive_idx]))
    else
        # Cholesky succeeded
        L_full = chol.L  # n×n lower triangular matrix

        # Determine effective rank from L by finding non-zero columns
        col_norms = [norm(L_full[:, k]) for k in 1:n]
        r_effective = sum(col_norms .> atol)

        # Extract V as n×r_effective matrix
        V = L_full[:, 1:r_effective]
    end

    return V
end

"""
    reconstruct_gram_matrix(V::Matrix)

Reconstruct the Gram matrix from internal degrees of freedom coefficients.

# Arguments
- `V::Matrix`: n×r matrix where V[i,k] = coefficient of photon i in internal mode k

# Returns
- `S::Matrix`: n×n reconstructed Gram matrix where S[i,j] = Σₖ V[i,k] * conj(V[j,k])

# Example
```julia
V = gram_to_coefficients(S)
S_reconstructed = reconstruct_gram_matrix(V)
@assert S ≈ S_reconstructed
```
"""
function reconstruct_gram_matrix(V::Matrix)
    return V * V'
end
