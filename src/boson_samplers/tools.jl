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
        return res
    end

end

"""
    williamson_decomp(V::AbstractMatrix)

Find the Williamson decomposition of the positive semidefinite matrix ``V``.
Returns ``D`` and ``S`` such that ``V = S D S^T
!!! Reference
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
