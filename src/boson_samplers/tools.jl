function williamson_decomp(V::AbstractMatrix)
    @argcheck isposdef(V) && size(V)[1] == size(V)[2] && size(V)[1] % 2 == 0

    n = div(size(V)[1], 2)
    J = symplectic_mat(n)
    v = sqrt(inv(V))
    id = Matrix{Float64}(I, 2n, 2n)

    F = schur(v*J*v)
    d = F.values
    select = [(i>0) for i in imag(d)]
    ordschur!(F, select)

    K = F.Z
    D = diagm(inv.(F.values))
    S = inv(transpose(v * K * sqrt(D)))

    T = real.(1/2 * S * transpose(S))
    W = real.(S * (D - 1/2 * id) * transpose(S))
    return T, W
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
