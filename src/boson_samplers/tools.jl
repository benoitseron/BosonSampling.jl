function LaplaceExpansion(perm_mat, full_mat)

    function collect_sub_mat_perm(A)

        if length(size(A)) == 1
            v = [(A...)...]
            return [v[setdiff(1:end)] for i in 1:length(v)]
        else
            sub_mat = [remove_row_col(A, [], [i]) for i in 1:size(A)[2]]
            # return [fast_glynn_perm(a) for a in sub_mat]
            return [ryser(a) for a in sub_mat]
        end

    end

    res = []
    v_perms = collect_sub_mat_perm(perm_mat)

    for i in 1:size(full_mat)[1]
        rowA = full_mat[i,:]
        ans = 0
        for j in 1:length(v_perms)
            ans += rowA[j] * v_perms[j]
        end
        push!(res, ans)
    end

    return [abs(r)^2 for r in res]

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
