function collect_sub_mat_perm(U)

    if size(U)[1] == 1 || size(U)[2] == 1
        return U
    else
        return [fast_glynn_perm(remove_row_col(U, [], [i])) for i in 1:size(U)[2]]
    end

end

function LaplaceExpansion(perm_mat, full_mat)
    permanents = collect_sub_mat_perm(perm_mat)
    row_ = [full_mat[i,:] for i in 1:size(full_mat)[1]]
    res = [dot(row, permanents) for row in row_]
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
