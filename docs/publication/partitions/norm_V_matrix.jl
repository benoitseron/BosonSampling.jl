n_iter = 100000
m = 20

for iter in n_iter


    V = RandHaar(m).U
    S = rand_gram_matrix(m)

    eigen(V).values

    if any(abs.(eigen(V .* S).values) .> 1)
        @show V
        @show S
        @show abs.(eigen(V .* S).values)

    end

end
