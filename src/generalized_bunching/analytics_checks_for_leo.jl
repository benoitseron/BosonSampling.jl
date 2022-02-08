function leo_ratio(n)

    numerator(n) = (factorial(n) + 0.25*n^2 * factorial(n-2))/(n^n)

    U = fourier_matrix(n)
    H = H_matrix(U, [1 for i in 1:n], [i <= 2 ? 1 : 0 for i in 1:n])
    omega_n = exp(2*pi*1im/n)

    pol_states = 1/sqrt(2) * ones(ComplexF64, 2,n)
    pol_states[2,:] .= 1/sqrt(2) * [omega_n^(i-1) for i in 1:n]

    S = pol_states' * pol_states

    @test is_a_gram_matrix(S)

    numerator(n)/clean_proba(ryser(H .* S))

end

n_min = 3
n_max = 27
n_array = big.(collect(n_min:n_max))

plotly();
plt = plot(n_array, leo_ratio.(n_array))

@show (n_array, leo_ratio.(n_array))
plt

savefig(plt,"src/generalized_bunching/analytics_checks_for_leo.html")
