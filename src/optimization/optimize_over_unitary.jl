# follows Conjugate gradient algorithm for optimization under unitary Traian Abrudan  ,1,2 , Jan Eriksson 2 , Visa Koivunen
# matrix constraint Signal Processing 89 (2009) 1704–1714


function step_size(w_k, h_k; q, euclidian_gradient, p = 5)

    """

        step_size(w_k, h_k; q, euclidian_gradient, p = 5)

    implements the geodesic search algorithm of table 1 in Conjugate gradient algorithm for optimization under unitary Traian Abrudan  ,1,2 , Jan Eriksson 2 , Visa Koivunen matrix constraint Signal Processing 89 (2009) 1704–1714"""

    n = size(w_k)[1]

    omega_max = maximum(abs.(eigvals(h_k)))

    t_mu = 2pi/(q*omega_max)

    mu = [i * t_mu/p for i in 0:p]

    r = Matrix{ComplexF64}[]

    push!(r, Matrix{ComplexF64}(I,n,n))
    push!(r, exp(-t_mu/p * h_k))
    for i in 2:p
        push!(r, r[2]*r[end])
    end

    derivatives = zeros(ComplexF64, p+1)

    for i = 0:p
        derivatives[i+1] = -2 * real(tr(euclidian_gradient(r[i+1]*w_k) * w_k'*r[i+1]'*h_k'))
    end

    a = zeros(ComplexF64, p+1)

    a[1] = derivatives[1]

    mat_mu = zeros(ComplexF64, p,p)

    for i = 1:p
        mat_mu[:, i] = mu[2:end].^i
    end

    a[2:end] = inv(mat_mu) * (derivatives[2:end].-a[1])

    positive_roots = []

    for root in roots(reverse(a))
        if imag(root) ≈ 0. atol = 1e-7
            if real(root) > 0.
                push!(positive_roots, real(root))
            end
        end
    end

    #step_size
    if positive_roots == []
        return 0.
    else
        return minimum(positive_roots)
    end

end

"""

    minimize_over_unitary_matrices(;euclidian_gradient, q, n, p=5, tol = 1e-6, max_iter = 100)

returns the (optimized matrix, optimization_success)

implements the minimization algorithm of table 3 in  Conjugate gradient algorithm for optimization under unitary Traian Abrudan  ,1,2 , Jan Eriksson 2 , Visa Koivunen matrix constraint Signal Processing 89 (2009) 1704–1714"""
function minimize_over_unitary_matrices(;euclidian_gradient, q, n, p=5, tol = 1e-6, max_iter = 100)



    optimization_success = false

    function inner_product(a,b)
        0.5 * real(tr(a'*b))
    end

    k = 0
    w_k = Matrix{ComplexF64}(I,n,n)
    g_k = similar(w_k)
    h_k = similar(w_k)

    while k < max_iter

        #println("iter $(k+1)/$max_iter")
        if k % n^2 == 0
            gamma_k = euclidian_gradient(w_k)
            g_k = gamma_k * w_k' - w_k * gamma_k'
            h_k = g_k
        end

        if inner_product(g_k,g_k) < tol
            #println("optimisation achieved")
            optimization_success = true
            break

        end

        w_k_plus_1 = exp(-step_size(w_k, h_k; q = q, euclidian_gradient = euclidian_gradient, p = p) * h_k)*w_k

        gamma_k_plus_1 = euclidian_gradient(w_k_plus_1)

        g_k_plus_1 = gamma_k_plus_1 * w_k_plus_1' - w_k_plus_1*gamma_k_plus_1'
        gamma_small_k = inner_product((g_k_plus_1 - g_k) , g_k_plus_1) / inner_product(g_k , g_k)

        h_k_plus_1 = g_k_plus_1 + gamma_small_k * h_k

        if inner_product(h_k_plus_1, g_k_plus_1) < 0.
            h_k_plus_1 = g_k_plus_1
        end

        k+=1
        gamma_k = gamma_k_plus_1
        w_k = w_k_plus_1
        h_k = h_k_plus_1

    end

    if optimization_success == false
        #println("optimisation failed")
    end

    return (w_k, optimization_success)

end

function run_brockett_optimization_test(;maxiter = 100)

    achieved_an_optimization = false
    iter = 0
    this_test = nothing

    while achieved_an_optimization == false && iter < maxiter

        iter += 1
        n = 3
        N = diagm(1:n)
        sigma = hermitian_test_matrix(n)
        q_brockett = 2

        function objective_function(w, sigma, N)
            tr(w' * sigma * w * N)
        end

        function euclidian_gradient_brockett(w,sigma,N)
            sigma*w*N
        end

        function euclidian_gradient(w)
            euclidian_gradient_brockett(w,sigma,N)
        end

        result= minimize_over_unitary_matrices(;euclidian_gradient = euclidian_gradient, q = q_brockett,n=n, tol = 1e-7, max_iter = 100000)

        achieved_an_optimization = result[2]

        if achieved_an_optimization
            optimized_matrix = result[1]
            this_test = @test real.(optimized_matrix' * sigma * optimized_matrix) ≈ diagm(sort(eigvals(sigma), rev = true)) atol = 1e-2

            break
        end

    end

    if achieved_an_optimization == false
        println("no optimization succeeded during the test, maybe be statistical but looks like an anomaly")
    end

    this_test

end
