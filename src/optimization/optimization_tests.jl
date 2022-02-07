n = 10
A = Diagonal(range(1, stop=2, length=n))
f(x) = dot(x,A*x)/2
g(x) = A*x
g!(stor,x) = copyto!(stor,g(x))
x0 = randn(n)

manif = Optim.Sphere()
Optim.optimize(f, g!, x0, Optim.ConjugateGradient(manifold=manif))

n = 10
A = Diagonal(range(1, stop=2, length=n))
f(x) = dot(x[1:n],A*x[1:n])/2 + dot(x[n+1:2n],A*x[n+1:2n])/2
x0 = randn(2n)

manif = Optim.ProductManifold(Optim.Sphere(), Optim.Sphere(), (n,), (n,))
res = Optim.optimize(f, x0, Optim.ConjugateGradient(manifold=manif))

v =Optim.minimizer(res)

norm(v[1:n])
norm(v[n+1:end])

norm(v)

#### iterative manifold ####

# r column vectors noramlized of dimension n

n = 10
r = n

manif = Optim.Sphere()
for layer = 1:r-1
    manif = Optim.ProductManifold(Optim.Sphere(), manif, (n,), (layer * n,))
end


function array_to_matrix(array, r, n)
    """gives a (r,n) matrix from the r, n-sized vectors concatenated in array"""

    mat = Matrix{eltype(array)}(undef, r, n)

    for line = 1 : r
        for col = 1 : n
            mat[line, col] = array[(line-1)*n+col]
        end
    end

    mat
end

function objective(array, r, n)

    M = array_to_matrix(array, r, n)

    gram_matrix = M'*M

    -abs(permanent_ryser(gram_matrix))

end


function gradient_permanent!(G, array, n, r)

    A = array_to_matrix(array, r, n)

    for i = 1:n
        for j = 1:n
            grad[i,j] = permanent_ryser(remove_row_col(A, i, j))
        end
    end

end


function foo(array)
    objective(array,r,n)
end

function g!(G, array)
    gradient_permanent!(G, array, n, r)
end

res = Optim.optimize(foo,rand_gram_matrix(n), Optim.GradientDescent(manifold=manif), Optim.Options(x_tol = 1e-16, iterations = 100000))
Float64(factorial(n))
v = Optim.minimizer(res)

println("###############")
for i = 1 : n
    println(real(norm(v[i,:])))
end
