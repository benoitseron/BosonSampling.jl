using Combinatorics
using LinearAlgebra

# abstract type inputMat{T} <: AbstractMatrix{T} end
# abstract type unitary{T} <: inputMat{T} end
# abstract type hsdp{T} <: inputMat{T} end

function permanent_def(U::AbstractMatrix{T}) where {T}

    """ computes the permanent using the definition """

    n = size(U)[1]
    sum_diag(perm, U) = prod(U[i,perm[i]] for i = 1:n)

    return sum(sum_diag(perm, U) for perm in permutations(collect(1:n)))
end

function ryser_slow(U::AbstractMatrix{T}) where {T}

    """ computes the permanent using Ryser formula O(n^2×2^n) """

    n = size(U)[1]

    function submatrix(k::Integer)::Matrix{typeof(U[1,1])}[]

        """ returns the set of n×k submatrices of U """

        set = Matrix{typeof(U[1,1])}[]
        x = vcat(ones(Int8, k), zeros(Int8, n-k))
        perms = unique(collect(permutations(x)))

        for i = 1:length(perms)
            count = 0
            submat = zeros(typeof(U[1,1]), n, k)
            for j = 1:length(perms[i])
                if perms[i][j]≠0
                    count+=1
                    submat[:,count] = U[:,j]
                end
            end
            push!(set, submat)
        end

        return set
    end

    sum_rows(k::Integer) = sum(prod(sum(u[j,l] for l = 1:k) for j = 1:n) for u in submatrix(U, k))

    return sum((-1)^k * sum_rows(U, n-k) for k = 0:n-1)
end

function ryser_fast(U::AbstractMatrix{T}) where {T}

    """ computes permanent by Ryser formula using Gray ordering O(2^n×n) """
    # code from https://numbersandshapes.net/posts/permanents_and_rysers_algorithm/

    n = size(U)[1]

    gd = zeros(Int64, 1, 2^n-1)
    gd[1] = 1

    v = vcat([(2^i,i+1,1) for i in 0:n-1], [(-2^i,i+1,-1) for i in 0:n-1])
    lut = Dict(x[1] => (x[2],x[3]) for x in v)

    r = U[1,:]
    s = -1
    pm = s*prod(r)

    for i in (2:2^n-1)
        if iseven(i)
            gd[i] = 2*gd[div(i,2)]
        else
            gd[i] = (-1)^((i-1)/2)
        end
        r += U[lut[gd[i]][1],:]*lut[gd[i]][2]
	   	s *= -1
	   	pm += s*prod(r)
	end

	return pm * (-1)^n
end

permanent_ryser = ryser_fast

function quantumInspired(Λ::AbstractMatrix{T}, ϵ=10^(-4), δ=10^(-4), C=2) where {T}

	""" Approximates the permanent of Λ with an error ϵ and a probability of
		failure δ in N steps """
		# https://arxiv.org/pdf/1609.02416.pdf

	if isposdef(Λ)
		n = size(Λ)[1]
		U = eigvecs(Λ)
		D = eigvals(Λ)
		λmax = maximum(D)

		Z = (C*λmax)^(2*n) / prod(C*λmax-D[i] for i = 1:n)
		N = trunc(Int, log(1/δ)*(Z^2*exp(-2*n))/(2*ϵ^2))
		pcs = Array{Float64}(undef, N)
		α = [1/sqrt(D[i]) for i = 1:n]

		for j = 1:N
		    α = α .* randn(eltype(α), n)
		    β = U'*α
		    pcs[j] = prod(exp(-abs(b)^2) * abs(b)^2 for b in β)
		end

		μ = sum(pcs[j]/n for j = 1:n)

		return Z * μ
	else
		throw(ArgumentError("permanents: argument must be a HPSD matrix"))
	end
end

function gurvits(U::AbstractMatrix{T}, m=10000) where {T}

	""" Approximate perm(U) using Glynn's formula """
	# https://arxiv.org/pdf/1212.0025.pdf

	n = size(U)[1]
	res = 0.0

	function glynn_estimator(x::AbstractArray)
		return prod(x[i] for i = 1:n) * prod(sum(U[i, :] .* x) for i = 1:n)
	end

	for i = 1:m
		x = collect(rand([-1,1], n))
		res += glynn_estimator(x)
	end

	return res/m
end
