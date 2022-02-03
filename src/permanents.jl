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

function ryser_fast(U::AbstractMatrix{T}) where T

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

function permanent_ryser(A::AbstractMatrix{T}) where T

	"""computes the permanent of A using ryser with Gray ordering"""
	# code from https://discourse.julialang.org/t/matrix-permanent/10766
	# see Combinatorial Algorithms for Computers and Calculators, Second Edition (Computer Science and Applied Mathematics) by Albert Nijenhuis, Herbert S. Wilf (z-lib.org)
	# chapter 23 for the permanent algorithm
	# chapter 1 for the gray code

    function grayBitToFlip(n::Int)
    	n1 = (n-1) ⊻ ((n-1)>>1)
    	n2 = n ⊻ (n>>1)
    	d = n1 ⊻ n2

    	j = 0
    	while d>0
    		d >>= 1
    		j += 1
    	end
    	j
    end

    n,m = size(A)
	if (n == m)
        AA = 2*A
		D = true
		rho = [true for i = 1:n]
		v = sum(A, dims = 2)
		p = prod(v)
		@inbounds for i = 1:(2^(n-1)-1)
			a = grayBitToFlip(i)+1
            if rho[a]
                @simd for j=1:n
                    v[j] -= AA[j,a]
                end
            else
                @simd for j=1:n
                    v[j] += AA[j,a]
                end
            end
            rho[a] = !rho[a]
            pv = one(typeof(p))
            @simd for j=1:n
                pv *= v[j] #most expensive part
            end
			D = !D
            if D
                p += pv
            else
                p -= pv
            end
		end

		return p * 2.0^(1-n)
	else
		throw(ArgumentError("perm: argument must be a square matrix"))
	end
end

function multi_dim_ryser(U, gram_matrix) # Still neeed to be verified

    # https://arxiv.org/pdf/1410.7687.pdf

    n = size(U)[1]
    nstring = collect(1:n)
    sub_nstring = collect(powerset(nstring))
    sub_nstring = sub_nstring[2:length(sub_nstring)]
    res = 0

    function delta_set(S1, S2)
        if S1 == S2
            return 1
        else
            return 0
        end
    end

    for r = 1:length(sub_nstring)
        for s = r:length(sub_nstring)
            R = sub_nstring[r]
            S = sub_nstring[s]

            t = prod(
                sum(
                    U[ss, j] * conj(U[rr, j]) * gram_matrix[rr, ss] for rr in R for ss in S
                ) for j = 1:n
            )
            res +=
                (2 - delta_set(S, R)) * (-1)^(length(S) + length(R)) * real(t)
        end
    end
    return res
end

# function quantumInspired(Λ::AbstractMatrix{T}, ϵ=10^(-4), δ=10^(-4), C=2) where {T}
#
# 	""" Approximates the permanent of Λ with an error ϵ and a probability of
# 		failure δ in N steps """
# 		# https://arxiv.org/pdf/1609.02416.pdf
#
# 	if isposdef(Λ)
# 		n = size(Λ)[1]
# 		U = eigvecs(Λ)
# 		D = eigvals(Λ)
# 		λmax = maximum(D)
#
# 		Z = (C*λmax)^(2*n) / prod(C*λmax-D[i] for i = 1:n)
# 		N = trunc(Int, log(1/δ)*(Z^2*exp(-2*n))/(2*ϵ^2))
# 		pcs = Array{Float64}(undef, N)
# 		α = [1/sqrt(D[i]) for i = 1:n]
#
# 		for j = 1:N
# 		    α = α .* randn(eltype(α), n)
# 		    β = U'*α
# 		    pcs[j] = prod(exp(-abs(b)^2) * abs(b)^2 for b in β)
# 		end
#
# 		μ = sum(pcs[j]/n for j = 1:n)
#
# 		return Z * μ
# 	else
# 		throw(ArgumentError("permanents: argument must be a HPSD matrix"))
# 	end
# end

function gurvits(U::AbstractMatrix{T}, ϵ=1e-3) where {T}

	""" Approximate perm(U) using Glynn's formula """
	# https://arxiv.org/pdf/1212.0025.pdf

	n = size(U)[1]
	m = 10^(trunc(Int, log(10, ϵ^(-2))))

	function glynn_estimator(x::AbstractArray)
		return prod(x[i] for i = 1:n) * prod(sum(U[i, :] .* x) for i = 1:n)
	end

	term_ = []
	dist_r = 1
	dist_i = 1

	for i = 1:m
		push!(term_, glynn_estimator(rand([-1,1], n)))
	end
	perm_i = sum(term_)/length(term_)

	while dist_r > ϵ && dist_i > ϵ

		push!(term_, glynn_estimator(rand([-1,1], n)))
		perm_f = sum(term_)/length(term_)

		dist_r = abs(real(perm_f)-real(perm_i))
		dist_i = abs(imag(perm_f)-imag(perm_i))
		# println(dist_r, dist_r)
		perm_i = perm_f

	end

	return perm_i

end

function fast_glynn_perm(U::AbstractMatrix{T}) where T

	size(U)[1] == size(U)[2] ? n=size(U)[1] : error("Non square matrix as input")

	row_ = [U[:,i] for i in 1:n]
	row_comb = [sum(row) for row in row_]

	res = 0
	old_gray = 0
	sign = +1
	binary_power_dict = Dict(2^k => k for k in 0:n)
	num_iter = 2^(n-1)

	for bin_index in 1:num_iter
		res += sign * reduce(*, row_comb)

		new_gray = bin_index ⊻ trunc(Int, bin_index/2)
		gray_diff = old_gray ⊻ new_gray
		gray_diff_index = binary_power_dict[gray_diff]

		new_vec = U[gray_diff_index+1,:]
		direction = 2 * cmp(old_gray, new_gray)

		for i in 1:n
			row_comb[i] += new_vec[i] * direction
		end

		sign = -sign
		old_gray = new_gray
	end

	return res/num_iter

end

# function jsv(A::AbstractMatrix{T}) where {T}
#
# 	""" using the fact that A with 0,1 entries can be seen as the adjancy matrix
# 	 	of graph G and its permanent is the number of perfect matchings in G """
#
# 	# https://faculty.cc.gatech.edu/~vigoda/Permanent.pdf
#
# 	A = real(make_matrix_more_readable(A))
# 	if size(A)[1] != size(A)[2]
# 		throw(ArgumentError("requiere square matrix as input"))
# 	end
#
# 	any(a->a < 0, A) ? throw(ArgumentError("requiere matrix with positive entries as input")) : nothing
#
# 	η = 1e-4 # overall failure probability
# 	η_ = 10.0^(trunc(Int, log(10, η/log(n^4*log(n))))) # failure probability at each estimation of the mean from a sample
# 	δ = 10.0^(trunc(Int, log(10, n^(-2)))) # error in variation distance
# 	sample_size = 10^(trunc(Int, log(10, n^2*log(1/η_))))
# 	number_step = 10^(trunc(Int, log(10, n^7*log(n))))
#
# 	max_a = maximum(A)
# 	min_a = minimum(A)
# 	n = size(A)[1]
#
# 	for i = 1:n
# 		for j = 1:n
# 			if A[i,j] == 0
# 				A[i,j] = min_a^n/factorial(n)
# 			end
# 		end
# 	end
#
# 	edge_activity = [[max_a for i = 1:n] for j = 1:n]
# 	hole_weigth = [[n*max_a for i = 1:n] for j = 1:n]
# 	ans_ = [edge_activity[i][j] > A[i,j] for (i,j) in product(1:n, 1:n)]
#
# 	function markov_chain(edge_activity, hole_weigth)
#
#
# 	while any(ans_)
#
#
#
# 	end
# end
