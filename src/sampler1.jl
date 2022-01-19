include("permanents.jl")
include("scattering.jl")
include("special_matrices.jl")


function sample(Nmodes::Int, Nphotons::Int, x=0.8, η=1)

    """ computes approximate & exact distribution from boson samplers
        x: distinguishability
        η: transmission probablity """

    inputmodes = random_occupancy(Nmodes, Nphotons)
    input_occupancy_mode = occupancy_vector_to_partition(inputmodes)
    U = rand_haar(Nmodes)

    Nsample = 10^5
    n = length(input_occupancy_mode)
    N = size(U)[1]
    k = Nphotons-1
    kmax = k-1

    Nlist = []
    for i in product(Base.OneTo.(tuple([N for j=1:k]...))...)
        push!(Nlist, collect(i))
    end

    Papprox = zeros(length(Nlist))
    Pexact = zeros(length(Nlist))
    combs = collect(combinations(input_occupancy_mode, k))

    function submat(rows, cols)
        u = Array{typeof(U[1,1])}(undef, length(rows), length(cols))
        for i = 1:length(rows)
            for j = 1:length(cols)
                u[i,j] = U[rows[i], cols[j]]
            end
        end
        return u
    end

    for i = 1:length(Nlist)
        for j = 1:length(combs)
            perm = collect(permutations(combs[j]))
            for l = 1:length(perm)
                count = sum(collect(Int(perm[l][ll]==combs[j][ll]) for ll=1:length(perm[l])))
                korder = k-count

                pterm = x^korder * ryser_fast(submat(combs[j], Nlist[i]) .* conj(submat(perm[l], Nlist[i])))
                pterm = real(pterm/(binomial(n, k)*factorial(k)))

                Pexact[i] += pterm
                if korder <= kmax
                    Papprox[i] += pterm
                end
            end
        end
    end

    return Pexact, Papprox
end

Pe, Pa = sample(10, 10, 1, 1)
plot(1:length(Pe), Pe)
plot!(1:length(Pa), Pa)
