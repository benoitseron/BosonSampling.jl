include("gaussian_partition.jl")


begin

    m = 20
    input_state = GeneralGaussian(m = m, r = 0.4 * ones(m))
    interferometer = RandHaar(m)
    part = equilibrated_partition(m, 2)

    # part = Partition(Subset(first_modes(1, m)))
    # part.subsets[1].subset

    n_max = 20
    mc = compute_probabilities_partition_gaussian(interferometer, part, input_state, n_max)

    # bar(real(pdf))

    sort_samples_total_photon_number_in_partition!(mc)

    @show mc

    bar(mc.proba)

end


### separating events by photon number ###

sorted_counts = sort_by_detected_photons(mc)

n_detected = 14

bar(sorted_counts[n_detected].proba)


### cutoff ###

m = 300
input_state = GeneralGaussian(m = m, r = 1.5 * ones(m))

"""

    find_cutoff(input_state::GeneralGaussian; atol = ATOL)

Given a value `atol`, finds the number of photons that is less likely to be detected than `atol`. Follows S3 in 2110.0156. 
"""
function find_cutoff(input_state::GeneralGaussian; atol = ATOL)

    n_in = input_state.m
    # check that all the elements of r are the same
    @argcheck all(input_state.r .== input_state.r[1]) "r must be the same for all modes"
    r = input_state.r[1]

    # @show n_in,r
    cutoff(α, n_in::Int, r::Real) = α * n_in * sech(r)^2

    bound(α) = exp(-0.5 *α * n_in * (1-1/α)^2)

    f = x -> bound(x) - atol

    α = find_zero(f, 1.) 

    bound(α)
    cutoff(α, n_in, r)

end

########### need to be modified - it's for photon pairs !

find_cutoff(input_state; atol = 1e-10)



m = 20
input_state = GeneralGaussian(m = m, r = 1.5 * ones(m))

function probability_n_photons(n::Int, input_state::GeneralGaussian)

    m = input_state.m

    @argcheck m % 2 == 0 "only implemented for even m, still holds but need to use gamma functions instead of factorials"
    @argcheck all(input_state.r .== input_state.r[1]) "r must be the same for all modes"
    r = input_state.r[1]
    
    p(k, m) = binomial(div(m, 2) + k - 1, k) * sech(r)^m * tanh(r)^(2k)

    if n % 2 == 0
        return p(div(n,2), m)
    else
        return 0.
    end

end

k_range = 0:200

foo(k) = probability_n_photons(k, input_state)

bar(k_range, foo.(k_range))
