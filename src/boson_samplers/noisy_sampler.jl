"""
    noisy_sampler(;input::Input, reflectivity::Real, interf::Interferometer)

Sample partially-distinguishable photons through a lossy interferometer, which
runs (at most) in ``O(n2^m + Poly(n,m))`` time.

!!! note "Reference"
    [https://arxiv.org/pdf/1907.00022.pdf](https://arxiv.org/pdf/1907.00022.pdf)
"""
function noisy_sampler(;input::Input, reflectivity::Real, interf::Interferometer)

    list_assignement = fill_arrangement(input.r.state)
    l = rand(Binomial(input.n, reflectivity))
    remaining_photons = Int.(zeros(input.m))
    remaining_subset = rand(collect(multiset_combinations(list_assignement, l)))
    remaining_photons[remaining_subset] .= 1

    list_assignement = fill_arrangement(remaining_photons)
    i = rand(Binomial(l, input.distinguishability_param))
    bosonic_input = Int.(zeros(input.m))
    bosonic_subset = rand(collect(multiset_combinations(list_assignement, i)))
    bosonic_input[bosonic_subset] .= 1

    classical_input = input.r.state .- bosonic_input
    classical_input = Input{Distinguishable}(ModeOccupation(classical_input))
    bosonic_input = Input{Bosonic}(ModeOccupation(bosonic_input))

    classical_output = fill_arrangement(classical_sampler(input=classical_input, interf=interf))
    bosonic_output = cliffords_sampler(input=bosonic_input, interf=interf)

    return sort(append!(classical_output, bosonic_output))

end
