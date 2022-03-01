function noisy_sampling(;input::Input, distinguishability::Real, reflectivity::Real, interf::Interferometer)

    """ choose state truncation s.t the runtime goes as O(2^k + Poly(m,n,k)).
        k increases linearly as n and the error ϵ goes as x^n
        https://arxiv.org/pdf/1907.00022.pdf """

    input_modes = input.r.state
    n = input.r.n
    m = input.r.m
    U = interf.U

    k = trunc(Int, distinguishability*reflectivity*n)
    k == 0 ? k=1 : nothing

    l = rand(Binomial(n, reflectivity))
    l > k ? l=k : nothing
    i = rand(Binomial(l, distinguishability))
    remaining_photons = Int.(zeros(m))
    bosonic_input = Int.(zeros(m))

    count = 0
    while count != l
        rand_idx = rand(1:m, l)
        remaining_photons[rand_idx] = input_modes[rand_idx]
        count = sum(remaining_photons)
    end

    count = 0
    while count != i
        rand_idx = rand(1:m, i)
        bosonic_input[rand_idx] = remaining_photons[rand_idx]
        count = sum(bosonic_input)
    end

    input = [input_modes[i]-bosonic_input[i] for i = 1:m]
    classical_input = Input{Distinguishable}(ModeOccupation(input))
    classical_output = fill_arrangement(classical_sampler(input=classical_input, interf=interf))
    input = Input{Bosonic}(ModeOccupation(bosonic_input))
    bosonic_output = cliffords_sampler(input=input, interf=interf)

    return append!(classical_output, bosonic_output)

end
