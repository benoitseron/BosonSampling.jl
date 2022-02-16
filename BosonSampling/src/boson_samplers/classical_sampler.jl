function classical_sampler(;input::Input, interf::Interferometer)

    """outputs a sample of output mode of the classical distribution
    (distinguishable particles) form the input of [1^n O^(m-n)]
    according to section 8.1 of https://arxiv.org/abs/1309.7460"""

    input_modes = input.r.state
    n = input.r.n
    m = input.r.m

    output_state = zeros(Int, m)
    output_modes = collect(1:m)

    U = interf.U

    for j in 1:n
        this_output_mode = wsample(output_modes, abs.(U[j,:]).^2)
        output_state[this_output_mode] += 1
    end

    output_state

end
