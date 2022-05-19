function classical_sampler(U, n, m)

    output_state = zeros(Int,m)
    output_modes = collect(1:m)

    #@warn "check U or U'"
    for j in 1:n
        # sample photons according to the distinguishable case, adds it to the output mode
        this_output_mode = wsample(output_modes, abs.(U[j,:]) .^2)
        output_state[this_output_mode] += 1
    end

    output_state

end

classical_sampler(input::Input, interf::Interferometer) = classical_sampler(interf.U, input.n, input.m)

classical_sampler(;input::Input, interf::Interferometer) = classical_sampler(interf.U, input.n, input.m)
