function classical_sampler(;U,n,m)
    """outputs a sample of output mode of the classical distribution
    (distinguishable particles) form the input of [1^n O^(m-n)]
    according to section 8.1 of https://arxiv.org/abs/1309.7460"""

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
