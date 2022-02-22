# implements a metropolis_independant_sampler that does standard boson sampling
# following https://arxiv.org/abs/1705.00686

# the paper is limited to collisionless events so we keep track of this in the
# following through iterate_until_collisionless


function metropolis_sampler(;target_pdf, known_pdf, known_sampler, starting_state, n_iter, n_burn = 100, n_thinning = 100)

    function transition_probability(target_pdf, known_pdf, new_state, previous_state)

        """metropolis_independant_sampler transition probability, see eq 1 in https://arxiv.org/abs/1705.00686"""
        min(1, target_pdf(new_state) * known_pdf(previous_state)/(target_pdf(previous_state) * known_pdf(new_state)))
    end

    function do_with_probability(p)

        """returns true with probability p, false with (1-p)"""

        rand() < p ? true : false

    end

    function update_state_parameters(target_pdf, known_pdf, new_state, previous_state)

        if do_with_probability(transition_probability(target_pdf, known_pdf, new_state, previous_state))

            return new_state
        else
            return previous_state
        end

    end

    function metropolis_iteration!(target_pdf, known_pdf, known_sampler, new_state, previous_state)
        new_state[:] = known_sampler()[:]
        new_state[:] = update_state_parameters(target_pdf, known_pdf, new_state, previous_state)[:]
        previous_state[:] = new_state[:]
    end

    samples = Array{eltype(starting_state)}(undef, n_iter, size(starting_state,1))

    previous_state = similar(starting_state) #x
    new_state = similar(previous_state)
    previous_state = starting_state

    if n_burn == n_thinning
        @showprogress for iter in 1:n_iter
            for i in 1:n_thinning
                metropolis_iteration!(target_pdf, known_pdf, known_sampler, new_state, previous_state)
            end
            samples[iter,:] = new_state[:]
        end
    else
        throw(Exception("n_burn != n_thinning not implemented"))
    end

    samples

end
