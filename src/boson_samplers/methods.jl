function sampling(;input::Input, interf::Interferometer, distinguishability::Real, reflectivity::Real, exact=nothing)

    if exact == nothing || exact == true

        if distinguishability == 1 && reflectivity == 1
            return cliffords_sampler(input=input, interf=interf)
        elseif distinguishability == 0 && reflectivity == 1
            return classical_sampler(input=input, interf=interf)
        else
            throw(ArgumentError("need approximation"))
        end

    else

        if distinguishability == 1 && reflectivity == 1
                  n = input.r.n
                  m = input.r.m
                  U = copy(interf.U)

                  # generate a collisionless state as a starting point
                  starting_state = iterate_until_collisionless(() -> random_occupancy(n,m))

                  known_pdf(state) = process_probability_distinguishable(U, input_state, state)
                  target_pdf(state) = process_probability(U, input_state, state)
                  known_sampler = () -> iterate_until_collisionless(() -> classical_sampler(U = U, m = m, n = n)) # gives a classical sampler

                  return metropolis_sampler(;target_pdf = target_pdf, known_pdf = known_pdf , known_sampler = known_sampler , starting_state = starting_state, n_iter = 100)
        else
            return noisy_sampling(input=input, distinguishability=distinguishability, reflectivity=reflectivity, interf=interf)
        end

    end

end
