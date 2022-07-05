using BosonSampling

abstract type NonLinearInterferometer end
# not set as a subtype of Interferometer so that functions do not use it as
# a linear interferometer

########## if you do NonLinearInterferometer <: Interferometer check that we don't make assumptions of unitarity or how to remove them!

struct SingleModeNonLinearPhaseShift
    m::Int
    V::Matrix
    W::Matrix
    ϕ::Real

    # if V and W are given
    SingleModeNonLinearPhaseShift(V::Matrix,W::Matrix,ϕ) = size(V) == size(W) ? new(size(V,1), V,W, ϕ) : error("V and W do not have the same size")

    # if they are not given, sample them from the haar measure
    SingleModeNonLinearPhaseShift(m::Int,ϕ) = new(m, rand_haar(m),rand_haar(m), ϕ)
end

n = 2
m = 4
i = Input{Bosonic}(first_modes(n,m))
o = FockDetection([1,0,1,0])
interf = SingleModeNonLinearPhaseShift(m, π)

println(".........................")
for config in ranked_partition_list(all_mode_configurations(n,m, only_photon_number_conserving = true))
    println(config)
end
############# not yet what we want

#
# function compute_probability!(ev::Event{TIn,TOut}) where {TIn<:Bosonic, TOut<:FockDetection}
#
# 	# warn users if a computation was already performed
# 	check_probability_empty(ev)
#
# 	# set precision parameters of the current computation as an exact computation
# 	ev.proba_params.precision = eps()
# 	ev.proba_params.failure_probability = 0
#
# 		ev.proba_params.probability = bosonic_probability(ev.interferometer.U, ev.input_state.r.state, ev.output_measurement.s.state)
