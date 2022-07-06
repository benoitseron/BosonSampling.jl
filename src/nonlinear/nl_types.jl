using BosonSampling

abstract type NonLinearInterferometer end
# not set as a subtype of Interferometer so that functions do not use it as
# a linear interferometer

########## if you do NonLinearInterferometer <: Interferometer check that we don't make assumptions of unitarity or how to remove them!

struct SingleModeNonLinearPhaseShift <: NonLinearInterferometer
    m::Int
    V::Matrix
    W::Matrix
    ϕ::Real
    x::Int #mode

    # if V and W are given
    SingleModeNonLinearPhaseShift(V::Matrix,W::Matrix,ϕ,x) = size(V) == size(W) ? new(size(V,1), V,W, ϕ, x) : error("V and W do not have the same size")

    # if they are not given, sample them from the haar measure
    SingleModeNonLinearPhaseShift(m::Int,ϕ,x) = new(m, rand_haar(m),rand_haar(m), ϕ, x)
end

n = 2
m = 4
i = Input{Bosonic}(first_modes(n,m))
o = FockDetection([1,0,1,0])

ϕ = π
x = 1
interf = SingleModeNonLinearPhaseShift(m, π, 1)

ev = Event(i,o,interf)


function single_mode_nl_phase_shift_probability(ev::Event{TIn, TInterf, TOut}) where {TIn<:InputType, TInterf<:SingleModeNonLinearPhaseShift, TOut<:OutputMeasurementType}

	@unpack m, V,W,x,ϕ = interf

	input_state = i.r.state
	output_state = o.s.state

	n = i.r.n

	inner_term(r) = exp(-1im * r[x]^2 * ϕ) * permanent(scattering_matrix(V, r, output_state)) * permanent(scattering_matrix(W, input_state, r))/sqrt(vector_factorial(input_state) * vector_factorial(output_state) * vector_factorial(r)^2)

	abs(sum([inner_term(r) for r in all_mode_configurations(n,m, only_photon_number_conserving = true)]))^2

end

function compute_probability!(ev::Event{TIn, TInterf, TOut}) where {TIn<:InputType, TInterf<:NonLinearInterferometer, TOut<:OutputMeasurementType}

	# warn users if a computation was already performed
	check_probability_empty(ev)

	# set precision parameters of the current computation as an exact computation
	ev.proba_params.precision = eps()
	ev.proba_params.failure_probability = 0

	ev.proba_params.probability = single_mode_nl_phase_shift_probability(ev)

end
