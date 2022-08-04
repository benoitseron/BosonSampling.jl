using BosonSampling
using Plots


###### checking if linear around bosonic ######

# in Totally Destructive Many-Particle Interference
# they claim that around the bosonic case, when p_bos = 0
# the perturbation is linear and proportional to p_d
# we believe this is false

n = 8
m = n

interf = Fourier(m)

ib = Input{OneParameterInterpolation}(first_modes(n,m),1)
id = Input{OneParameterInterpolation}(first_modes(n,m),0)

# write a suppressed event
state = zeros(Int, n)
state[1] = n - 1
state[2] = 1
state = ModeOccupation(state)
o = FockDetection(state)

evb = Event(ib, o, interf)
evd = Event(id, o, interf)

events = [evb, evd]

for event in events
    compute_probability!(event)
end

p_b = evb.proba_params.probability
p_d = evd.proba_params.probability

function p(x)

    i = Input{OneParameterInterpolation}(first_modes(n,m),x)
    ev = Event(i, o, interf)
    compute_probability!(ev)

    ev.proba_params.probability

end

p_claim(x) = n*(1-x)* p_d

x_array = [x for x in range(0.9,1,length = 100)]
p_x_array = p.(x_array)
p_claim_array = p_claim.(x_array)

plt = plot(x_array, p_x_array, label = "p_x")
plot!(x_array, p_claim_array, label = "p_claim")
title!("Fourier, n = $(n)")

savefig(plt, "src/certification/images/check_dittel_approximation_fourier_n=$(n).png")
