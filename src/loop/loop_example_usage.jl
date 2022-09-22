begin
    using Revise

    using BosonSampling
    using Plots
    using ProgressMeter
    using Distributions
    using Random
    using Test
    using ArgCheck
    using StatsBase
    using ColorSchemes
    using Interpolations
    using Dierckx
    using LinearAlgebra
    using PrettyTables
    using LaTeXStrings
    using JLD
    using AutoHashEquals
    using LinearRegression

    using DataStructures
end

###### one dimension example ######

# n = 1
# m = 1
# i = Input{Bosonic}(first_modes(n,m))
# o = FockDetection(first_modes(n,m))
# η_loss = 0.9
#
# circuit = LossyCircuit(1)
# interf = LossyLine(η_loss)
# target_modes = [1]
#
# add_element!(circuit, interf, target_modes = target_modes)

# function to_lossy(target_modes)

n = 1
m = 1

function lossy_line_example(η_loss)

    circuit = LossyCircuit(1)
    interf = LossyLine(η_loss)
    target_modes = [1]

    add_element_lossy!(circuit, interf, target_modes = target_modes)
    circuit

end

lossy_line_example(0.9)

transmission_amplitude_loss_array = 0:0.1:1
output_proba = []

i = Input{Bosonic}(to_lossy(first_modes(n,m)))
o = FockDetection(to_lossy(first_modes(n,m)))

for transmission in transmission_amplitude_loss_array

    ev = Event(i,o, lossy_line_example(transmission))
    @show compute_probability!(ev)
    push!(output_proba, ev.proba_params.probability)
end

plot(transmission_amplitude_loss_array, output_proba)
ylabel!("p no lost")
xlabel!("transmission amplitude")

### the same with autoconversion of the input and output dimensions ###

i = Input{Bosonic}(first_modes(n,m))
o = FockDetection(first_modes(n,m))

for transmission in transmission_amplitude_loss_array

    ev = Event(i,o, lossy_line_example(transmission))
    @show compute_probability!(ev)
    push!(output_proba, ev.proba_params.probability)
end

### 2d HOM with loss example ###

n = 2
m = 2
i = Input{Bosonic}(first_modes(n,m))
o = FockDetection(ModeOccupation([2,0])) # detecting bunching, should be 0.5 in probability if there was no loss
transmission_amplitude_loss_array = 0:0.1:1
output_proba = []

function lossy_bs_example(η_loss)

    circuit = LossyCircuit(2)
    interf = LossyBeamSplitter(1/sqrt(2), η_loss)
    target_modes = [1,2]

    add_element_lossy!(circuit, interf, target_modes = target_modes)
    circuit

end

for transmission in transmission_amplitude_loss_array

    ev = Event(i,o, lossy_bs_example(transmission))
    @show compute_probability!(ev)
    push!(output_proba, ev.proba_params.probability)
end

plot(transmission_amplitude_loss_array, output_proba)
ylabel!("p bunching top mode")
xlabel!("transmission amplitude")


### building the loop ###

n = 3
m = n
η = 1/sqrt(2) .* ones(m-1)
# 1/sqrt(2) .* [1,0] #ones(m-1) # see selection of target_modes = [i, i+1] for m-1
# [1/sqrt(2), 1] #1/sqrt(2) .* ones(m-1) # see selection of target_modes = [i, i+1] for m-1

# η_loss = 1. .* ones(m-1)

circuit = LosslessCircuit(m) #LossyCircuit(m)

for mode in 1:m-1
    @show mode
    interf = BeamSplitter(η[mode]) #LossyBeamSplitter(reflectivities[mode], η_loss[mode])
    target_modes_in = [mode, mode+1]
    target_modes_out = [mode, mode+1]
    add_element!(circuit, interf, target_modes_in = target_modes_in, target_modes_out = target_modes_out)

end

BeamSplitter(1/sqrt(2))

pretty_table(circuit.U)

i = Input{Bosonic}(first_modes(n,m))

#outputs compatible with two photons top mode
o1 = FockDetection(ModeOccupation([2,1,0]))
o2 = FockDetection(ModeOccupation([2,0,1]))

o_array = [o1,o2]

p_two_photon_first_mode = 0

for o in o_array
    ev = Event(i,o, circuit)
    @show compute_probability!(ev)
    p_two_photon_first_mode += ev.proba_params.probability
end

p_two_photon_first_mode

o3 = FockDetection(ModeOccupation([3,0,0]))
ev = Event(i,o3, circuit)
@show compute_probability!(ev)


### random phase shifter ###

d = Uniform(0, 2pi)

RandomPhaseShifter(d)
abs(RandomPhaseShifter(d).U[1,1])




@show
# ###### unitary of Motes et al. ######
# TO BE DEBUGGED
# U = zeros(Complex, (n+1,n+1))
# reflectivities = 0.3 .* ones(n+1)
# reflectivities[1] = 0
# reflectivities[n+1] = 1
#
# bs = BeamSplitter.(reflectivities)
#
# for i in 1:n
#     for j in 1:n
#
#         if i>j+1
#             U[i,j] = 0
#         elseif i == j+1
#             U[i,j] = bs[i].U[2,1]
#         else
#
#             p = 1 # intermediary product in the expression
#
#             if length((i+1):j) != 0
#                 p = prod([bs[k].U[1,2] for k in (i+1):j])
#             end
#             U[i,j] = bs[i].U[2,2] * bs[j+1].U[1,1] * p
#         end
#     end
# end
#
# U
