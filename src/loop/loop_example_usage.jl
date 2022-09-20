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

    add_element!(circuit, interf, target_modes = target_modes)
    circuit

end

lossy_line_example(0.9)

transmission_amplitude_loss_array = 0:0.1:1
output_proba = []

########### need conversion in to_lossy for Inputs etc at the construction of Events

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

### 2d examples ###

n = 2
m = 2
i = Input{Bosonic}(first_modes(n,m))
o = FockDetection(first_modes(n,m))
η_loss = 0.9

circuit = LossyCircuit(m)
interf = LossyBeamSplitter(1/sqrt(2), η_loss)
target_modes = [1,2]

add_element!(circuit, interf, target_modes = target_modes)
