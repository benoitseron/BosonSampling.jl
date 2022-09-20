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

n = 1
m = 1
i = Input{Bosonic}(first_modes(n,m))
o = FockDetection(first_modes(n,m))
η_loss = 0.9

circuit = LossyCircuit(1)
interf = LossyLine(η_loss)
target_modes = [1]

add_element_lossy!(circuit, interf, target_modes = target_modes)

typeof([1, 1])

length(target_modes) != interf.m_real

isa(interf, LossyInterferometer)

if !(isa(interf, LossyInterferometer))
    # convert to a LossyInterferometer any lossless element just for size requirements and consistency
    @show "..."
    interf = to_lossy(interf)

end

target_modes = ModeList([1,2],4)

convert(ModeOccupation, target_modes)

# function to_lossy(target_modes)

function lossy_line_example(transmission_amplitude_loss)

    interf = Circuit(2)
    add_element!(interf, LossyLine(transmission_amplitude_loss), target_modes = [1,2])

    interf
end

transmission_amplitude_loss_array = 0:0.1:1
output_proba = []

for transmission in transmission_amplitude_loss_array
    ev = Event(i,o, lossy_line_example(transmission))
    @show compute_probability!(ev)
    push!(output_proba, ev.proba_params.probability)
end

plot(transmission_amplitude_loss_array, output_proba)
ylabel!("p no lost")
xlabel!("transmission amplitude")
