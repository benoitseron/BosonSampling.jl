"""
Clifford Sampler Playground

Simple workspace to experiment with the Clifford sampler and validate it using Bayesian method.
"""

using BosonSampling
using Random

include("clifford_final.jl")

#%% Setup
Random.seed!(42)

n = 3  # photons
m = 6  # modes
n_samples = 1000

#%% Generate samples
input_state = Input{Bosonic}(first_modes(n, m))
interf = RandHaar(m)

samples = []
for i in 1:n_samples
    sample = clifford_sampler_input_final(input_state, interf)
    push!(samples, sample)
end

println("Generated $n_samples samples")
println("Example: $(samples[1])")

#%% Bayesian validation
events = []
for sample in samples
    # Convert occupancy vector to ModeOccupation  
    mode_occ = ModeOccupation(sample)
    output = FockDetection(mode_occ)
    ev = Event(input_state, output, interf)
    push!(events, ev)
end

p_bosonic = HypothesisFunction(BosonSampling.p_B)
p_distinguishable = HypothesisFunction(BosonSampling.p_D)

certifier = Bayesian(events, p_bosonic, p_distinguishable)
certify!(certifier)

confidence = certifier.confidence
println("Bayesian confidence: $(round(confidence, digits=3))")
println(confidence > 0.8 ? "✓ Working well" : "⚠ Check implementation")