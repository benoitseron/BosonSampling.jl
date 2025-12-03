"""
Test script for integrated Householder-based partial distinguishability sampler

This tests the sample!() interface with UserDefinedGramMatrix inputs.
"""

using Revise
using BosonSampling
using LinearAlgebra

println("="^60)
println("Testing Householder Sampler Integration")
println("="^60)

# Parameters
n = 3  # number of photons
m = 5  # number of modes
r = 2  # rank of Gram matrix

println("\nSetup:")
println("  n = $n photons")
println("  m = $m modes")
println("  r = $r (rank of Gram matrix)")

# Generate random Gram matrix
S = rand_gram_matrix_from_orthonormal_basis(n, r)

println("\nGram matrix S:")
display(S)
println()

# Create input state with partial distinguishability
input = Input{UserDefinedGramMatrix}(first_modes(n, m), S)

# Create interferometer
interf = RandHaar(m)

# Create event
output = FockSample()
ev = Event(input, output, interf)

println("\n" * "="^60)
println("Sampling using sample!() with UserDefinedGramMatrix")
println("="^60)

# Sample using the integrated interface
sample!(ev)

# Get result
sampled_output = ev.output_measurement.s

println("\nSampled output:")
println("  ", sampled_output)
println("  Total photons: ", sum(sampled_output.state))

# Run multiple samples to verify it works
println("\n" * "="^60)
println("Running 10 samples to verify consistency")
println("="^60)

samples = []
for i in 1:10
    ev_test = Event(input, FockSample(), interf)
    sample!(ev_test)
    push!(samples, ev_test.output_measurement.s.state)
    println("Sample $i: ", ev_test.output_measurement.s.state, " (total: ", sum(ev_test.output_measurement.s.state), ")")
end

println("\n" * "="^60)
println("✓ Integration test completed successfully!")
println("="^60)
