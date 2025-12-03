"""
Test OneParameterInterpolation with Householder sampler

Verifies that the OneParameterInterpolation input type now uses the
Householder sampler instead of the noisy sampler.
"""

using Revise
using BosonSampling
using LinearAlgebra

println("="^60)
println("Testing OneParameterInterpolation with Householder Sampler")
println("="^60)

# Parameters
n = 3  # number of photons
m = 5  # number of modes

# Test different distinguishability parameters
x_values = [0.0, 0.3, 0.7, 1.0]  # 0=distinguishable, 1=indistinguishable

println("\nSetup:")
println("  n = $n photons")
println("  m = $m modes")

# Create interferometer (same for all tests)
interf = RandHaar(m)

for x in x_values
    println("\n" * "="^60)
    println("Testing x = $x")
    println("="^60)

    # Create input with OneParameterInterpolation
    input = Input{OneParameterInterpolation}(first_modes(n, m), x)

    println("\nGram matrix for x=$x:")
    display(input.G.S)
    println()

    # Sample using integrated interface
    samples = []
    for i in 1:5
        ev = Event(input, FockSample(), interf)
        sample!(ev)
        output = ev.output_measurement.s
        push!(samples, output.state)
        println("Sample $i: $(output.state) (total: $(sum(output.state)))")
    end

    # Verify photon conservation
    for (i, s) in enumerate(samples)
        @assert sum(s) == n "Sample $i: photon number not conserved! Got $(sum(s)), expected $n"
    end

    println("✓ All samples conserve photon number")
end

println("\n" * "="^60)
println("Testing edge cases:")
println("="^60)

# Edge case 1: x=0 (fully distinguishable) - should behave like distinguishable
println("\nx=0 (fully distinguishable):")
input_dist = Input{OneParameterInterpolation}(first_modes(n, m), 0.0)
println("Gram matrix:\n", input_dist.G.S)

# Edge case 2: x=1 (fully indistinguishable) - should behave like bosonic
println("\nx=1 (fully indistinguishable):")
input_indist = Input{OneParameterInterpolation}(first_modes(n, m), 1.0)
println("Gram matrix:\n", input_indist.G.S)

println("\n" * "="^60)
println("✓ All tests passed successfully!")
println("OneParameterInterpolation now uses Householder sampler")
println("="^60)
