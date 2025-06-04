"""
Test the built-in Clifford sampler against the full distribution
"""

using BosonSampling
using Random

Random.seed!(42)

# Use the same comparison functions
include("compare_full_distribution.jl")

"""
    compute_builtin_sampled_distribution(input_state, interf, n_samples)

Generate empirical distribution using the built-in sampler.
"""
function compute_builtin_sampled_distribution(input_state, interf, n_samples)
    println("Generating $n_samples samples using BUILT-IN sampler...")
    
    sample_counts = Dict{Vector{Int}, Int}()
    
    for i in 1:n_samples
        # Use built-in sampler
        ev = Event(input_state, FockSample(), interf)
        BosonSampling.sample!(ev)
        sample = ev.output_measurement.s.state
        
        sample_counts[sample] = get(sample_counts, sample, 0) + 1
        
        if i % (n_samples ÷ 20) == 0
            print(".")
        end
    end
    println(" done!")
    
    # Convert to probabilities
    sampled_dist = Dict{Vector{Int}, Float64}()
    for (config, count) in sample_counts
        sampled_dist[config] = count / n_samples
    end
    
    return sampled_dist
end

function test_builtin_distribution(n, m, n_samples=5000)
    println("TESTING BUILT-IN CLIFFORD SAMPLER")
    println("="^50)
    println("System: $n photons, $m modes")
    println("Samples: $n_samples")
    println()
    
    # Create the same setup
    input_state = Input{Bosonic}(first_modes(n, m))
    interf = RandHaar(m)
    
    # Compute theoretical distribution (same as before)
    theoretical = compute_theoretical_distribution(input_state, interf)
    
    # Generate sampled distribution using BUILT-IN
    sampled = compute_builtin_sampled_distribution(input_state, interf, n_samples)
    
    # Compare
    tvd = compare_distributions(theoretical, sampled)
    
    println("="^50)
    return tvd
end

# Test built-in sampler
println("Testing BUILT-IN Clifford sampler...")
builtin_tvd = test_builtin_distribution(2, 3, 5000)
println("\nBuilt-in TVD: $(round(builtin_tvd, digits=4))")