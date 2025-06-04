"""
Validate Clifford sampler by comparing exact vs sampled distributions

This is the proper way to test: compute exact probabilities for all possible
outputs and compare with the frequency of samples from our algorithm.
"""

using BosonSampling
using Random
using LinearAlgebra

include("clifford_final.jl")

Random.seed!(123)

"""
    all_possible_outputs(n, m)

Generate all possible output configurations for n photons in m modes.
"""
function all_possible_outputs(n, m)
    configs = Vector{Int}[]
    
    function generate_config(remaining_photons, current_mode, current_config)
        if current_mode > m
            if remaining_photons == 0
                push!(configs, copy(current_config))
            end
            return
        end
        
        # Try all possible photon counts in current mode
        for count in 0:remaining_photons
            push!(current_config, count)
            generate_config(remaining_photons - count, current_mode + 1, current_config)
            pop!(current_config)
        end
    end
    
    generate_config(n, 1, Int[])
    return configs
end

"""
    compute_exact_probabilities(input_state, interf)

Compute exact probability for each possible output configuration.
"""
function compute_exact_probabilities(input_state, interf)
    n = input_state.n
    m = input_state.m
    
    # Get all possible outputs
    all_outputs = all_possible_outputs(n, m)
    
    # Compute exact probability for each
    exact_probs = Dict{Vector{Int}, Float64}()
    
    for output_config in all_outputs
        # Create event
        output_mode_occ = ModeOccupation(output_config)
        output = FockDetection(output_mode_occ)
        ev = Event(input_state, output, interf)
        
        # Compute exact probability
        BosonSampling.compute_probability!(ev)
        prob = real(ev.proba_params.probability)
        
        exact_probs[output_config] = prob
    end
    
    return exact_probs
end

"""
    validate_clifford_sampler(n, m, n_samples=10000)

Main validation function comparing exact vs sampled distributions.
"""
function validate_clifford_sampler(n, m, n_samples=10000)
    println("="^60)
    println("VALIDATION: EXACT vs SAMPLED DISTRIBUTIONS")
    println("="^60)
    println("System: n=$n photons, m=$m modes")
    println("Samples: $n_samples")
    
    # Create input and fixed interferometer
    input_state = Input{Bosonic}(first_modes(n, m))
    interf = RandHaar(m)  # Use fixed random interferometer
    
    println("\nComputing exact probabilities...")
    exact_probs = compute_exact_probabilities(input_state, interf)
    n_configs = length(exact_probs)
    println("Found $n_configs possible output configurations")
    
    println("\nGenerating samples...")
    sample_counts = Dict{Vector{Int}, Int}()
    
    for i in 1:n_samples
        sample = clifford_sampler_input_final(input_state, interf)
        sample_counts[sample] = get(sample_counts, sample, 0) + 1
        
        if i % (n_samples ÷ 10) == 0
            print(".")
        end
    end
    println(" done!")
    
    # Convert counts to probabilities
    sampled_probs = Dict{Vector{Int}, Float64}()
    for (config, count) in sample_counts
        sampled_probs[config] = count / n_samples
    end
    
    # Compare distributions
    println("\nCOMPARISON:")
    println("Configuration | Exact Prob | Sampled Prob | Ratio")
    println(repeat("-", 55))
    
    total_variation = 0.0
    max_ratio_error = 0.0
    
    # Sort by exact probability (highest first)
    sorted_configs = sort(collect(exact_probs), by=x->x[2], rev=true)
    
    for (config, exact_prob) in sorted_configs[1:min(10, length(sorted_configs))]
        sampled_prob = get(sampled_probs, config, 0.0)
        
        if exact_prob > 1e-10  # Only show non-negligible probabilities
            ratio = sampled_prob > 0 ? sampled_prob / exact_prob : 0.0
            ratio_error = abs(ratio - 1.0)
            max_ratio_error = max(max_ratio_error, ratio_error)
            
            println("$config | $(round(exact_prob, digits=4)) | $(round(sampled_prob, digits=4)) | $(round(ratio, digits=2))")
        end
        
        total_variation += abs(exact_prob - sampled_prob)
    end
    
    # Check for missing configurations
    missing_configs = 0
    for (config, exact_prob) in exact_probs
        if exact_prob > 1e-6 && !haskey(sampled_probs, config)
            missing_configs += 1
        end
    end
    
    # Check for spurious configurations  
    spurious_configs = 0
    for config in keys(sampled_probs)
        if !haskey(exact_probs, config)
            spurious_configs += 1
        end
    end
    
    # Summary
    total_variation /= 2  # TVD normalization
    
    println("\nSUMMARY:")
    println("Total Variation Distance: $(round(total_variation, digits=4))")
    println("Max ratio error: $(round(max_ratio_error, digits=4))")
    println("Missing configs (prob > 1e-6): $missing_configs")
    println("Spurious configs: $spurious_configs")
    
    # Verdict
    if total_variation < 0.05 && max_ratio_error < 0.2
        println("✓ EXCELLENT: Sampler appears to work correctly")
    elseif total_variation < 0.1 && max_ratio_error < 0.5
        println("✓ GOOD: Sampler works reasonably well")
    elseif total_variation < 0.2
        println("⚠ ACCEPTABLE: Some discrepancies but roughly correct")
    else
        println("✗ POOR: Significant discrepancies detected")
    end
    
    println("="^60)
    
    return exact_probs, sampled_probs, total_variation
end

# Test with very small system first
println("Testing with tiny system...")
exact, sampled, tvd = validate_clifford_sampler(2, 3, 5000)

# Test slightly larger if small one works
if tvd < 0.1
    println("\nSmall system looks good, testing larger...")
    validate_clifford_sampler(3, 4, 10000)
end