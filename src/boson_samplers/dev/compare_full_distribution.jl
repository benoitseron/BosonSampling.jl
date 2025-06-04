"""
Compare Clifford sampler against the full theoretical boson sampling distribution

This computes the complete theoretical probability distribution for all possible
outputs and compares it with the empirical distribution from sampling.
"""

using BosonSampling
using Random
using LinearAlgebra

include("clifford_exact.jl")

Random.seed!(42)

"""
    generate_all_configurations(n, m)

Generate all possible n-photon output configurations in m modes.
"""
function generate_all_configurations(n, m)
    configs = Vector{Int}[]
    
    function generate_recursive(remaining_photons, start_mode, current_config)
        if start_mode > m
            if remaining_photons == 0
                push!(configs, copy(current_config))
            end
            return
        end
        
        for photons_in_mode in 0:remaining_photons
            new_config = copy(current_config)
            push!(new_config, photons_in_mode)
            generate_recursive(remaining_photons - photons_in_mode, start_mode + 1, new_config)
        end
    end
    
    generate_recursive(n, 1, Int[])
    return configs
end

"""
    compute_theoretical_distribution(input_state, interf)

Compute the full theoretical boson sampling distribution.
"""
function compute_theoretical_distribution(input_state, interf)
    n = input_state.n
    m = input_state.m
    
    println("Computing full theoretical distribution...")
    
    # Generate all possible configurations
    all_configs = generate_all_configurations(n, m)
    println("Total configurations: $(length(all_configs))")
    
    # Compute probability for each configuration
    theoretical_dist = Dict{Vector{Int}, Float64}()
    total_prob = 0.0
    
    for (i, config) in enumerate(all_configs)
        # Create event for this configuration
        mode_occ = ModeOccupation(config)
        output = FockDetection(mode_occ)
        ev = Event(input_state, output, interf)
        
        # Compute exact probability using full distribution method
        BosonSampling.compute_probability!(ev)
        prob = real(ev.proba_params.probability)
        
        theoretical_dist[config] = prob
        total_prob += prob
        
        if i % 10 == 0
            print(".")
        end
    end
    println(" done!")
    
    println("Total probability mass: $(round(total_prob, digits=6))")
    
    # Normalize (should be ≈ 1 already)
    if abs(total_prob - 1.0) > 1e-6
        println("Warning: Total probability = $total_prob (not 1.0)")
        for config in keys(theoretical_dist)
            theoretical_dist[config] /= total_prob
        end
    end
    
    return theoretical_dist
end

"""
    compute_sampled_distribution(input_state, interf, n_samples)

Generate empirical distribution from sampling.
"""
function compute_sampled_distribution(input_state, interf, n_samples)
    println("Generating $n_samples samples...")
    
    sample_counts = Dict{Vector{Int}, Int}()
    
    for i in 1:n_samples
        sample = clifford_sampler_exact(input_state, interf)
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

"""
    compare_distributions(theoretical, sampled)

Compare theoretical and sampled distributions.
"""
function compare_distributions(theoretical, sampled)
    println("\nDISTRIBUTION COMPARISON:")
    println("="^70)
    
    # Get all configurations that appear in either distribution
    all_configs = unique([keys(theoretical)..., keys(sampled)...])
    
    # Sort by theoretical probability (highest first)
    sorted_configs = sort(all_configs, by=x->get(theoretical, x, 0.0), rev=true)
    
    println("Configuration       | Theoretical | Sampled   | Difference | Ratio")
    println(repeat("-", 70))
    
    total_variation = 0.0
    max_error = 0.0
    significant_configs = 0
    
    for config in sorted_configs
        theo_prob = get(theoretical, config, 0.0)
        samp_prob = get(sampled, config, 0.0)
        difference = samp_prob - theo_prob
        ratio = samp_prob > 0 && theo_prob > 0 ? samp_prob / theo_prob : NaN
        
        # Only show configurations with non-negligible probability
        if theo_prob > 1e-4 || samp_prob > 1e-4
            significant_configs += 1
            ratio_str = isnan(ratio) ? "---" : "$(round(ratio, digits=2))"
            println("$(rpad(string(config), 19)) | $(rpad(round(theo_prob, digits=4), 11)) | $(rpad(round(samp_prob, digits=4), 9)) | $(rpad(round(difference, sigdigits=3), 10)) | $ratio_str")
        end
        
        total_variation += abs(difference)
        max_error = max(max_error, abs(difference))
    end
    
    total_variation /= 2  # Standard TVD normalization
    
    println(repeat("-", 70))
    println("SUMMARY:")
    println("Configurations shown: $significant_configs")
    println("Total Variation Distance: $(round(total_variation, digits=4))")
    println("Maximum single error: $(round(max_error, digits=4))")
    
    # Quality assessment
    if total_variation < 0.01
        println("✓ EXCELLENT: Distributions match very well")
    elseif total_variation < 0.05
        println("✓ GOOD: Distributions match reasonably well")
    elseif total_variation < 0.1
        println("⚠ ACCEPTABLE: Some discrepancies but roughly correct")
    elseif total_variation < 0.2
        println("⚠ POOR: Significant discrepancies")
    else
        println("✗ VERY POOR: Major distribution mismatch")
    end
    
    return total_variation
end

"""
    full_validation_test(n, m, n_samples=10000)

Complete validation test comparing full distributions.
"""
function full_validation_test(n, m, n_samples=10000)
    println("FULL DISTRIBUTION VALIDATION")
    println("="^50)
    println("System: $n photons, $m modes")
    println("Samples: $n_samples")
    println()
    
    # Create fixed input and interferometer
    input_state = Input{Bosonic}(first_modes(n, m))
    interf = RandHaar(m)
    
    # Compute theoretical distribution
    theoretical = compute_theoretical_distribution(input_state, interf)
    
    # Generate sampled distribution  
    sampled = compute_sampled_distribution(input_state, interf, n_samples)
    
    # Compare
    tvd = compare_distributions(theoretical, sampled)
    
    println("="^50)
    return theoretical, sampled, tvd
end

# Test with small system first
println("Testing 2 photons, 3 modes...")
theo1, samp1, tvd1 = full_validation_test(2, 3, 5000)

if tvd1 < 0.1
    println("\nSmall system acceptable, testing larger...")
    theo2, samp2, tvd2 = full_validation_test(3, 4, 10000)
else
    println("\nSmall system already shows issues.")
end