"""
Convergence Test for Clifford Algorithm

Test that TVD converges to zero as we increase sample size.
Start with small circuits (m=3) and show convergence.
"""

using BosonSampling
using Random
using Plots

include("clifford_paper_implementation.jl")
include("compare_full_distribution.jl")

Random.seed!(42)

"""
    convergence_test(n, m, sample_sizes)

Test convergence of TVD vs sample size for fixed system.
"""
function convergence_test(n, m, sample_sizes)
    println("="^60)
    println("CONVERGENCE TEST: n=$n, m=$m")
    println("="^60)
    
    # Fixed setup for consistency
    input_state = Input{Bosonic}(first_modes(n, m))
    interf = RandHaar(m)  # Fixed random interferometer
    
    # Compute exact theoretical distribution once
    println("Computing exact theoretical distribution...")
    theoretical = compute_theoretical_distribution(input_state, interf)
    
    # Test convergence for different sample sizes
    tvds = Float64[]
    
    for n_samples in sample_sizes
        println("\nTesting with $n_samples samples...")
        
        # Generate samples
        sample_counts = Dict{Vector{Int}, Int}()
        for i in 1:n_samples
            sample = clifford_sampler_paper(input_state, interf)
            sample_counts[sample] = get(sample_counts, sample, 0) + 1
            
            if i % max(1, n_samples ÷ 10) == 0
                print(".")
            end
        end
        print(" done!")
        
        # Convert to probabilities
        sampled_dist = Dict{Vector{Int}, Float64}()
        for (config, count) in sample_counts
            sampled_dist[config] = count / n_samples
        end
        
        # Compute TVD
        tvd = 0.0
        all_configs = unique([keys(theoretical)..., keys(sampled_dist)...])
        for config in all_configs
            theo_prob = get(theoretical, config, 0.0)
            samp_prob = get(sampled_dist, config, 0.0)
            tvd += abs(theo_prob - samp_prob)
        end
        tvd /= 2
        
        push!(tvds, tvd)
        println(" TVD = $(round(tvd, digits=4))")
    end
    
    return sample_sizes, tvds
end

"""
    test_multiple_systems()

Test convergence for multiple small systems.
"""
function test_multiple_systems()
    println("TESTING CONVERGENCE FOR MULTIPLE SMALL SYSTEMS")
    println("="^60)
    
    # Sample sizes to test
    sample_sizes = [500, 1000, 2000, 5000, 10000, 20000, 50000]
    
    systems = [
        (2, 3, "2 photons, 3 modes"),
        (2, 4, "2 photons, 4 modes"), 
        (3, 4, "3 photons, 4 modes")
    ]
    
    results = Dict()
    
    for (n, m, description) in systems
        println("\n" * "="^40)
        println("SYSTEM: $description")
        println("="^40)
        
        sizes, tvds = convergence_test(n, m, sample_sizes)
        results[(n, m)] = (sizes, tvds)
        
        # Show trend
        println("\nConvergence trend:")
        println("Samples | TVD")
        println("-" * 20)
        for (size, tvd) in zip(sizes, tvds)
            println("$(rpad(size, 7)) | $(round(tvd, digits=4))")
        end
        
        # Check if converging
        if length(tvds) >= 3
            recent_trend = tvds[end] - tvds[end-2]
            if recent_trend < 0
                println("✓ TVD is decreasing (good trend)")
            else
                println("⚠ TVD not clearly decreasing")
            end
        end
    end
    
    return results
end

"""
    detailed_convergence_small_system()

Very detailed convergence test on smallest system.
"""
function detailed_convergence_small_system()
    println("\n" * "="^60)
    println("DETAILED CONVERGENCE: 2 photons, 3 modes")
    println("="^60)
    
    n, m = 2, 3
    input_state = Input{Bosonic}(first_modes(n, m))
    interf = RandHaar(m)
    
    # Very detailed sample sizes
    sample_sizes = [100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000]
    
    println("Testing convergence with sample sizes: $sample_sizes")
    
    # Exact distribution
    theoretical = compute_theoretical_distribution(input_state, interf)
    n_configs = length(theoretical)
    println("System has $n_configs possible configurations")
    
    # Show theoretical distribution
    println("\nTheoretical distribution:")
    sorted_theo = sort(collect(theoretical), by=x->x[2], rev=true)
    for (config, prob) in sorted_theo
        println("  $config: $(round(prob, digits=4))")
    end
    
    println("\nTesting convergence...")
    tvds = Float64[]
    max_errors = Float64[]
    
    for n_samples in sample_sizes
        print("$n_samples samples...")
        
        # Generate samples
        sample_counts = Dict{Vector{Int}, Int}()
        for i in 1:n_samples
            sample = clifford_sampler_paper(input_state, interf)
            sample_counts[sample] = get(sample_counts, sample, 0) + 1
        end
        
        # Convert to probabilities
        sampled_dist = Dict{Vector{Int}, Float64}()
        for (config, count) in sample_counts
            sampled_dist[config] = count / n_samples
        end
        
        # Compute TVD and max error
        tvd = 0.0
        max_error = 0.0
        for config in keys(theoretical)
            theo_prob = theoretical[config]
            samp_prob = get(sampled_dist, config, 0.0)
            error = abs(theo_prob - samp_prob)
            tvd += error
            max_error = max(max_error, error)
        end
        tvd /= 2
        
        push!(tvds, tvd)
        push!(max_errors, max_error)
        
        println(" TVD=$(round(tvd, digits=4)), MaxErr=$(round(max_error, digits=4))")
    end
    
    # Summary
    println("\nCONVERGENCE SUMMARY:")
    println("Samples   | TVD    | Max Error | Improvement")
    println("-" * 45)
    for i in 1:length(sample_sizes)
        if i == 1
            improvement = "---"
        else
            improvement = "$(round((tvds[i-1] - tvds[i])/tvds[i-1] * 100, digits=1))%"
        end
        println("$(rpad(sample_sizes[i], 9)) | $(rpad(round(tvds[i], digits=4), 6)) | $(rpad(round(max_errors[i], digits=4), 9)) | $improvement")
    end
    
    # Check theoretical limit
    expected_std = sqrt(1/(4*sample_sizes[end]))  # Rough estimate for TVD standard error
    println("\nFor $(sample_sizes[end]) samples:")
    println("Final TVD: $(round(tvds[end], digits=4))")
    println("Expected statistical error: ~$(round(expected_std, digits=4))")
    
    if tvds[end] < 2*expected_std
        println("✓ TVD is consistent with statistical noise")
    else
        println("⚠ TVD still above statistical noise level")
    end
    
    return sample_sizes, tvds, max_errors
end

# Run the tests
println("Starting convergence tests...")

# Test 1: Multiple systems with moderate sample sizes
results = test_multiple_systems()

# Test 2: Detailed convergence on smallest system
detailed_results = detailed_convergence_small_system()

println("\n" * "="^60)
println("ALL CONVERGENCE TESTS COMPLETED")
println("="^60)