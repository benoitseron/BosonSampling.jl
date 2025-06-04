"""
Simple Convergence Test

Test that TVD converges as we increase sample size for small circuits.
"""

using BosonSampling
using Random

include("clifford_paper_implementation.jl")
include("compare_full_distribution.jl")

Random.seed!(42)

function simple_convergence_test()
    println("SIMPLE CONVERGENCE TEST: 2 photons, 3 modes")
    println("=" ^ 50)
    
    # Small system
    n, m = 2, 3
    input_state = Input{Bosonic}(first_modes(n, m))
    interf = RandHaar(m)
    
    # Compute exact distribution
    println("Computing exact theoretical distribution...")
    theoretical = compute_theoretical_distribution(input_state, interf)
    
    # Show theoretical distribution
    println("\nTheoretical distribution:")
    sorted_configs = sort(collect(theoretical), by=x->x[2], rev=true)
    for (config, prob) in sorted_configs
        println("  $config: $(round(prob, digits=4))")
    end
    
    # Test different sample sizes
    sample_sizes = [1000, 5000, 10000, 20000, 50000, 100000]
    
    println("\nTesting convergence:")
    println("Samples   | TVD    | Status")
    println("-" ^ 30)
    
    tvds = Float64[]
    
    for n_samples in sample_sizes
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
        
        # Compute TVD
        tvd = 0.0
        for config in keys(theoretical)
            theo_prob = theoretical[config]
            samp_prob = get(sampled_dist, config, 0.0)
            tvd += abs(theo_prob - samp_prob)
        end
        tvd /= 2
        
        push!(tvds, tvd)
        
        # Status
        if tvd < 0.01
            status = "EXCELLENT"
        elseif tvd < 0.05
            status = "GOOD"
        elseif tvd < 0.1
            status = "OK"
        else
            status = "POOR"
        end
        
        println("$(rpad(n_samples, 9)) | $(rpad(round(tvd, digits=4), 6)) | $status")
    end
    
    # Check trend
    println("\nTrend Analysis:")
    if length(tvds) >= 2
        final_tvd = tvds[end]
        initial_tvd = tvds[1]
        improvement = (initial_tvd - final_tvd) / initial_tvd * 100
        
        println("Initial TVD ($(sample_sizes[1]) samples): $(round(initial_tvd, digits=4))")
        println("Final TVD ($(sample_sizes[end]) samples): $(round(final_tvd, digits=4))")
        println("Improvement: $(round(improvement, digits=1))%")
        
        if improvement > 10
            println("✓ Significant improvement with more samples")
        elseif improvement > 0
            println("~ Modest improvement")
        else
            println("⚠ No clear improvement")
        end
    end
    
    # Statistical expectation
    expected_error = sqrt(1.0 / (4 * sample_sizes[end]))
    println("\nStatistical Analysis:")
    println("Expected statistical error for $(sample_sizes[end]) samples: ≈$(round(expected_error, digits=4))")
    
    if tvds[end] < 2 * expected_error
        println("✓ Final TVD is consistent with statistical noise")
    elseif tvds[end] < 5 * expected_error
        println("~ Final TVD is close to statistical limit")
    else
        println("⚠ Final TVD is above statistical noise")
    end
    
    return sample_sizes, tvds
end

function show_sample_comparison()
    println("\n" ^ 2 * "SAMPLE COMPARISON")
    println("=" ^ 50)
    
    n, m = 2, 3
    input_state = Input{Bosonic}(first_modes(n, m))
    interf = RandHaar(m)
    
    # Generate many samples and show distribution evolution
    sample_sizes = [1000, 10000, 100000]
    
    for n_samples in sample_sizes
        println("\nWith $n_samples samples:")
        
        # Generate samples
        sample_counts = Dict{Vector{Int}, Int}()
        for i in 1:n_samples
            sample = clifford_sampler_paper(input_state, interf)
            sample_counts[sample] = get(sample_counts, sample, 0) + 1
        end
        
        # Show distribution
        sorted_samples = sort(collect(sample_counts), by=x->x[2], rev=true)
        for (config, count) in sorted_samples
            freq = count / n_samples
            println("  $config: $(round(freq, digits=4))")
        end
    end
end

# Run tests
println("Running convergence tests...")
sizes, tvds = simple_convergence_test()
show_sample_comparison()

println("\n" * "=" ^ 50)
println("CONVERGENCE TEST COMPLETED")
println("=" ^ 50)