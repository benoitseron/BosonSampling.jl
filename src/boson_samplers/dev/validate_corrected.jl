"""
Validate Corrected Clifford Algorithm

Proper validation without using problematic constructors.
"""

using BosonSampling
using Random
include("corrected_clifford_algorithm.jl")

function compute_exact_probability(input_state, output_config, interf)
    """Compute exact probability for a given output configuration"""
    # Create output measurement
    mode_occ = ModeOccupation(output_config)
    output = FockDetection(mode_occ)
    ev = Event(input_state, output, interf)
    BosonSampling.compute_probability!(ev)
    return real(ev.proba_params.probability)
end

function test_tvd(n, m, n_samples=5000)
    """Test Total Variation Distance for given system size"""
    Random.seed!(42)
    
    input_state = Input{Bosonic}(first_modes(n, m))
    interf = RandHaar(m)
    
    println("Testing $(n)x$(m) system with $n_samples samples...")
    
    # Generate samples
    counts = Dict{Vector{Int}, Int}()
    for i in 1:n_samples
        if i % 1000 == 0
            print(".")
        end
        sample = corrected_clifford_sampler(input_state, interf)
        counts[sample] = get(counts, sample, 0) + 1
    end
    println()
    
    # Compute TVD
    tvd = 0.0
    n_configs = length(counts)
    
    for (i, (config, count)) in enumerate(counts)
        if i % 10 == 0
            print(".")
        end
        emp_prob = count / n_samples
        exact_prob = compute_exact_probability(input_state, config, interf)
        tvd += abs(emp_prob - exact_prob)
    end
    println()
    tvd /= 2
    
    println("  Configurations found: $n_configs")
    println("  TVD: $(round(tvd, digits=4))")
    
    return tvd, counts
end

function comprehensive_test()
    println("COMPREHENSIVE VALIDATION OF CORRECTED ALGORITHM")
    println("=" ^ 60)
    
    # Test configurations
    test_configs = [
        (2, 2, "Baseline: should work"),
        (2, 3, "Previously biased case"),
        (2, 4, "Higher m/n ratio"),
        (3, 3, "Larger n=m case"),
        (3, 4, "Larger m>n case")
    ]
    
    results = []
    
    for (n, m, description) in test_configs
        println("\n" * "=" ^ 40)
        println("Test: $n photons, $m modes")
        println("Description: $description")
        println("=" ^ 40)
        
        try
            tvd, counts = test_tvd(n, m, 5000)
            
            # Analyze results
            if tvd < 0.01
                status = "✅ EXCELLENT"
                color = "green"
            elseif tvd < 0.05
                status = "✅ GOOD"  
                color = "green"
            elseif tvd < 0.1
                status = "⚠️  MARGINAL"
                color = "yellow"
            else
                status = "❌ POOR"
                color = "red"
            end
            
            println("Result: $status (TVD = $(round(tvd, digits=4)))")
            
            # Show top configurations for problematic cases
            if tvd > 0.05
                println("\nTop configurations (empirical vs exact):")
                sorted_counts = sort(collect(counts), by=x->x[2], rev=true)
                for (config, count) in sorted_counts[1:min(3, length(sorted_counts))]
                    emp_prob = count / 5000
                    exact_prob = compute_exact_probability(Input{Bosonic}(first_modes(n, m)), config, RandHaar(m))
                    ratio = emp_prob / exact_prob
                    println("  $config: $(round(emp_prob, digits=3)) vs $(round(exact_prob, digits=3)) (ratio $(round(ratio, digits=2)))")
                end
            end
            
            push!(results, (n, m, tvd, status))
            
        catch e
            println("❌ ERROR: Test failed with $e")
            push!(results, (n, m, NaN, "ERROR"))
        end
    end
    
    # Summary
    println("\n" * "=" ^ 60)
    println("SUMMARY")
    println("=" ^ 60)
    
    println("System | TVD     | Status")
    println("-" ^ 30)
    for (n, m, tvd, status) in results
        tvd_str = isnan(tvd) ? "ERROR" : "$(round(tvd, digits=4))"
        println("$(n)x$(m)    | $(rpad(tvd_str, 7)) | $status")
    end
    
    # Overall assessment
    good_results = count(r -> r[3] < 0.05, filter(r -> !isnan(r[3]), results))
    total_results = length(filter(r -> !isnan(r[3]), results))
    
    println("\nOverall: $good_results/$total_results tests passed (TVD < 0.05)")
    
    if good_results == total_results
        println("🎉 COMPLETE SUCCESS: All tests passed!")
        println("   The corrected algorithm appears to fix the bias issues.")
    elseif good_results > total_results / 2
        println("✅ MOSTLY SUCCESSFUL: Most tests passed")
        println("   The algorithm is improved but may need minor adjustments.")
    else
        println("❌ NEEDS MORE WORK: Many tests still failing")
        println("   The algorithm needs further debugging.")
    end
    
    return results
end

# Run comprehensive test
results = comprehensive_test()