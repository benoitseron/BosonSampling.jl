"""
Test Corrected Clifford Algorithm

Validate that the corrected implementation fixes the m > n bias issue.
"""

using BosonSampling
using Random
include("corrected_clifford_algorithm.jl")
include("compare_full_distribution.jl")

function test_corrected_algorithm()
    println("TESTING CORRECTED CLIFFORD ALGORITHM")
    println("=" ^ 50)
    
    Random.seed!(42)
    
    # Test configurations that previously showed bias
    test_configs = [
        (2, 2, "2x2 system (should work)"),
        (2, 3, "2x3 system (previously biased)"),
        (2, 4, "2x4 system (previously biased)"),
        (3, 4, "3x4 system (previously biased)")
    ]
    
    for (n, m, description) in test_configs
        println("\nTesting: $description")
        
        input_state = Input{Bosonic}(first_modes(n, m))
        interf = RandHaar(m)
        
        # Generate samples with corrected algorithm
        n_samples = 10000
        counts = Dict{Vector{Int}, Int}()
        
        for _ in 1:n_samples
            sample = corrected_clifford_sampler(input_state, interf)
            counts[sample] = get(counts, sample, 0) + 1
        end
        
        # Compute TVD
        tvd = 0.0
        for (config, count) in counts
            emp_prob = count / n_samples
            mode_occ = ModeOccupation(config)
            output = FockDetection(mode_occ)
            ev = Event(input_state, output, interf)
            BosonSampling.compute_probability!(ev)
            exact_prob = real(ev.proba_params.probability)
            tvd += abs(emp_prob - exact_prob)
        end
        tvd /= 2
        
        status = tvd < 0.05 ? "✅ FIXED!" : "❌ Still biased"
        println("  TVD: $(round(tvd, digits=4)) $status")
        
        if tvd > 0.05
            println("  ⚠️  Algorithm still has bias - needs further debugging")
        end
    end
end

function compare_algorithms()
    println("\nCOMPARING OLD vs CORRECTED ALGORITHMS")
    println("=" ^ 50)
    
    # Test on problematic 2x3 system
    n, m = 2, 3
    Random.seed!(42)
    
    input_state = Input{Bosonic}(first_modes(n, m))
    interf = RandHaar(m)
    
    n_samples = 5000
    
    # Old algorithm (if available)
    println("Testing 2x3 system:")
    
    try
        include("clifford_paper_implementation.jl")
        
        # Old algorithm
        old_counts = Dict{Vector{Int}, Int}()
        for _ in 1:n_samples
            sample = clifford_sampler_paper(input_state, interf)
            old_counts[sample] = get(old_counts, sample, 0) + 1
        end
        
        old_tvd = 0.0
        for (config, count) in old_counts
            emp_prob = count / n_samples
            mode_occ = ModeOccupation(config)
            output = FockDetection(mode_occ)
            ev = Event(input_state, output, interf)
            BosonSampling.compute_probability!(ev)
            exact_prob = real(ev.proba_params.probability)
            old_tvd += abs(emp_prob - exact_prob)
        end
        old_tvd /= 2
        
        println("Old algorithm TVD: $(round(old_tvd, digits=4))")
        
    catch e
        println("Could not test old algorithm: $e")
    end
    
    # Reset seed for fair comparison
    Random.seed!(42)
    
    # Corrected algorithm
    corrected_counts = Dict{Vector{Int}, Int}()
    for _ in 1:n_samples
        sample = corrected_clifford_sampler(input_state, interf)
        corrected_counts[sample] = get(corrected_counts, sample, 0) + 1
    end
    
    corrected_tvd = 0.0
    for (config, count) in corrected_counts
        emp_prob = count / n_samples
        mode_occ = ModeOccupation(config)
        output = FockDetection(mode_occ)
        ev = Event(input_state, output, interf)
        BosonSampling.compute_probability!(ev)
        exact_prob = real(ev.proba_params.probability)
        corrected_tvd += abs(emp_prob - exact_prob)
    end
    corrected_tvd /= 2
    
    println("Corrected algorithm TVD: $(round(corrected_tvd, digits=4))")
    
    improvement = corrected_tvd < 0.05 ? "✅ SIGNIFICANT IMPROVEMENT!" : "❌ Still needs work"
    println("Result: $improvement")
end

function bayesian_validation_test()
    println("\nBAYESIAN VALIDATION TEST")
    println("=" ^ 50)
    
    # Test corrected algorithm with Bayesian validation
    n, m = 2, 3
    Random.seed!(123)
    
    input_state = Input{Bosonic}(first_modes(n, m))
    interf = RandHaar(m)
    
    # Generate samples
    n_samples = 2000
    samples = []
    
    for _ in 1:n_samples
        sample = corrected_clifford_sampler(input_state, interf)
        push!(samples, sample)
    end
    
    # Bayesian validation
    try
        include("../../../certification/bayesian.jl")
        
        println("Running Bayesian hypothesis test...")
        confidence = bayesian_test(samples, input_state, interf)
        println("Bayesian confidence: $(round(confidence * 100, digits=2))%")
        
        if confidence > 0.95
            println("✅ EXCELLENT: High confidence in sampler correctness")
        elseif confidence > 0.8
            println("✅ GOOD: Reasonable confidence in sampler")
        else
            println("❌ POOR: Low confidence - sampler likely has issues")
        end
        
    catch e
        println("Could not run Bayesian test: $e")
    end
end

# Run all tests
test_corrected_algorithm()
compare_algorithms()
bayesian_validation_test()

println("\n" * "=" ^ 50)
println("CORRECTED ALGORITHM TESTING COMPLETE")
println("=" ^ 50)