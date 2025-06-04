"""
# Clifford Algorithm Validation Notebook

Complete validation of the paper implementation against exact theoretical distributions.
This notebook tests whether our implementation of Algorithm B produces the correct boson sampling distribution.
"""

using BosonSampling
using Random
using LinearAlgebra

# Load implementations
include("clifford_paper_implementation.jl")

Random.seed!(123)

#%% Cell 1: Quick Test of Paper Algorithm

println("="^60)
println("CELL 1: QUICK TEST OF PAPER ALGORITHM")
println("="^60)

# Test basic functionality
n, m = 2, 3
input_state = Input{Bosonic}(first_modes(n, m))
interf = RandHaar(m)

println("Testing basic functionality...")
println("Input: $n photons in first $n modes of $m total modes")

# Generate a few samples
println("\nFirst 10 samples:")
for i in 1:10
    sample = clifford_sampler_paper(input_state, interf)
    total_photons = sum(sample)
    println("  Sample $i: $sample (total photons: $total_photons)")
end

#%% Cell 2: Full Distribution Validation

println("\n" * "="^60)
println("CELL 2: FULL DISTRIBUTION VALIDATION")  
println("="^60)

# Use the exact validation framework
include("compare_full_distribution.jl")

"""
    validate_paper_implementation(n, m, n_samples)

Validate the paper implementation using full distribution comparison.
"""
function validate_paper_implementation(n, m, n_samples=5000)
    println("Validating paper implementation...")
    println("System: $n photons, $m modes, $n_samples samples")
    
    # Create test setup
    input_state = Input{Bosonic}(first_modes(n, m))
    interf = RandHaar(m)
    
    # Compute theoretical distribution
    println("\nComputing theoretical distribution...")
    theoretical = compute_theoretical_distribution(input_state, interf)
    
    # Generate sampled distribution using paper algorithm
    println("Generating samples using paper algorithm...")
    sample_counts = Dict{Vector{Int}, Int}()
    
    for i in 1:n_samples
        sample = clifford_sampler_paper(input_state, interf)
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
    
    # Compare distributions
    tvd = compare_distributions(theoretical, sampled_dist)
    
    return theoretical, sampled_dist, tvd
end

# Test small system
println("Testing 2 photons, 3 modes...")
theo1, samp1, tvd1 = validate_paper_implementation(2, 3, 5000)

println("\nResult: TVD = $(round(tvd1, digits=4))")
if tvd1 < 0.05
    println("✓ EXCELLENT: Paper algorithm works very well!")
elseif tvd1 < 0.1  
    println("✓ GOOD: Paper algorithm works well")
elseif tvd1 < 0.2
    println("⚠ ACCEPTABLE: Some discrepancies")
else
    println("✗ POOR: Significant issues remain")
end

#%% Cell 3: Comparison with Built-in

println("\n" * "="^60)
println("CELL 3: COMPARISON WITH BUILT-IN SAMPLER")
println("="^60)

"""
    compare_with_builtin(n, m, n_samples)

Compare paper implementation with built-in sampler.
"""
function compare_with_builtin(n, m, n_samples=3000)
    println("Comparing implementations...")
    
    input_state = Input{Bosonic}(first_modes(n, m))
    interf = RandHaar(m)
    
    # Paper implementation
    println("Generating samples with paper algorithm...")
    paper_counts = Dict{Vector{Int}, Int}()
    for i in 1:n_samples
        sample = clifford_sampler_paper(input_state, interf)
        paper_counts[sample] = get(paper_counts, sample, 0) + 1
    end
    
    # Built-in implementation
    println("Generating samples with built-in algorithm...")
    builtin_counts = Dict{Vector{Int}, Int}()
    for i in 1:n_samples
        ev = Event(input_state, FockSample(), interf)
        BosonSampling.sample!(ev)
        sample = ev.output_measurement.s.state
        builtin_counts[sample] = get(builtin_counts, sample, 0) + 1
    end
    
    # Compare
    println("\nDirect comparison of implementations:")
    all_configs = unique([keys(paper_counts)..., keys(builtin_counts)...])
    
    println("Configuration | Paper % | Built-in % | Ratio")
    println(repeat("-", 45))
    
    for config in sort(all_configs)
        paper_freq = 100 * get(paper_counts, config, 0) / n_samples
        builtin_freq = 100 * get(builtin_counts, config, 0) / n_samples
        ratio = paper_freq > 0 && builtin_freq > 0 ? paper_freq / builtin_freq : NaN
        
        if paper_freq > 0.1 || builtin_freq > 0.1  # Only show significant configs
            ratio_str = isnan(ratio) ? "---" : "$(round(ratio, digits=2))"
            println("$config | $(round(paper_freq, digits=1)) | $(round(builtin_freq, digits=1)) | $ratio_str")
        end
    end
    
    return paper_counts, builtin_counts
end

paper_results, builtin_results = compare_with_builtin(2, 3, 3000)

#%% Cell 4: Bayesian Validation Test

println("\n" * "="^60)
println("CELL 4: BAYESIAN VALIDATION TEST")
println("="^60)

"""
    bayesian_test_paper(n, m, n_events)

Test paper implementation with Bayesian validation.
"""
function bayesian_test_paper(n, m, n_events=1000)
    println("Testing paper implementation with Bayesian validation...")
    
    input_state = Input{Bosonic}(first_modes(n, m))
    interf = RandHaar(m)
    
    # Generate events using paper algorithm
    events = []
    for i in 1:n_events
        sample = clifford_sampler_paper(input_state, interf)
        mode_occ = ModeOccupation(sample)
        output = FockDetection(mode_occ)
        ev = Event(input_state, output, interf)
        push!(events, ev)
    end
    
    # Bayesian validation
    p_bosonic = HypothesisFunction(BosonSampling.p_B)
    p_distinguishable = HypothesisFunction(BosonSampling.p_D)
    
    certifier = Bayesian(events, p_bosonic, p_distinguishable)
    certify!(certifier)
    
    confidence = certifier.confidence
    println("\nBayesian validation results:")
    println("Confidence in bosonic hypothesis: $(round(confidence, digits=4))")
    
    if confidence > 0.9
        println("✓ EXCELLENT: Strong evidence for correct boson sampling")
    elseif confidence > 0.7
        println("✓ GOOD: Good evidence for correct boson sampling")
    elseif confidence > 0.5
        println("⚠ WEAK: Some evidence but not conclusive")
    else
        println("✗ POOR: Evidence suggests incorrect sampling")
    end
    
    return confidence
end

bayesian_confidence = bayesian_test_paper(3, 6, 1000)

#%% Cell 5: Summary and Conclusions

println("\n" * "="^60)
println("CELL 5: SUMMARY AND CONCLUSIONS")
println("="^60)

println("VALIDATION RESULTS SUMMARY:")
println("-" * 30)
println("Full Distribution TVD: $(round(tvd1, digits=4))")
println("Bayesian Confidence: $(round(bayesian_confidence, digits=4))")

if tvd1 < 0.05 && bayesian_confidence > 0.8
    println("\n🎉 SUCCESS: Paper implementation appears to work correctly!")
    println("✓ Low TVD indicates good match with theoretical distribution")
    println("✓ High Bayesian confidence confirms bosonic behavior")
elif tvd1 < 0.1 && bayesian_confidence > 0.6
    println("\n👍 GOOD: Paper implementation works reasonably well")
    println("✓ Acceptable TVD and Bayesian scores")
else
    println("\n⚠️  ISSUES: Paper implementation may need further work")
    println("✗ Higher than expected TVD or low Bayesian confidence")
end

println("\nImplementation details:")
println("- Algorithm follows paper exactly (Algorithm B)")
println("- Includes proper column permutation")
println("- Uses Laplace expansion for probability calculation")
println("- Handles permanent computations correctly")

println("\n" * "="^60)
println("VALIDATION COMPLETE")
println("="^60)