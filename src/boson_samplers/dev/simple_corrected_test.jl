"""
Simple Test of Corrected Clifford Algorithm
"""

using BosonSampling
using Random
include("corrected_clifford_algorithm.jl")

function simple_test()
    println("SIMPLE TEST OF CORRECTED ALGORITHM")
    println("=" ^ 40)
    
    Random.seed!(42)
    
    # Test the problematic 2x3 system
    n, m = 2, 3
    input_state = Input{Bosonic}(first_modes(n, m))
    interf = RandHaar(m)
    
    println("Testing 2 photons, 3 modes...")
    
    # Generate samples
    n_samples = 5000
    counts = Dict{Vector{Int}, Int}()
    
    for i in 1:n_samples
        if i % 1000 == 0
            println("  Generated $i samples...")
        end
        sample = corrected_clifford_sampler(input_state, interf)
        counts[sample] = get(counts, sample, 0) + 1
    end
    
    # Compute TVD
    tvd = 0.0
    println("\nComputing exact probabilities...")
    
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
    
    println("\nResults:")
    println("TVD: $(round(tvd, digits=4))")
    
    if tvd < 0.05
        println("✅ SUCCESS: Algorithm appears to be working correctly!")
        println("   (TVD < 0.05 indicates good agreement with exact probabilities)")
    else
        println("❌ FAILURE: Algorithm still has significant bias")
        println("   (TVD = $(round(tvd, digits=4)) is too high)")
        
        # Show some examples
        println("\nTop configurations (empirical vs exact):")
        sorted_counts = sort(collect(counts), by=x->x[2], rev=true)[1:min(3, length(counts))]
        for (config, count) in sorted_counts
            emp_prob = count / n_samples
            mode_occ = ModeOccupation(config)
            output = FockDetection(mode_occ)
            ev = Event(input_state, output, interf)
            BosonSampling.compute_probability!(ev)
            exact_prob = real(ev.proba_params.probability)
            ratio = emp_prob / exact_prob
            println("  $config: $(round(emp_prob, digits=3)) vs $(round(exact_prob, digits=3)) (ratio $(round(ratio, digits=2)))")
        end
    end
    
    return tvd < 0.05
end

function test_2x2_baseline()
    println("\n2x2 BASELINE TEST")
    println("=" ^ 40)
    
    Random.seed!(42)
    
    # Test 2x2 system that should work
    n, m = 2, 2
    input_state = Input{Bosonic}(first_modes(n, m))
    interf = RandHaar(m)
    
    println("Testing 2 photons, 2 modes (baseline)...")
    
    # Generate samples
    n_samples = 5000
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
    
    println("2x2 TVD: $(round(tvd, digits=4))")
    
    if tvd < 0.05
        println("✅ 2x2 system works correctly")
    else
        println("❌ Even 2x2 system has issues - fundamental problem")
    end
    
    return tvd < 0.05
end

# Run tests
baseline_ok = test_2x2_baseline()
corrected_ok = simple_test()

println("\n" * "=" ^ 40)
println("SUMMARY")
println("=" ^ 40)
println("2x2 baseline: $(baseline_ok ? "✅ PASS" : "❌ FAIL")")
println("2x3 corrected: $(corrected_ok ? "✅ PASS" : "❌ FAIL")")

if baseline_ok && corrected_ok
    println("\n🎉 SUCCESS: Corrected algorithm appears to fix the bias!")
elseif baseline_ok && !corrected_ok
    println("\n⚠️ PARTIAL: Algorithm still needs work on m > n cases")
else
    println("\n❌ MAJOR ISSUES: Fundamental problems with implementation")
end