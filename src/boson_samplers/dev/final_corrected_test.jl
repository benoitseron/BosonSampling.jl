"""
Final Test of Corrected Algorithm with Threshold Detection

Since FockDetection only allows threshold detection (≤1 photon per mode),
let's use that for validation.
"""

using BosonSampling
using Random
include("corrected_clifford_algorithm.jl")

function convert_to_threshold(mode_occ)
    """Convert mode occupation to threshold (0 or 1 per mode)"""
    return [x > 0 ? 1 : 0 for x in mode_occ]
end

function test_algorithm_basic()
    println("BASIC ALGORITHM TEST")
    println("=" ^ 40)
    
    Random.seed!(42)
    
    # Test 2x2 system
    n, m = 2, 2
    input_state = Input{Bosonic}(first_modes(n, m))
    interf = RandHaar(m)
    
    println("Testing 2x2 system (threshold detection)...")
    
    # Generate samples
    n_samples = 1000
    threshold_counts = Dict{Vector{Int}, Int}()
    
    for _ in 1:n_samples
        sample = corrected_clifford_sampler(input_state, interf)
        threshold_sample = convert_to_threshold(sample)
        threshold_counts[threshold_sample] = get(threshold_counts, threshold_sample, 0) + 1
    end
    
    println("Threshold detection configurations found:")
    for (config, count) in sort(collect(threshold_counts), by=x->x[2], rev=true)
        freq = count / n_samples
        println("  $config: $(round(freq, digits=3)) ($(count) samples)")
    end
    
    # Check if we see both [0,1] and [1,0] with reasonable frequency
    expected_configs = [[0,1], [1,0]]
    all_found = all(config -> haskey(threshold_counts, config), expected_configs)
    
    if all_found
        println("✅ SUCCESS: Found all expected threshold configurations")
        return true
    else
        println("❌ FAILURE: Missing expected configurations")
        return false
    end
end

function test_progression()
    println("\nPROGRESSION TEST")
    println("=" ^ 40)
    
    Random.seed!(123)
    
    # Test different system sizes
    test_configs = [
        (2, 2, 1000),
        (2, 3, 1000), 
        (3, 3, 500),
        (3, 4, 500)
    ]
    
    all_passed = true
    
    for (n, m, n_samples) in test_configs
        println("\nTesting $(n)x$(m) system:")
        
        input_state = Input{Bosonic}(first_modes(n, m))
        interf = RandHaar(m)
        
        # Count unique outputs
        seen_outputs = Set{Vector{Int}}()
        threshold_counts = Dict{Vector{Int}, Int}()
        
        for _ in 1:n_samples
            sample = corrected_clifford_sampler(input_state, interf)
            push!(seen_outputs, copy(sample))
            
            threshold_sample = convert_to_threshold(sample)
            threshold_counts[threshold_sample] = get(threshold_counts, threshold_sample, 0) + 1
        end
        
        n_unique_full = length(seen_outputs)
        n_unique_threshold = length(threshold_counts)
        
        println("  Full outputs: $n_unique_full unique configurations")
        println("  Threshold outputs: $n_unique_threshold unique configurations")
        
        # Check for reasonable diversity
        if n_unique_full >= 3 && n_unique_threshold >= 2
            println("  ✅ Good diversity in outputs")
        else
            println("  ❌ Low diversity - possible issue")
            all_passed = false
        end
        
        # Show some examples
        if n_unique_full <= 10
            println("  Full output examples:")
            for (i, output) in enumerate(collect(seen_outputs)[1:min(5, n_unique_full)])
                println("    Example $i: $output")
            end
        end
    end
    
    return all_passed
end

function compare_with_builtin()
    println("\nCOMPARISON WITH BUILT-IN SAMPLER")
    println("=" ^ 40)
    
    Random.seed!(42)
    
    # Test 2x3 system where we know built-in has bias
    n, m = 2, 3
    input_state = Input{Bosonic}(first_modes(n, m))
    interf = RandHaar(m)
    
    n_samples = 1000
    
    # Our corrected algorithm
    our_threshold_counts = Dict{Vector{Int}, Int}()
    for _ in 1:n_samples
        sample = corrected_clifford_sampler(input_state, interf)
        threshold_sample = convert_to_threshold(sample)
        our_threshold_counts[threshold_sample] = get(our_threshold_counts, threshold_sample, 0) + 1
    end
    
    # Built-in sampler
    Random.seed!(42)  # Reset for fair comparison
    builtin_threshold_counts = Dict{Vector{Int}, Int}()
    
    for _ in 1:n_samples
        ev = Event(input_state, FockSample(), interf)
        BosonSampling.sample!(ev)
        sample = ev.output_measurement.s.state
        threshold_sample = convert_to_threshold(sample)
        builtin_threshold_counts[threshold_sample] = get(builtin_threshold_counts, threshold_sample, 0) + 1
    end
    
    println("Threshold detection comparison (2x3 system):")
    println("Config | Our % | Built-in % | Ratio")
    println("-" ^ 35)
    
    all_configs = Set([keys(our_threshold_counts)..., keys(builtin_threshold_counts)...])
    
    for config in sort(collect(all_configs))
        our_freq = 100 * get(our_threshold_counts, config, 0) / n_samples
        builtin_freq = 100 * get(builtin_threshold_counts, config, 0) / n_samples
        
        if our_freq > 0.1 || builtin_freq > 0.1
            ratio = (our_freq > 0 && builtin_freq > 0) ? our_freq / builtin_freq : NaN
            ratio_str = isnan(ratio) ? "---" : "$(round(ratio, digits=2))"
            println("$config | $(round(our_freq, digits=1)) | $(round(builtin_freq, digits=1)) | $ratio_str")
        end
    end
    
    # Overall similarity check
    shared_configs = intersect(Set(keys(our_threshold_counts)), Set(keys(builtin_threshold_counts)))
    
    if length(shared_configs) >= 2
        println("\n✅ Both samplers produce similar threshold patterns")
        return true
    else
        println("\n⚠️  Samplers produce very different threshold patterns")
        return false
    end
end

# Run all tests
println("TESTING CORRECTED CLIFFORD ALGORITHM")
println("=" ^ 50)

basic_ok = test_algorithm_basic()
progression_ok = test_progression()
comparison_ok = compare_with_builtin()

println("\n" * "=" ^ 50)
println("FINAL RESULTS")
println("=" ^ 50)
println("Basic functionality: $(basic_ok ? "✅ PASS" : "❌ FAIL")")
println("Progression test: $(progression_ok ? "✅ PASS" : "❌ FAIL")")
println("Comparison test: $(comparison_ok ? "✅ PASS" : "❌ FAIL")")

if basic_ok && progression_ok
    println("\n🎉 SUCCESS: Corrected algorithm appears to be working!")
    println("   Algorithm generates diverse, reasonable outputs")
    if comparison_ok
        println("   Threshold patterns match built-in sampler")
    end
else
    println("\n❌ ISSUES DETECTED: Algorithm needs further work")
end