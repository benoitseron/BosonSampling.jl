"""
Ultra Simple Test - Just check if algorithm runs without errors
"""

using BosonSampling
using Random
include("corrected_clifford_algorithm.jl")

function ultra_simple_test()
    println("ULTRA SIMPLE TEST")
    println("=" ^ 30)
    
    Random.seed!(42)
    
    # Simplest test: 2 photons, 2 modes
    n, m = 2, 2
    input_state = Input{Bosonic}(first_modes(n, m))
    interf = RandHaar(m)
    
    println("Testing corrected algorithm basic functionality...")
    
    try
        # Just try to generate some samples
        samples = []
        for i in 1:10
            sample = corrected_clifford_sampler(input_state, interf)
            push!(samples, sample)
            print(".")
        end
        println()
        
        println("✅ SUCCESS: Algorithm runs without errors")
        println("Sample outputs:")
        for (i, sample) in enumerate(samples[1:5])
            println("  Sample $i: $sample")
        end
        
        # Check if samples look reasonable
        all_zero = all(s -> all(s .== 0), samples)
        all_same = all(s -> s == samples[1], samples)
        
        if all_zero
            println("❌ WARNING: All samples are zero - likely bug")
        elseif all_same
            println("❌ WARNING: All samples identical - no randomness")
        else
            println("✅ Samples look reasonable (varied, non-zero)")
        end
        
        return true
        
    catch e
        println("❌ FAILURE: Algorithm crashed with error:")
        println("  $e")
        return false
    end
end

# Run test
success = ultra_simple_test()
println("\nResult: $(success ? "PASS" : "FAIL")")