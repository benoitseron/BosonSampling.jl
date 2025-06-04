using BosonSampling
using Test
using Statistics
using Random
using ProgressMeter

# Test 1: Verify Clifford sampler produces valid output configurations
function test_clifford_output_validity()
    println("\n=== Test: Clifford Output Validity ===")
    
    n = 4  # photons
    m = 6  # modes
    
    interf = RandHaar(m)
    input_state = Input{Bosonic}(first_modes(n, m))
    
    n_samples = 1000
    all_valid = true
    
    for i in 1:n_samples
        sample = clifford_sampler_unoptimised(input_state, interf, occupancy_vector=true)
        
        # Check: correct number of photons
        if sum(sample) != n
            println("ERROR: Sample $i has $(sum(sample)) photons, expected $n")
            all_valid = false
        end
        
        # Check: valid mode occupation (non-negative integers)
        if any(s -> s < 0 || !isinteger(s), sample)
            println("ERROR: Sample $i has invalid mode occupation: $sample")
            all_valid = false
        end
        
        # Check: correct length
        if length(sample) != m
            println("ERROR: Sample $i has length $(length(sample)), expected $m")
            all_valid = false
        end
    end
    
    @test all_valid
    println("All $n_samples samples are valid: $all_valid")
    
    return all_valid
end

# Test 2: Statistical test - check if Clifford produces bunching
function test_clifford_bunching()
    println("\n=== Test: Clifford Bunching Behavior ===")
    
    n = 6
    m = 10
    n_samples = 10000
    
    interf = RandHaar(m)
    input_state = Input{Bosonic}(first_modes(n, m))
    
    # Collect samples
    samples = []
    for i in 1:n_samples
        sample = clifford_sampler_unoptimised(input_state, interf, occupancy_vector=true)
        push!(samples, sample)
    end
    
    # Calculate bunching metric: average number of modes with >1 photon
    bunching_scores = [sum(s .> 1) for s in samples]
    avg_bunching = mean(bunching_scores)
    
    # For bosons, we expect more bunching than classical particles
    # With n=6, m=10, classical would give very few collisions
    classical_expected = n * (n-1) / (2*m)  # Rough approximation
    
    println("Average modes with >1 photon: $avg_bunching")
    println("Classical expectation (approx): $classical_expected")
    println("Bunching enhancement: $(avg_bunching / classical_expected)x")
    
    @test avg_bunching > classical_expected
    
    return avg_bunching, classical_expected
end

# Test 3: Compare Clifford with exact distribution for small system
function test_clifford_exact_comparison()
    println("\n=== Test: Clifford vs Exact Distribution ===")
    
    # Small system for exact comparison
    n = 3
    m = 3
    n_samples = 50000
    
    interf = Fourier(m)  # Use Fourier for more predictable structure
    input_state = Input{Bosonic}(first_modes(n, m))
    
    # Generate Clifford samples
    println("Generating Clifford samples...")
    clifford_counts = Dict{Vector{Int}, Int}()
    
    @showprogress for i in 1:n_samples
        sample = clifford_sampler_unoptimised(input_state, interf, occupancy_vector=true)
        key = Vector{Int}(sample)
        clifford_counts[key] = get(clifford_counts, key, 0) + 1
    end
    
    # Compute exact probabilities
    println("Computing exact probabilities...")
    exact_probs = Dict{Vector{Int}, Float64}()
    
    all_configs = all_mode_configurations(n, m, only_photon_number_conserving=true)
    
    for config in all_configs
        o = FockDetection(ModeOccupation(config))
        ev = Event(input_state, o, interf)
        compute_probability!(ev)
        exact_probs[config] = ev.proba_params.probability
    end
    
    # Compare distributions
    println("\nTop 5 configurations:")
    println("Config\t\tExact Prob\tClifford Freq\tRatio")
    
    sorted_configs = sort(collect(exact_probs), by=x->x[2], rev=true)
    
    total_tvd = 0.0
    for (i, (config, exact_p)) in enumerate(sorted_configs[1:min(5, length(sorted_configs))])
        clifford_freq = get(clifford_counts, config, 0) / n_samples
        ratio = clifford_freq / exact_p
        println("$config\t$(round(exact_p, digits=4))\t$(round(clifford_freq, digits=4))\t$(round(ratio, digits=2))")
        
        total_tvd += abs(exact_p - clifford_freq)
    end
    
    # Add remaining configs to TVD
    for (config, exact_p) in exact_probs
        if !(config in [c[1] for c in sorted_configs[1:min(5, length(sorted_configs))]])
            clifford_freq = get(clifford_counts, config, 0) / n_samples
            total_tvd += abs(exact_p - clifford_freq)
        end
    end
    
    println("\nTotal Variation Distance: $total_tvd")
    
    @test total_tvd < 0.05  # Should be close to exact distribution
    
    return total_tvd, clifford_counts, exact_probs
end

# Test 4: Verify Clifford respects interferometer symmetries
function test_clifford_symmetries()
    println("\n=== Test: Clifford Symmetry Preservation ===")
    
    n = 4
    m = 4
    n_samples = 10000
    
    # Use Fourier transform which has known symmetries
    interf = Fourier(m)
    input_state = Input{Bosonic}(first_modes(n, m))
    
    # Collect output mode statistics
    mode_counts = zeros(m)
    
    for i in 1:n_samples
        sample = clifford_sampler_unoptimised(input_state, interf, occupancy_vector=true)
        mode_counts .+= sample
    end
    
    mode_probs = mode_counts / (n * n_samples)
    
    println("Output mode probabilities:")
    for (i, p) in enumerate(mode_probs)
        println("Mode $i: $(round(p, digits=4))")
    end
    
    # For Fourier transform with uniform input, output should be uniform
    expected_prob = 1/m
    max_deviation = maximum(abs.(mode_probs .- expected_prob))
    
    println("\nExpected probability per mode: $expected_prob")
    println("Maximum deviation: $max_deviation")
    
    @test max_deviation < 0.02  # Should be approximately uniform
    
    return mode_probs, max_deviation
end

# Test 5: Performance comparison
function test_clifford_performance()
    println("\n=== Test: Clifford Performance ===")
    
    sizes = [(4,4), (6,8), (8,12), (10,20)]
    n_samples = 1000
    
    println("System Size\tTime per sample (ms)")
    
    for (n, m) in sizes
        interf = RandHaar(m)
        input_state = Input{Bosonic}(first_modes(n, m))
        
        # Warmup
        clifford_sampler_unoptimised(input_state, interf)
        
        # Time the sampling
        t_start = time()
        for i in 1:n_samples
            clifford_sampler_unoptimised(input_state, interf)
        end
        t_elapsed = time() - t_start
        
        time_per_sample = t_elapsed / n_samples * 1000  # Convert to ms
        println("n=$n, m=$m\t$(round(time_per_sample, digits=3))")
    end
end

# Run all tests
function run_all_tests()
    Random.seed!(42)
    
    println("Running Clifford Sampler Test Suite")
    println("=" ^ 40)
    
    @testset "Clifford Sampler Tests" begin
        @test test_clifford_output_validity()
        
        avg_bunch, classical = test_clifford_bunching()
        @test avg_bunch > classical
        
        tvd, _, _ = test_clifford_exact_comparison()
        @test tvd < 0.05
        
        _, max_dev = test_clifford_symmetries()
        @test max_dev < 0.02
    end
    
    test_clifford_performance()
    
    println("\n" * "=" ^ 40)
    println("All tests completed!")
end

# Example usage:
# run_all_tests()

# Or run individual tests:
# test_clifford_output_validity()
# test_clifford_bunching()
# test_clifford_exact_comparison()
# test_clifford_symmetries()
# test_clifford_performance()