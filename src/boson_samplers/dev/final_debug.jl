"""
Final Debug: Find the exact source of bias

We know:
- 2x2 system: TVD = 0.0012 (perfect)  
- 2x3 system: TVD = 0.14 (biased)

The issue appears when m > n. Let's find exactly where.
"""

using BosonSampling
using Random

include("clifford_paper_implementation.jl")

function test_progressive_systems()
    println("PROGRESSIVE SYSTEM SIZE TEST")
    println("=" ^ 40)
    
    # Test increasing system sizes to find where bias appears
    test_configs = [
        (2, 2, "2 photons, 2 modes"),
        (2, 3, "2 photons, 3 modes"), 
        (2, 4, "2 photons, 4 modes"),
        (3, 3, "3 photons, 3 modes"),
        (3, 4, "3 photons, 4 modes")
    ]
    
    Random.seed!(42)
    n_samples = 5000
    
    for (n, m, description) in test_configs
        println("\nTesting: $description")
        
        input_state = Input{Bosonic}(first_modes(n, m))
        interf = RandHaar(m)
        
        # Generate samples
        counts = Dict{Vector{Int}, Int}()
        for _ in 1:n_samples
            sample = clifford_sampler_paper(input_state, interf)
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
        
        status = tvd < 0.05 ? "✓ GOOD" : "⚠ BIASED"
        println("  TVD: $(round(tvd, digits=4)) $status")
        
        # Show top configurations if biased
        if tvd > 0.05
            println("  Top configurations (empirical vs exact):")
            sorted_counts = sort(collect(counts), by=x->x[2], rev=true)[1:3]
            for (config, count) in sorted_counts
                emp_prob = count / n_samples
                mode_occ = ModeOccupation(config)
                output = FockDetection(mode_occ)
                ev = Event(input_state, output, interf)
                BosonSampling.compute_probability!(ev)
                exact_prob = real(ev.proba_params.probability)
                ratio = emp_prob / exact_prob
                println("    $config: $(round(emp_prob, digits=3)) vs $(round(exact_prob, digits=3)) (ratio $(round(ratio, digits=2)))")
            end
        end
    end
end

function hypothesis_test()
    println("\nHYPOTHESIS TEST")
    println("=" ^ 40)
    
    println("Hypothesis: The bias appears when there are unused output modes")
    println("Testing this by comparing systems with same n but different m")
    
    Random.seed!(123)  # Different seed for independence
    n_samples = 10000
    
    # Test 2 photons with different numbers of modes
    for m in [2, 3, 4, 5]
        println("\n2 photons, $m modes:")
        
        input_state = Input{Bosonic}(first_modes(2, m))
        interf = RandHaar(m)
        
        # Quick TVD calculation
        counts = Dict{Vector{Int}, Int}()
        for _ in 1:n_samples
            sample = clifford_sampler_paper(input_state, interf)
            counts[sample] = get(counts, sample, 0) + 1
        end
        
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
        
        println("  TVD: $(round(tvd, digits=4))")
        
        # Count how many modes actually get photons
        used_modes = Set{Int}()
        for (config, count) in counts
            for (mode_idx, photon_count) in enumerate(config)
                if photon_count > 0
                    push!(used_modes, mode_idx)
                end
            end
        end
        println("  Modes that received photons: $(length(used_modes))/$m")
    end
end

function check_built_in_vs_paper()
    println("\nBUILT-IN vs PAPER COMPARISON")
    println("=" ^ 40)
    
    # Direct comparison on problematic 2x3 system
    n, m = 2, 3
    Random.seed!(42)
    
    input_state = Input{Bosonic}(first_modes(n, m))
    interf = RandHaar(m)
    
    n_samples = 5000
    
    # Paper implementation
    paper_counts = Dict{Vector{Int}, Int}()
    for _ in 1:n_samples
        sample = clifford_sampler_paper(input_state, interf)
        paper_counts[sample] = get(paper_counts, sample, 0) + 1
    end
    
    # Built-in implementation  
    builtin_counts = Dict{Vector{Int}, Int}()
    for _ in 1:n_samples
        ev = Event(input_state, FockSample(), interf)
        BosonSampling.sample!(ev)
        sample = ev.output_measurement.s.state
        builtin_counts[sample] = get(builtin_counts, sample, 0) + 1
    end
    
    # Compare the two implementations
    println("Paper vs Built-in sampler comparison:")
    println("Config | Paper % | Built-in % | Ratio")
    println("_" ^ 35)
    
    all_configs = unique([keys(paper_counts)..., keys(builtin_counts)...])
    
    for config in sort(all_configs)
        paper_freq = 100 * get(paper_counts, config, 0) / n_samples
        builtin_freq = 100 * get(builtin_counts, config, 0) / n_samples
        
        if paper_freq > 0.1 || builtin_freq > 0.1
            ratio = paper_freq > 0 && builtin_freq > 0 ? paper_freq / builtin_freq : NaN
            ratio_str = isnan(ratio) ? "---" : "$(round(ratio, digits=2))"
            println("$config | $(round(paper_freq, digits=1)) | $(round(builtin_freq, digits=1)) | $ratio_str")
        end
    end
    
    # If they're similar, both have the same bias
    # If they're different, our implementation has a unique bug
end

# Run all tests
test_progressive_systems()
hypothesis_test()
check_built_in_vs_paper()

println("\n" * "=" ^ 40)
println("FINAL DEBUG COMPLETE")
println("=" ^ 40)