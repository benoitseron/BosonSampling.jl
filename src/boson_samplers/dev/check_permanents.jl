"""
Check Permanent Calculations

The 2x2 system works perfectly (TVD = 0.0012) but 2x3 has bias (TVD = 0.14).
This suggests the issue is in permanent calculations for larger systems.
"""

using BosonSampling
using Random
using LinearAlgebra

include("clifford_paper_implementation.jl")

Random.seed!(42)

function test_permanent_accuracy()
    println("TESTING PERMANENT CALCULATION ACCURACY")
    println("=" ^ 50)
    
    # Test with matrices of different sizes
    test_sizes = [1, 2, 3, 4]
    
    for n in test_sizes
        println("\nTesting $(n)x$(n) matrices:")
        
        # Generate random complex matrices
        for trial in 1:5
            M = randn(ComplexF64, n, n)
            
            # Our calculation
            our_perm = permanent_fast(M)
            
            # Built-in calculation  
            builtin_perm = BosonSampling.permanent(M)
            
            # Relative error
            if abs(builtin_perm) > 1e-12
                rel_error = abs(our_perm - builtin_perm) / abs(builtin_perm)
            else
                rel_error = abs(our_perm - builtin_perm)
            end
            
            println("  Trial $trial: rel_error = $(round(rel_error, sigdigits=3))")
            
            if rel_error > 1e-10
                println("    ❌ Large error detected!")
                println("    Our result: $our_perm")
                println("    Built-in: $builtin_perm")
            end
        end
    end
end

function test_specific_problematic_cases()
    println("\nTESTING SPECIFIC PROBLEMATIC CASES")
    println("=" ^ 50)
    
    # Test the exact scenario from 2x3 system
    n, m = 2, 3
    input_state = Input{Bosonic}(first_modes(n, m))
    interf = RandHaar(m)
    
    U = interf.U
    input_arrangement = BosonSampling.fill_arrangement(input_state)
    A = U[:, input_arrangement]
    
    println("Testing with 2x3 system matrix:")
    println("A = ")
    display(A)
    
    # Test each step of the algorithm for potential issues
    println("\nStep-by-step permanent calculations:")
    
    # Permute columns
    perm = randperm(n)
    A_perm = permute_columns(A, n)
    println("After permutation $perm:")
    display(A_perm)
    
    # Simulate first sampling step
    w1 = abs2.(A_perm[:, 1])
    println("\nFirst step weights: $w1")
    println("Sum: $(sum(w1)) (should be ≈ 1.0)")
    
    # Sample first mode (let's say mode 1 for testing)
    r = [1]
    
    # Second step permanent calculations
    println("\nSecond step (k=2, r=$r):")
    B_k_x = A_perm[r, 1:2]
    println("B_k_x = ")
    display(B_k_x)
    
    # Compute permanents for second step
    permanents = compute_permanents_lemma2(B_k_x, 2)
    println("Permanents: $permanents")
    
    # Check each permanent manually
    for ℓ in 1:2
        cols_to_keep = [1:(ℓ-1); (ℓ+1):2]
        println("\nPermanent $ℓ (remove col $ℓ, keep cols $cols_to_keep):")
        
        if length(cols_to_keep) > 0
            submatrix = B_k_x[:, cols_to_keep]
            println("  Submatrix:")
            display(submatrix)
            
            our_perm = permanent_fast(submatrix)
            builtin_perm = BosonSampling.permanent(submatrix)
            
            println("  Our permanent: $our_perm")
            println("  Built-in permanent: $builtin_perm")
            println("  Match: $(isapprox(our_perm, builtin_perm))")
        else
            println("  Empty matrix, permanent = 1")
        end
    end
    
    # Check weight calculation
    println("\nWeight calculations:")
    w = zeros(Float64, m)
    for i in 1:m
        amplitude = zero(ComplexF64)
        for ℓ in 1:2
            term = A_perm[i, ℓ] * permanents[ℓ]
            amplitude += term
            println("  Mode $i, term $ℓ: $(A_perm[i, ℓ]) * $(permanents[ℓ]) = $term")
        end
        w[i] = abs2(amplitude)
        println("  Mode $i: total amplitude = $amplitude, weight = $(w[i])")
    end
    
    println("Final weights: $w")
    println("Sum: $(sum(w))")
end

function compare_2x2_vs_2x3()
    println("\nCOMPARING 2x2 vs 2x3 SYSTEMS")
    println("=" ^ 50)
    
    # Test both systems with same random seed
    Random.seed!(42)
    
    println("2x2 system:")
    n, m = 2, 2
    input_state = Input{Bosonic}(first_modes(n, m))
    interf = RandHaar(m)
    
    n_samples = 5000
    counts_2x2 = Dict{Vector{Int}, Int}()
    
    for _ in 1:n_samples
        sample = clifford_sampler_paper(input_state, interf)
        counts_2x2[sample] = get(counts_2x2, sample, 0) + 1
    end
    
    # Compute TVD for 2x2
    tvd_2x2 = 0.0
    for (config, count) in counts_2x2
        emp_prob = count / n_samples
        mode_occ = ModeOccupation(config)
        output = FockDetection(mode_occ)
        ev = Event(input_state, output, interf)
        BosonSampling.compute_probability!(ev)
        exact_prob = real(ev.proba_params.probability)
        tvd_2x2 += abs(emp_prob - exact_prob)
    end
    tvd_2x2 /= 2
    
    println("TVD for 2x2: $(round(tvd_2x2, digits=4))")
    
    # Reset and test 2x3
    Random.seed!(42)
    
    println("\n2x3 system:")
    n, m = 2, 3
    input_state = Input{Bosonic}(first_modes(n, m))
    interf = RandHaar(m)
    
    counts_2x3 = Dict{Vector{Int}, Int}()
    
    for _ in 1:n_samples
        sample = clifford_sampler_paper(input_state, interf)
        counts_2x3[sample] = get(counts_2x3, sample, 0) + 1
    end
    
    # Compute TVD for 2x3
    tvd_2x3 = 0.0
    for (config, count) in counts_2x3
        emp_prob = count / n_samples
        mode_occ = ModeOccupation(config)
        output = FockDetection(mode_occ)
        ev = Event(input_state, output, interf)
        BosonSampling.compute_probability!(ev)
        exact_prob = real(ev.proba_params.probability)
        tvd_2x3 += abs(emp_prob - exact_prob)
    end
    tvd_2x3 /= 2
    
    println("TVD for 2x3: $(round(tvd_2x3, digits=4))")
    
    println("\nConclusion:")
    if tvd_2x2 < 0.01 && tvd_2x3 > 0.1
        println("✓ Confirms issue appears when m > n")
        println("  Likely problem in permanent calculations for larger systems")
    elseif tvd_2x2 < 0.01 && tvd_2x3 < 0.01
        println("✓ Both systems work well - issue may be random seed dependent")
    else
        println("⚠ Both systems have issues - fundamental algorithm problem")
    end
end

# Run tests
test_permanent_accuracy()
test_specific_problematic_cases()
compare_2x2_vs_2x3()