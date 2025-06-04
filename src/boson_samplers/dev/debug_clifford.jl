"""
Debug Clifford Algorithm Implementation

Systematically debug each step of the algorithm to find the source of the bias.
"""

using BosonSampling
using Random
using LinearAlgebra

include("clifford_paper_implementation.jl")

Random.seed!(42)

"""
    debug_single_run()

Debug a single run of the algorithm step by step.
"""
function debug_single_run()
    println("="^60)
    println("DEBUGGING SINGLE ALGORITHM RUN")
    println("="^60)
    
    # Simple 2x2 system for easier debugging
    n, m = 2, 2
    input_state = Input{Bosonic}(first_modes(n, m))
    interf = RandHaar(m)
    
    U = interf.U
    input_arrangement = BosonSampling.fill_arrangement(input_state)
    A = U[:, input_arrangement]
    
    println("Unitary matrix U:")
    display(U)
    println("\nInput arrangement: $input_arrangement")
    println("Submatrix A:")
    display(A)
    
    # Step-by-step algorithm execution
    println("\n" * "="^40)
    println("ALGORITHM EXECUTION")
    println("="^40)
    
    # Step 1: r ← ∅
    r = Int[]
    println("Step 1: r = $r")
    
    # Step 2: A ← PERMUTE(A)
    perm = randperm(n)
    A_orig = copy(A)
    A = permute_columns(A, n)
    println("\nStep 2: Column permutation = $perm")
    println("Original A:")
    display(A_orig)
    println("Permuted A:")
    display(A)
    
    # Step 3: wi ← |ai,1|², i ∈ [m]
    w = abs2.(A[:, 1])
    println("\nStep 3: Weights w = $w")
    println("Sum of weights: $(sum(w))")
    
    # Step 4: x ← SAMPLE(w)
    x = wsample(1:m, Weights(w))
    println("\nStep 4: Sampled x = $x")
    
    # Step 5: r ← (r, x)
    push!(r, x)
    println("\nStep 5: r = $r")
    
    # Step 6-12: FOR k ← 2 TO n
    for k in 2:n
        println("\n" * "-"^30)
        println("Loop iteration k = $k")
        println("-"^30)
        
        # Step 7: B_k^x ← A^[k]
        B_k_x = A[r, 1:k]
        println("Step 7: B_k_x (rows $r, cols 1:$k):")
        display(B_k_x)
        
        # Step 8: COMPUTE {Per B_{k,ℓ}^x, ℓ ∈ [k]}
        permanents = compute_permanents_lemma2(B_k_x, k)
        println("\nStep 8: Permanents = $permanents")
        
        # Debug: Show each permanent calculation
        for ℓ in 1:k
            cols_to_keep = [1:(ℓ-1); (ℓ+1):k]
            println("  Permanent $ℓ (remove col $ℓ, keep cols $cols_to_keep):")
            if length(cols_to_keep) > 0
                submatrix = B_k_x[:, cols_to_keep]
                println("    Submatrix:")
                display(submatrix)
                perm_val = BosonSampling.permanent(submatrix)
                println("    Permanent = $perm_val")
            else
                println("    Empty matrix, permanent = 1")
            end
        end
        
        # Step 9: wi ← |∑_{ℓ=1}^k ai,ℓ Per B_{k,ℓ}^x|², i ∈ [m]
        w = zeros(Float64, m)
        println("\nStep 9: Computing weights...")
        for i in 1:m
            amplitude = zero(ComplexF64)
            for ℓ in 1:k
                term = A[i, ℓ] * permanents[ℓ]
                amplitude += term
                println("  Mode $i, term $ℓ: A[$i,$ℓ] * perm[$ℓ] = $(A[i, ℓ]) * $(permanents[ℓ]) = $term")
            end
            w[i] = abs2(amplitude)
            println("  Mode $i: amplitude = $amplitude, weight = $(w[i])")
        end
        
        println("Final weights: $w")
        println("Sum of weights: $(sum(w))")
        
        # Step 10: x ← SAMPLE(w)
        x = wsample(1:m, Weights(w))
        println("\nStep 10: Sampled x = $x")
        
        # Step 11: r ← (r, x)
        push!(r, x)
        println("\nStep 11: r = $r")
    end
    
    # Step 13: z ← INCSORT(r)
    z = sort(r)
    println("\n" * "="^30)
    println("Step 13: Final result z = $z")
    
    # Convert to occupancy vector
    occupancy = zeros(Int, m)
    for mode in z
        occupancy[mode] += 1
    end
    println("Occupancy vector: $occupancy")
    
    return z, occupancy
end

"""
    compare_with_exact_probability()

Compare our sampled probability with exact calculation.
"""
function compare_with_exact_probability()
    println("\n" * "="^60)
    println("COMPARING WITH EXACT PROBABILITY")
    println("="^60)
    
    # Use same system
    n, m = 2, 2
    input_state = Input{Bosonic}(first_modes(n, m))
    interf = RandHaar(m)
    
    # Generate one sample and compute its exact probability
    sample = clifford_sampler_paper(input_state, interf)
    println("Generated sample: $sample")
    
    # Create event and compute exact probability
    mode_occ = ModeOccupation(sample)
    output = FockDetection(mode_occ)
    ev = Event(input_state, output, interf)
    BosonSampling.compute_probability!(ev)
    exact_prob = real(ev.proba_params.probability)
    
    println("Exact probability: $exact_prob")
    
    # Generate many samples and see empirical frequency
    n_samples = 10000
    counts = Dict{Vector{Int}, Int}()
    
    for _ in 1:n_samples
        s = clifford_sampler_paper(input_state, interf)
        counts[s] = get(counts, s, 0) + 1
    end
    
    empirical_prob = get(counts, sample, 0) / n_samples
    println("Empirical probability: $empirical_prob")
    println("Ratio (empirical/exact): $(empirical_prob/exact_prob)")
    
    # Show all probabilities
    println("\nAll configurations:")
    println("Config | Exact | Empirical | Ratio")
    println("-" * 35)
    
    for (config, count) in sort(collect(counts), by=x->x[2], rev=true)
        emp_prob = count / n_samples
        
        # Compute exact probability for this config
        mode_occ = ModeOccupation(config)
        output = FockDetection(mode_occ)
        ev = Event(input_state, output, interf)
        BosonSampling.compute_probability!(ev)
        ex_prob = real(ev.proba_params.probability)
        
        ratio = emp_prob / ex_prob
        println("$config | $(round(ex_prob, digits=4)) | $(round(emp_prob, digits=4)) | $(round(ratio, digits=2))")
    end
end

"""
    check_probability_normalization()

Check if the probabilities in each step are properly normalized.
"""
function check_probability_normalization()
    println("\n" * "="^60)
    println("CHECKING PROBABILITY NORMALIZATION")
    println("="^60)
    
    # Test system
    n, m = 2, 3
    input_state = Input{Bosonic}(first_modes(n, m))
    interf = RandHaar(m)
    
    U = interf.U
    input_arrangement = BosonSampling.fill_arrangement(input_state)
    A = U[:, input_arrangement]
    
    # Check if columns of U are properly normalized
    println("Checking unitary matrix properties:")
    println("U† U =")
    display(U' * U)
    
    println("\nColumn norms of U:")
    for i in 1:size(U, 2)
        norm_val = norm(U[:, i])
        println("  Column $i: norm = $norm_val")
    end
    
    # Check our submatrix A
    println("\nSubmatrix A properties:")
    println("A:")
    display(A)
    
    println("\nColumn norms of A:")
    for i in 1:size(A, 2)
        norm_val = norm(A[:, i])
        println("  Column $i: norm = $norm_val")
    end
    
    # Test first step normalization
    w1 = abs2.(A[:, 1])
    println("\nFirst step weights: $w1")
    println("Sum: $(sum(w1))")
    
    # Test if this should equal 1
    expected_sum = sum(abs2.(U[:, input_arrangement[1]]))
    println("Expected sum (from unitary property): $expected_sum")
end

"""
    test_permanent_calculations()

Verify our permanent calculations are correct.
"""
function test_permanent_calculations()
    println("\n" * "="^60)
    println("TESTING PERMANENT CALCULATIONS")
    println("="^60)
    
    # Test small matrices
    test_matrices = [
        [1.0+0.0im],
        [1.0+0.0im 0.5+0.0im; 0.5+0.0im 1.0+0.0im],
        [1.0+0.0im 0.5+0.0im 0.2+0.0im; 0.3+0.0im 1.0+0.0im 0.4+0.0im; 0.1+0.0im 0.6+0.0im 1.0+0.0im]
    ]
    
    for (i, M) in enumerate(test_matrices)
        println("\nTest matrix $i:")
        display(M)
        
        # Our calculation
        our_perm = permanent_fast(M)
        
        # Built-in calculation
        builtin_perm = BosonSampling.permanent(M)
        
        println("Our permanent: $our_perm")
        println("Built-in permanent: $builtin_perm")
        println("Match: $(isapprox(our_perm, builtin_perm))")
        
        if !isapprox(our_perm, builtin_perm, rtol=1e-10)
            println("❌ MISMATCH DETECTED!")
        end
    end
end

# Run all debugging functions
println("STARTING COMPREHENSIVE DEBUG")
println("="^60)

# debug_single_run()
compare_with_exact_probability()
check_probability_normalization()
test_permanent_calculations()