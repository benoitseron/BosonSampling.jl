"""
Simple Debug of Clifford Algorithm

Focus on finding the source of the systematic bias.
"""

using BosonSampling
using Random
using LinearAlgebra

include("clifford_paper_implementation.jl")

Random.seed!(42)

function test_2x2_system()
    println("DEBUGGING 2x2 SYSTEM")
    println("=" ^ 40)
    
    # Simplest possible system: 2 photons, 2 modes
    n, m = 2, 2
    input_state = Input{Bosonic}(first_modes(n, m))
    interf = RandHaar(m)
    
    println("Testing 2 photons in 2 modes...")
    
    # Generate many samples
    n_samples = 10000
    counts = Dict{Vector{Int}, Int}()
    
    for _ in 1:n_samples
        sample = clifford_sampler_paper(input_state, interf)
        counts[sample] = get(counts, sample, 0) + 1
    end
    
    # Compare with exact probabilities
    println("\nEmpirical vs Exact comparison:")
    println("Config | Empirical | Exact | Ratio")
    println("_" ^ 35)
    
    total_tvd = 0.0
    
    for (config, count) in sort(collect(counts), by=x->x[2], rev=true)
        emp_prob = count / n_samples
        
        # Compute exact probability
        mode_occ = ModeOccupation(config)
        output = FockDetection(mode_occ)
        ev = Event(input_state, output, interf)
        BosonSampling.compute_probability!(ev)
        exact_prob = real(ev.proba_params.probability)
        
        ratio = emp_prob / exact_prob
        error = abs(emp_prob - exact_prob)
        total_tvd += error
        
        println("$config | $(round(emp_prob, digits=4)) | $(round(exact_prob, digits=4)) | $(round(ratio, digits=2))")
    end
    
    total_tvd /= 2
    println("\nTVD: $(round(total_tvd, digits=4))")
    
    return counts
end

function test_manual_calculation()
    println("\nMANUAL CALCULATION TEST")
    println("=" ^ 40)
    
    # Create a simple known unitary
    θ = π/4  # 45 degree rotation
    U = [cos(θ) -sin(θ); sin(θ) cos(θ)]
    
    println("Using simple rotation matrix:")
    println("U = ")
    for i in 1:2
        println("  [$(round(U[i,1], digits=3)) $(round(U[i,2], digits=3))]")
    end
    
    # Input: [1,1] (one photon in each of first 2 modes)
    # For 2 photons in 2 modes, there are 3 possible outputs:
    # [2,0], [1,1], [0,2]
    
    # Manual calculation of exact probabilities
    println("\nManual exact probability calculation:")
    
    # For [2,0]: permanent of [[U[1,1], U[1,2]], [U[1,1], U[1,2]]]
    submatrix_20 = [U[1,1] U[1,2]; U[1,1] U[1,2]]
    prob_20 = abs2(BosonSampling.permanent(submatrix_20))
    
    # For [1,1]: permanent of [[U[1,1], U[1,2]], [U[2,1], U[2,2]]]
    submatrix_11 = [U[1,1] U[1,2]; U[2,1] U[2,2]]
    prob_11 = abs2(BosonSampling.permanent(submatrix_11))
    
    # For [0,2]: permanent of [[U[2,1], U[2,2]], [U[2,1], U[2,2]]]
    submatrix_02 = [U[2,1] U[2,2]; U[2,1] U[2,2]]
    prob_02 = abs2(BosonSampling.permanent(submatrix_02))
    
    println("[2,0]: permanent = $(BosonSampling.permanent(submatrix_20)), prob = $(round(prob_20, digits=4))")
    println("[1,1]: permanent = $(BosonSampling.permanent(submatrix_11)), prob = $(round(prob_11, digits=4))")
    println("[0,2]: permanent = $(BosonSampling.permanent(submatrix_02)), prob = $(round(prob_02, digits=4))")
    println("Total: $(round(prob_20 + prob_11 + prob_02, digits=4))")
    
    # Now test our algorithm with this same unitary
    println("\nTesting our algorithm with this unitary:")
    
    # Create interferometer with our matrix
    interf = Interferometer(U)
    input_state = Input{Bosonic}(first_modes(2, 2))
    
    # Generate samples
    n_samples = 10000
    counts = Dict{Vector{Int}, Int}()
    
    for _ in 1:n_samples
        sample = clifford_sampler_paper(input_state, interf)
        counts[sample] = get(counts, sample, 0) + 1
    end
    
    println("\nComparison:")
    println("Config | Manual Exact | Algorithm Empirical | Ratio")
    println("_" ^ 50)
    
    exact_probs = Dict([2,0] => prob_20, [1,1] => prob_11, [0,2] => prob_02)
    
    for config in [[2,0], [1,1], [0,2]]
        exact = exact_probs[config]
        empirical = get(counts, config, 0) / n_samples
        ratio = empirical / exact
        
        println("$config | $(round(exact, digits=4)) | $(round(empirical, digits=4)) | $(round(ratio, digits=2))")
    end
end

function check_first_step_probabilities()
    println("\nCHECKING FIRST STEP PROBABILITIES")
    println("=" ^ 40)
    
    # Test if the first step sampling is correct
    n, m = 2, 3
    input_state = Input{Bosonic}(first_modes(n, m))
    interf = RandHaar(m)
    
    U = interf.U
    input_arrangement = BosonSampling.fill_arrangement(input_state)
    A = U[:, input_arrangement]
    
    println("Unitary submatrix A:")
    for i in 1:m
        println("  Row $i: [$(round(A[i,1], digits=3)) $(round(A[i,2], digits=3))]")
    end
    
    # First step: sample according to |A[i,1]|²
    w1 = abs2.(A[:, 1])
    println("\nFirst photon weights: $w1")
    println("Sum: $(sum(w1))")
    
    # This should equal 1 if A comes from a unitary
    expected_sum = sum(abs2.(U[:, input_arrangement[1]]))
    println("Expected sum (unitary property): $expected_sum")
    
    # Test empirically
    n_samples = 10000
    first_mode_counts = zeros(Int, m)
    
    for _ in 1:n_samples
        # Just do first step of algorithm
        A_perm = permute_columns(A, n)
        w = abs2.(A_perm[:, 1])
        x = wsample(1:m, Weights(w))
        first_mode_counts[x] += 1
    end
    
    println("\nEmpirical first mode frequencies:")
    for i in 1:m
        emp_freq = first_mode_counts[i] / n_samples
        expected_freq = w1[i]
        println("  Mode $i: empirical=$(round(emp_freq, digits=4)), expected=$(round(expected_freq, digits=4))")
    end
end

# Run debugging
println("SIMPLE CLIFFORD DEBUG")
println("=" ^ 40)

test_2x2_system()
test_manual_calculation()
check_first_step_probabilities()

println("\nDEBUG COMPLETED")
println("=" ^ 40)