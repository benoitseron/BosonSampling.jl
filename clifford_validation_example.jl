using BosonSampling
using Plots
using ProgressMeter
using Statistics
using StatsBase
using Random

# Example 1: TVD Validation of Clifford Sampler
# This example shows how to validate that the Clifford sampler produces samples
# following the correct bosonic distribution by computing the Total Variation Distance

function validate_clifford_with_tvd(n::Int, m::Int, n_events::Int=10000; interf_type="RandHaar")
    println("\n=== TVD Validation of Clifford Sampler ===")
    println("Parameters: n=$n photons, m=$m modes, $n_events events")
    
    # Create interferometer
    if interf_type == "RandHaar"
        interf = RandHaar(m)
    elseif interf_type == "Fourier"
        interf = Fourier(m)
    else
        error("Unknown interferometer type")
    end
    
    # Input state: n photons in first n modes
    input_state = Input{Bosonic}(first_modes(n, m))
    
    # Generate samples using the standard sampler (for experimental data)
    println("Generating samples...")
    events = Event{Bosonic, FockDetection}[]
    
    @showprogress for i in 1:n_events
        # Create event with FockSample output
        ev = Event(input_state, FockSample(), interf)
        BosonSampling.sample!(ev)
        
        # Convert to FockDetection for probability computation
        ev_detection = convert(Event{Bosonic, FockDetection}, ev)
        push!(events, ev_detection)
    end
    
    # Compute TVD between sampled and exact distribution
    tvd_value, sampled_dist, exact_dist = tvd_sampled_versus_exact_distribution(events)
    
    println("Total Variation Distance: $tvd_value")
    println("TVD < 0.05? $(tvd_value < 0.05)")
    
    # Plot comparison
    p = bar(sampled_dist.proba, label="Sampled", alpha=0.5, 
            title="Sampled vs Exact Distribution (TVD=$(@sprintf("%.4f", tvd_value)))")
    bar!(p, exact_dist.proba, label="Exact", alpha=0.5)
    xlabel!(p, "Output configuration index")
    ylabel!(p, "Probability")
    
    return tvd_value, p, sampled_dist, exact_dist
end

# Example 2: Bayesian Validation - Bosonic vs Distinguishable
# This validates that the samples are consistent with bosonic statistics
# rather than distinguishable particle statistics

function validate_bosonic_vs_distinguishable(n::Int, m::Int, n_events::Int=1000)
    println("\n=== Bayesian Validation: Bosonic vs Distinguishable ===")
    println("Parameters: n=$n photons, m=$m modes, $n_events events")
    
    # Create interferometer and input
    interf = RandHaar(m)
    input_state = Input{Bosonic}(first_modes(n, m))
    
    # Generate experimental data (bosonic sampling)
    events = generate_experimental_data(n_events=n_events, n=n, m=m, 
                                      interf=interf, TIn=Bosonic)
    
    # Set up hypothesis functions
    p_bosonic = HypothesisFunction(p_B)  # Bosonic hypothesis
    p_distinguishable = HypothesisFunction(p_D)  # Distinguishable hypothesis
    
    # Create Bayesian certifier
    certif = Bayesian(events, p_bosonic, p_distinguishable)
    
    # Run certification
    BosonSampling.certify!(certif, max_χ=Inf)
    
    println("Confidence in bosonic hypothesis: $(certif.confidence)")
    println("Bosonic hypothesis strongly favored? $(certif.confidence > 0.99)")
    
    # Plot probability evolution
    p = plot(certif.probabilities, label="P(Bosonic|data)", 
            title="Bayesian Validation: Bosonic vs Distinguishable",
            xlabel="Event number", ylabel="Posterior probability",
            ylim=(0, 1.05), legend=:bottomright)
    hline!(p, [0.5], label="Equal probability", linestyle=:dash, color=:gray)
    hline!(p, [0.99], label="99% confidence", linestyle=:dash, color=:green)
    
    return certif.confidence, p, certif
end

# Example 3: Validate Clifford sampler specifically
# This creates samples using the Clifford algorithm and validates them

function validate_clifford_algorithm(n::Int, m::Int, n_samples::Int=10000)
    println("\n=== Direct Clifford Algorithm Validation ===")
    println("Parameters: n=$n photons, m=$m modes, $n_samples samples")
    
    # Create interferometer
    interf = RandHaar(m)
    input_state = Input{Bosonic}(first_modes(n, m))
    
    # Generate samples using Clifford algorithm
    println("Generating Clifford samples...")
    clifford_samples = []
    
    @showprogress for i in 1:n_samples
        # Use the unoptimised Clifford sampler (the optimized one is disabled)
        sample = clifford_sampler_unoptimised(input_state, interf, occupancy_vector=true)
        push!(clifford_samples, ModeOccupation(sample))
    end
    
    # Count occurrences
    sample_counts = countmap(clifford_samples)
    
    # Compute exact probabilities for comparison
    println("Computing exact probabilities...")
    exact_probs = Dict{ModeOccupation, Float64}()
    
    all_configs = all_mode_configurations(n, m, only_photon_number_conserving=true)
    
    @showprogress for config in all_configs
        mode_occ = ModeOccupation(config)
        o = FockDetection(mode_occ)
        ev = Event(input_state, o, interf)
        compute_probability!(ev)
        exact_probs[mode_occ] = ev.proba_params.probability
    end
    
    # Compare distributions using TVD
    sampled_probs = Float64[]
    exact_probs_ordered = Float64[]
    
    for mode_occ in keys(exact_probs)
        push!(exact_probs_ordered, exact_probs[mode_occ])
        push!(sampled_probs, get(sample_counts, mode_occ, 0) / n_samples)
    end
    
    tvd_value = tvd(sampled_probs, exact_probs_ordered)
    
    println("TVD between Clifford samples and exact: $tvd_value")
    println("TVD < 0.05? $(tvd_value < 0.05)")
    
    # Plot top probabilities
    sorted_exact = sort(collect(exact_probs), by=x->x[2], rev=true)
    top_n = min(20, length(sorted_exact))
    
    exact_top = [p[2] for p in sorted_exact[1:top_n]]
    sampled_top = [get(sample_counts, p[1], 0)/n_samples for p in sorted_exact[1:top_n]]
    
    p = bar(1:top_n, exact_top, label="Exact", alpha=0.5,
            title="Top $top_n Output Probabilities (TVD=$(@sprintf("%.4f", tvd_value)))")
    bar!(p, 1:top_n, sampled_top, label="Clifford sampled", alpha=0.5)
    xlabel!(p, "Output rank")
    ylabel!(p, "Probability")
    
    return tvd_value, p, sample_counts, exact_probs
end

# Example 4: Partition-based validation
# This validates using partition statistics which are more experimentally accessible

function validate_with_partitions(n::Int, m::Int, n_events::Int=1000, n_subsets::Int=2)
    println("\n=== Partition-based Bayesian Validation ===")
    println("Parameters: n=$n, m=$m, $n_events events, $n_subsets partitions")
    
    # Generate events
    interf = RandHaar(m)
    events = generate_experimental_data(n_events=n_events, n=n, m=m, 
                                      interf=interf, TIn=Bosonic)
    
    # Create equilibrated partition
    part = equilibrated_partition(m, n_subsets)
    
    # Create Bayesian partition certifier
    certif = BayesianPartition(events, Bosonic(), Distinguishable(), part)
    
    # Certify
    certify!(certif, max_χ=Inf, min_χ=0.0001)
    
    println("Partition-based confidence: $(certif.confidence)")
    
    p = plot(certif.probabilities, label="P(Bosonic|partition data)",
            title="Partition-based Validation",
            xlabel="Event number", ylabel="Posterior probability", 
            ylim=(0, 1.05))
    hline!(p, [0.99], label="99% confidence", linestyle=:dash, color=:green)
    
    return certif.confidence, p, certif
end

# Run all validation examples
function run_all_validations()
    Random.seed!(42)  # For reproducibility
    
    # Small system for exact validation
    n, m = 4, 4
    
    # 1. TVD validation
    tvd_val, p1, _, _ = validate_clifford_with_tvd(n, m, 10000)
    
    # 2. Bayesian validation
    conf, p2, _ = validate_bosonic_vs_distinguishable(n, m, 1000)
    
    # 3. Direct Clifford validation
    tvd_cliff, p3, _, _ = validate_clifford_algorithm(n, m, 10000)
    
    # 4. Partition validation
    conf_part, p4, _ = validate_with_partitions(n, m, 1000, 2)
    
    # Combine plots
    plot(p1, p2, p3, p4, layout=(2,2), size=(1000, 800))
end

# Example usage:
# run_all_validations()

# Individual examples:
# validate_clifford_with_tvd(3, 3, 5000)
# validate_bosonic_vs_distinguishable(4, 4, 500)
# validate_clifford_algorithm(3, 5, 5000)
# validate_with_partitions(5, 10, 1000, 3)