using BosonSampling
using Plots
using Statistics
using ProgressMeter

# Include the validation functions file that contains tvd_sampled_versus_exact_distribution
include("src/boson_samplers/validating_samplers.jl")

"""
Simple example showing how to validate the Clifford sampler
against the theoretical bosonic distribution using TVD.
"""
function validate_clifford_tvd(n=4, m=4, n_events=10000)
    println("Validating Clifford sampler with TVD")
    println("Parameters: n=$n photons, m=$m modes, $n_events events")
    
    # Create a random Haar interferometer
    interf = RandHaar(m)
    
    # Input: n photons in first n modes
    input_state = Input{Bosonic}(first_modes(n, m))
    
    # Generate events using the standard sampler
    println("Generating $n_events samples...")
    events = Event{Bosonic, FockDetection}[]
    
    @showprogress for i in 1:n_events
        # Create sampling event
        ev = Event(input_state, FockSample(), interf)
        
        # Sample using the default sampler
        BosonSampling.sample!(ev)
        
        # Convert to FockDetection for probability calculation
        ev_detect = Event(input_state, 
                         FockDetection(ev.output_measurement.s), 
                         interf)
        push!(events, ev_detect)
    end
    
    # Calculate TVD between sampled and exact distributions
    tvd_value, sampled_dist, exact_dist = tvd_sampled_versus_exact_distribution(events)
    
    println("\nResults:")
    println("Total Variation Distance: $tvd_value")
    println("TVD < 0.05 (good agreement): $(tvd_value < 0.05)")
    
    # Plot the distributions
    p = bar(1:length(sampled_dist.proba), sampled_dist.proba, 
            label="Sampled", alpha=0.7, color=:blue)
    bar!(1:length(exact_dist.proba), exact_dist.proba, 
         label="Exact", alpha=0.7, color=:red)
    title!("Sampled vs Exact Distribution (TVD = $(round(tvd_value, digits=4)))")
    xlabel!("Output configuration")
    ylabel!("Probability")
    
    return tvd_value, p
end

"""
Test the Clifford sampler directly and validate it produces
correct bosonic statistics.
"""
function test_clifford_sampler_directly(n=3, m=5, n_samples=10000)
    println("\nTesting Clifford sampler directly")
    println("Parameters: n=$n photons, m=$m modes, $n_samples samples")
    
    interf = RandHaar(m)
    input_state = Input{Bosonic}(first_modes(n, m))
    
    # Collect Clifford samples
    println("Generating Clifford samples...")
    clifford_events = Event{Bosonic, FockDetection}[]
    
    @showprogress for i in 1:n_samples
        # Get sample from Clifford algorithm
        sample_vec = clifford_sampler_unoptimised(input_state, interf, 
                                                  occupancy_vector=true)
        
        # Create event with this sample
        mode_occ = ModeOccupation(sample_vec)
        ev = Event(input_state, FockDetection(mode_occ), interf)
        push!(clifford_events, ev)
    end
    
    # Validate using TVD
    tvd_value, sampled_dist, exact_dist = tvd_sampled_versus_exact_distribution(clifford_events)
    
    println("\nClifford Sampler Results:")
    println("Total Variation Distance: $tvd_value")
    println("TVD < 0.05 (good agreement): $(tvd_value < 0.05)")
    
    return tvd_value
end

"""
Run Bayesian validation: test if samples are consistent with
bosonic statistics rather than distinguishable particles.
"""
function bayesian_validation_example(n=4, m=4, n_events=1000)
    println("\nBayesian validation: Bosonic vs Distinguishable")
    println("Parameters: n=$n photons, m=$m modes, $n_events events")
    
    # Generate experimental data
    events = generate_experimental_data(n_events=n_events, n=n, m=m, 
                                      interf=RandHaar(m), TIn=Bosonic)
    
    # Set up Bayesian test
    p_bosonic = HypothesisFunction(p_B)
    p_distinguishable = HypothesisFunction(p_D)
    
    certif = Bayesian(events, p_bosonic, p_distinguishable)
    BosonSampling.certify!(certif, max_χ=Inf)
    
    println("Final confidence in bosonic hypothesis: $(certif.confidence)")
    
    # Plot evolution
    p = plot(certif.probabilities, linewidth=2,
            label="P(Bosonic | data)", 
            title="Bayesian Validation",
            xlabel="Event number", 
            ylabel="Posterior probability",
            ylim=(0, 1.05))
    hline!([0.5], linestyle=:dash, color=:gray, label="50%")
    hline!([0.99], linestyle=:dash, color=:green, label="99% confidence")
    
    return certif.confidence, p
end

# Example usage
println("=== Clifford Sampler Validation Examples ===\n")

# 1. Basic TVD validation
tvd1, plot1 = validate_clifford_tvd(3, 3, 5000)

# 2. Direct Clifford sampler test
tvd2 = test_clifford_sampler_directly(4, 6, 5000)

# 3. Bayesian validation
conf, plot2 = bayesian_validation_example(4, 4, 500)

# Display results
println("\n=== Summary ===")
println("Standard sampler TVD: $(round(tvd1, digits=4))")
println("Clifford sampler TVD: $(round(tvd2, digits=4))")
println("Bayesian confidence: $(round(conf, digits=4))")

# Show plots
display(plot1)
display(plot2)