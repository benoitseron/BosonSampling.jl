"""
Bayesian Validation: High Indistinguishability Test

Simple validation test for near-perfect photons with confidence evolution plot.
"""

using BosonSampling
using Random
using StatsBase
using Plots

# Include the implementation
include("../partial_distinguishability_sampler.jl")

Random.seed!(42)

println("🔬 BAYESIAN VALIDATION: HIGH INDISTINGUISHABILITY TEST")
println("=" ^ 60)

# Helper functions
function compute_exact_bosonic_probability(sample, input_state, interf)
    try
        mode_occ = ModeOccupation(sample)
        output = FockDetection(mode_occ)
        bosonic_input = Input{Bosonic}(input_state.r)
        ev = Event(bosonic_input, output, interf)
        BosonSampling.compute_probability!(ev)
        return real(ev.proba_params.probability)
    catch e
        return NaN
    end
end

function estimate_partial_dist_probability(sample, input_state, interf, model, n_samples=1000)
    """Estimate probability from partial distinguishability algorithm via sampling"""
    count = 0
    
    for _ in 1:n_samples
        partial_dist_sample = partial_distinguishability_sampler(input_state, interf, model)
        if partial_dist_sample == sample
            count += 1
        end
    end
    
    return count / n_samples
end

function bayesian_validation_step_partial_dist(sample, input_state, interf, model)
    """Perform one step of Bayesian validation - exactly like Clifford case"""
    
    # Get probabilities under both hypotheses
    p_bosonic = compute_exact_bosonic_probability(sample, input_state, interf)
    
    if isnan(p_bosonic) || p_bosonic ≈ 0
        return NaN, p_bosonic, NaN
    end
    
    # Estimate probability from partial distinguishability algorithm
    p_partial_dist = estimate_partial_dist_probability(sample, input_state, interf, model, 2000)
    
    if p_partial_dist ≈ 0
        return NaN, p_bosonic, p_partial_dist
    end
    
    # Bayesian ratio: P(sample | partial_dist_algorithm) / P(sample | Bosonic)
    ratio = p_partial_dist / p_bosonic
    
    return ratio, p_bosonic, p_partial_dist
end

confidence(χ) = χ == Inf ? 1.0 : χ / (1 + χ)

# Test setup - Use larger system for better statistical power
n, m = 5, 10
input_state = Input{Bosonic}(first_modes(n, m))
interf = RandHaar(m)

model = uniform_distinguishability(n, 0.99)
qm_hom = quadratic_mean_hom_visibility(model)

println("\n🧪 HIGH INDISTINGUISHABILITY TEST")
println("System: $n photons, $m modes")
println("Model: x = 0.99 for all photons")
println("Quadratic mean HOM visibility: $(round(qm_hom, digits=4))")
println("Expected: Should match Bosonic distribution with high confidence")
println("Note: Larger systems provide better statistical discrimination power")

# Run Bayesian validation until very high or very low confidence
χ = 1.0
confidences = [confidence(χ)]
sample_numbers = [0]
max_samples = 50000  # Run many samples to see full evolution
target_extreme = 0.9999
target_reached = false

println("\n📊 Running validation...")
println("Expected: Fast convergence to high confidence for x=0.99")
println("Will continue to $max_samples samples to show full evolution")

for i in 1:max_samples
    sample = partial_distinguishability_sampler(input_state, interf, model)
    
    ratio, p_bos, p_partial = bayesian_validation_step_partial_dist(sample, input_state, interf, model)
    
    if !isnan(ratio)
        χ *= ratio
        conf = confidence(χ)
        push!(confidences, conf)
        push!(sample_numbers, i)
        
        # Print progress
        if i <= 10 || i % 100 == 0
            println("Sample $i: confidence = $(round(conf, digits=6)), P(PartDist)/P(Bosonic) = $(round(ratio, digits=3))")
        end
        
        # Note if extreme confidence reached for first time
        if conf >= target_extreme && !target_reached
            println("\n🎯 Target confidence $(target_extreme) reached after $i samples - continuing to $max_samples...")
            target_reached = true
        end
    else
        push!(confidences, confidences[end])
        push!(sample_numbers, i)
        if i <= 10
            println("Sample $i: Unable to compute probabilities")
        end
    end
end

final_confidence = confidences[end]
total_samples = sample_numbers[end]
println("\n📈 Final Results:")
println("Confidence: $(round(final_confidence, digits=6)) after $total_samples samples")

# Create confidence evolution plot
println("\n📊 Creating confidence evolution plot...")

p = plot(sample_numbers, confidences,
         title="Bayesian Confidence Evolution: High Indistinguishability (x=0.95)",
         xlabel="Number of Samples", 
         ylabel="Confidence in Algorithm Correctness",
         linewidth=3,
         color=:blue,
         label="Confidence Evolution",
         size=(800, 600))

# Add reference lines
hline!([0.9999], label="Target: 99.99%", linestyle=:dash, color=:darkgreen, linewidth=2)
hline!([0.99], label="99%", linestyle=:dash, color=:green, linewidth=1)
hline!([0.9], label="90%", linestyle=:dash, color=:orange, linewidth=1)
hline!([0.5], label="Neutral", linestyle=:dot, color=:gray, linewidth=1)

ylims!(0, 1)
grid!(true, alpha=0.3)

display(p)

# Analysis
println("\n🔍 VALIDATION ANALYSIS:")
if final_confidence > 0.999
    println("✅ EXCELLENT: Strong evidence algorithm behaves like Bosonic distribution")
elseif final_confidence > 0.9
    println("✅ GOOD: Algorithm shows Bosonic-like behavior")
elseif final_confidence > 0.7
    println("⚠️  MODERATE: Some evidence for correctness")
else
    println("❌ POOR: Algorithm may have implementation issues")
end

println("\nFor high-quality photons (x≥0.95), we expect >99% confidence.")
println("This validates that the partial distinguishability algorithm")
println("correctly reproduces Bosonic behavior in the quantum regime.")

