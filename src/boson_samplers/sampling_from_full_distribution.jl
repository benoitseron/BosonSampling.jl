using Revise
using BosonSampling
using Plots
using Distributions

### generating the full distribution ###

n = 3
m = 12

interf = RandHaar(m)
TIn = Bosonic
input_state = Input{TIn}(first_modes(n,m))
i = input_state
output_measurement = ThresholdFockSample() # this is a way to encode the type of detector


if isa(output_measurement, FockSample)

    o = BosonSamplingDistribution()
    # println("number of computations: $(binomial(n+m-1,n)) (no loss)")

elseif isa(output_measurement, ThresholdFockSample)

    o = BosonSamplingThresholdDistribution()
    
else

    error("$(typeof(output_measurement)) not implemented")

end

ev_full_distribution_threshold = Event(i, o, interf)
compute_probability!(ev_full_distribution_threshold)

mc_full_proba = ev_full_distribution_threshold.proba_params.probability

bar(mc_full_proba.proba, label = "full", alpha = 0.5)
ylabel!("probability")
xlabel!("outcome number")

### getting a sample ###

n_samples = 100000
samples = wsample(mc_full_proba.counts, mc_full_proba.proba, n_samples) 
samples_clicks = [sample.state for sample in samples] 