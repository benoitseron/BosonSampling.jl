sample_number_modes_occupied(params::SamplingParameters) = number_modes_occupied(BosonSampling.sample!(params))

"""

    mean_number_modes_occupied(params::SamplingParameters, n_samples::Int)
    mean_number_modes_occupied(full_dist::MultipleCounts)
    mean_number_modes_occupied(params::SamplingParameters)

Mean occupied number of modes from full distribution. When giving a number of samples, approximates the boson sampling distribution. When giving `SamplingParameters`, computes the entire distribution exactly.
"""

mean_number_modes_occupied(params::SamplingParameters, n_samples::Int) =  mean([sample_number_modes_occupied(params) for i in 1:n_samples])


function mean_number_modes_occupied(full_dist::MultipleCounts)

    @argcheck eltype(full_dist.counts) == ModeOccupation

    sum([number_modes_occupied(full_dist.counts[i]) * full_dist.proba[i] for i in 1:length(full_dist.counts)])

end

function mean_number_modes_occupied(params::SamplingParameters)

    full_dist = full_distribution(params)
    mean_number_modes_occupied(full_dist)

end
