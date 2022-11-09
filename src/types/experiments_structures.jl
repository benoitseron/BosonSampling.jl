

"""
    @with_kw mutable struct OneLoopData

Contains experimental data telling everything necessary about an experiment running a one loop boson sampler.

For high number of detections, using a `Vector{ThresholdModeOccupation}` is inefficient, so `samples` can also be a `MultipleCounts`.
"""
@with_kw mutable struct OneLoopData

    params::LoopSamplingParameters
    samples::Union{Vector{ThresholdModeOccupation}, MultipleCounts}
    date::DateTime = now()
    name::String = string(date)
    extra_info

end

"""
    save(data::OneLoopData; path_to_file::String = "data/one_loop/")

Saves a OneLoopData.
"""
function JLD.save(data::OneLoopData; path_to_file::String = "data/one_loop/")

    save(path_to_file * "$(data.name).jld", "data", data)

end
