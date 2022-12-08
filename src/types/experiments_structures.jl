

"""
    @with_kw mutable struct OneLoopData

Contains experimental data telling everything necessary about an experiment running a one loop boson sampler.

For high number of detections, using a `Vector{ThresholdModeOccupation}` is inefficient, so `samples` can also be a `MultipleCounts`.
"""
@with_kw mutable struct OneLoopData

    params::LoopSamplingParameters
    samples::Union{Vector{ThresholdModeOccupation}, MultipleCounts}
    sample_type::Union{Nothing, ModeOccupation, ThresholdModeOccupation} = ThresholdModeOccupation
    date::DateTime = now()
    name::String = string(date)
    extra_info

end

"""
    save(data::OneLoopData; path_to_file::String = "data/one_loop/", recompile_interferometer = true)

Saves a OneLoopData. Recompiles the interferometer to make sure it is the right one.
"""
function JLD.save(data::OneLoopData; path_to_file::String = "data/one_loop/")
 
    recompile_interferometer ? build_loop!(data) : nothing
    save(path_to_file * "$(data.name).jld", "data", data)

end
