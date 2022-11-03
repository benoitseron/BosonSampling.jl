"""
    mutable struct ThresholdModeOccupation

Holds threshold detector clicks. Example

    ThresholdModeOccupation(ModeList([1,2,4], 4))

"""
@auto_hash_equals mutable struct ThresholdModeOccupation

    m::Int
    clicks::Vector{Int}

    function ThresholdModeOccupation(ml::ModeList)

        clicks = convert(ModeOccupation, ml).state

        if !all(clicks[:] .>= 0)
            error("negative mode clicks")
        elseif !all(clicks[:] .<= 1)
            error("clicks can be at most one")
        else
            new(ml.m, clicks)
        end
    end

end

# example ThresholdModeOccupation(ModeList([1,2,4], 4))


"""
    @with_kw mutable struct OneLoopData

Contains experimental data telling everything necessary about an experiment running a one loop boson sampler.
"""
@with_kw mutable struct OneLoopData

    params::LoopSamplingParameters
    samples::Vector{ThresholdModeOccupation}
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
