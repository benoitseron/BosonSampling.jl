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

@with_kw mutable struct OneLoopData

    params::LoopSamplingParameters
    samples::Vector{ThresholdModeOccupation}
    extra_info

end
