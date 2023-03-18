abstract type ExperimentalData end

"""
    @with_kw mutable struct OneLoopData

Contains experimental data telling everything necessary about an experiment running a one loop boson sampler.

For high number of detections, using a `Vector{ThresholdModeOccupation}` is inefficient, so `samples` can also be a `MultipleCounts`.
"""
@with_kw mutable struct OneLoopData <: ExperimentalData

    params::LoopSamplingParameters
    samples::Union{Vector{ThresholdModeOccupation}, MultipleCounts}
    sample_type::Type{T} where {T<:Union{Nothing, ModeOccupation, ThresholdModeOccupation}} = ThresholdModeOccupation
    date::DateTime = now()
    name::String = string(date)
    certification_data::Union{Nothing, Certifier} = nothing
    extra_info

end


function build_loop!(data::OneLoopData)
    build_loop!(data.params)
end

"""

    set_uniform_losses!(η_lines, η_bs, params)

Set the losses to be uniform for all the lines and beamsplitters. Rebuild the loop. Used for estimation of the losses.

"""
function set_uniform_losses!(η_lines, η_bs, params)
    params.η_loss_lines = η_lines * ones(params.m)
    params.η_loss_bs = η_bs * ones(params.m-1)

    build_loop!(params)
end

"""
    save(data::OneLoopData; path_to_file::String = "data/one_loop/", recompile_interferometer = true)

Saves a OneLoopData. Recompiles the interferometer to make sure it is the right one.
"""
function JLD.save(data::OneLoopData; path_to_file::String = "data/one_loop/")
 
    build_loop!(data) # making sure the interferometer is up to date

    if data.certification_data != nothing
        @warn "Error in saving the certification data. Setting it to nothing. Saving is not possible for the moment. Certifier contains functions, which cannot be saved by JLD because of an internal bug, see https://github.com/JuliaIO/JLD.jl/issues/57. I think that it is possible to fix this by removing the HypothesisFunction's, but it is quite annoying... "
    end

    data.certification_data = nothing

    save(path_to_file * "$(data.name).jld", "data", data)

end
