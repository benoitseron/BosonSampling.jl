### measurements ###

abstract type OutputMeasurementType end

"""
    StateMeasurement

Type trait to know which kind of state the detectors will measure, such as Fock or Gaussian.
"""
abstract type StateMeasurement end
struct FockStateMeasurement <: StateMeasurement end
struct PartitionMeasurement <: StateMeasurement end
struct GaussianStateMeasurement <: StateMeasurement end
struct CompleteDistribution <: StateMeasurement end

"""
    FockDetection(s::ModeOccupation)

Measuring the probability of getting the [`ModeOccupation`](@ref) `s` at the output.

    Fields:
        - s::ModeOccupation
"""
mutable struct FockDetection <: OutputMeasurementType
    s::ModeOccupation
    FockDetection(s::ModeOccupation) = new(s) #at_most_one_photon_per_bin(s) ? new(s) : error("more than one detector per more")
end

StateMeasurement(::Type{FockDetection}) = FockStateMeasurement()


"""
    ThresholdFockDetection(s::ThresholdModeOccupation)

Measuring the probability of getting the [`ThresholdModeOccupation`](@ref) `s` at the output.

    Fields:
        - s::ThresholdModeOccupation
"""
mutable struct ThresholdFockDetection <: OutputMeasurementType
    s::ThresholdModeOccupation
    ThresholdFockDetection(s::ThresholdModeOccupation) = new(s) #at_most_one_photon_per_bin(s) ? new(s) : error("more than one detector per more")

    ThresholdFockDetection(v::Vector{Int}) = ThresholdFockDetection(ThresholdModeOccupation(v))
end

StateMeasurement(::Type{ThresholdFockDetection}) = FockStateMeasurement()

Base.convert(::Type{FockDetection}, tmo::ThresholdFockDetection) = FockDetection(ModeOccupation(tmo.state)) 

Base.convert(::Type{ThresholdFockDetection}, mo::FockDetection) = ThresholdFockDetection(to_threshold(mo.s.state)) 


# write possible_treshold_detections for a ThresholdFockDetection
# extract the ThresholdModeOccupation and use the previous function

function possible_threshold_detections(n, state::ThresholdFockDetection; lossy = false)

    possible_threshold_detections(n,state.s, lossy = lossy)

end

# # write possible_threshold_detections for an Event
# # extract n from the input_state

# function possible_threshold_detections(ev::Event)

#     # check that the output_measurement is a ThresholdFockDetection

#     @argcheck ev.output_measurement isa ThresholdFockDetection

#     lossy = is_lossy(ev.interferometer)

#     n = ev.input_state.r.n

#     possible_threshold_detections(n,ev.output_measurement, lossy = lossy)

# end


"""
    PartitionCount(part_occupancy::PartitionOccupancy)

Measuring the probability of getting a specific count for a given partition `part_occupancy`.

    Fields:
        - part_occupancy::PartitionOccupancy
"""

struct PartitionCount <: OutputMeasurementType
    part_occupancy::PartitionOccupancy
    PartitionCount(part_occupancy::PartitionOccupancy) = new(part_occupancy)
end

StateMeasurement(::Type{PartitionCount}) = PartitionMeasurement()

"""
    PartitionCountsAll(part::Partition)

Measuring all possible counts probabilities in the partition `part`.

    Fields:
        - part::Partition
"""
struct PartitionCountsAll <: OutputMeasurementType
    part::Partition
    PartitionCountsAll(part::Partition) = new(part)
end

StateMeasurement(::Type{PartitionCountsAll}) = PartitionMeasurement()

struct OutputMeasurement{T<:OutputMeasurementType}

    # as in most of boson sampling literature, this is for detectors blind to
    # internal degrees of freedom

    s::Union{ModeOccupation,Nothing} # position of the detectors for Fock measurement

    # function OutputMeasurement{T}(s::ModeOccupation) where {T<:OutputMeasurementType}
    #     if T == FockDetection
    #         return at_most_one_photon_per_bin(s) ? new(s, nothing) : error("more than one detector per more")
    #     else
    #         return error(T, " not implemented")
    #     end
    # end

    function OutputMeasurement{FockDetection}(s::ModeOccupation)

        @warn "OutputMeasurement{FockDetection} obsolete, replace with FockDetection"

        at_most_one_photon_per_bin(s) ? new(s) : error("more than one detector per more")

    end
    OutputMeasurement(s::ModeOccupation) = OutputMeasurement{FockDetection}(s::ModeOccupation)

end

"""
    FockSample <: OutputMeasurementType

Container holding a sample from typical boson sampler.
"""
mutable struct FockSample <: OutputMeasurementType
    s::Union{ModeOccupation, Nothing}
    FockSample() = new(nothing)
    FockSample(s::Vector) = FockSample(ModeOccupation(s))
    FockSample(s::ModeOccupation) = new(s)
end

StateMeasurement(::Type{FockSample}) = FockStateMeasurement()

Base.convert(::Type{FockDetection}, fs::FockSample) = FockDetection(fs.s)



"""
    DarkCountFockSample(p)

Same as [`FockSample`](@ref) but each output mode has an extra probability `p` of giving a positive reading no matter if there is genuinely a photon.
"""
mutable struct DarkCountFockSample <: OutputMeasurementType

    s::Union{ModeOccupation, Nothing} # observed output, possibly undefined
    p::Real # probability of a dark count in each mode

    DarkCountFockSample(p::Real) = isa_probability(p) ? new(nothing, p) : error("invalid probability")
    # instantiate if no known output
end

StateMeasurement(::Type{DarkCountFockSample}) = FockStateMeasurement()

"""
    RealisticDetectorsFockSample(p_dark::Real, p_no_count::Real)

Same as [`DarkCountFockSample`](@ref) with the added possibility that no reading is observed although there is a photon. This same probability also removes dark counts (first a dark count sample is generated then readings are discarded with probability `p_no_count`).
"""
mutable struct RealisticDetectorsFockSample <: OutputMeasurementType

    s::Union{ModeOccupation, Nothing} # observed output, possibly undefined
    p_dark::Real # probability of a dark count in each mode
    p_no_count::Real # probability that there is a photon but it is not seen

    # instantiate if no known output
    RealisticDetectorsFockSample(p_dark::Real, p_no_count::Real) = begin
        if isa_probability(p_dark) && isa_probability(p_no_count)
            new(nothing, p_dark, p_no_count)
        else
            error("invalid probability")
        end
    end

end

StateMeasurement(::Type{RealisticDetectorsFockSample}) = FockStateMeasurement()




"""
    PartitionSample <: OutputMeasurementType

Container holding a sample from `Partition` photon count.
"""
mutable struct PartitionSample <: OutputMeasurementType
    part_occ::Union{PartitionOccupancy, Nothing}
    PartitionSample() = new(nothing)
    PartitionSample(p::PartitionOccupancy) = new(p)
end


StateMeasurement(::Type{PartitionSample}) = PartitionMeasurement()


mutable struct ThresholdDetection <: OutputMeasurementType
    s::Union{Vector{Int64}, Nothing}
    ThresholdDetection() = new(nothing)
    ThresholdDetection(s::Vector{Int64}) = new(s)
end

"""
	MultipleCounts()
	MultipleCounts(counts, proba)

Holds something like the photon counting probabilities with their respective
probability (in order to use them as a single observation). Can be declared
empty as a placeholder.

    Fields:
    - counts::Union{Nothing, Vector{ModeOccupation}, Vector{PartitionOccupancy}, Vector{ThresholdModeOccupation}},
    - proba::Union{Nothing,Vector{Real},Vector{Int}}

!!! Warning:

    The naming of `proba` and `counts` was done at a much earlier stage of the project. Understand `counts` as detector readings. `proba` can hold either probabilites or also number of times an event was observed.


"""
mutable struct MultipleCounts

	counts::Union{Nothing, Vector{ModeOccupation}, Vector{PartitionOccupancy}, Vector{ThresholdModeOccupation}}
	proba::Union{Nothing,Vector{TReal},Vector{TInt}} where {TReal <: Real, TInt <: Int}

	MultipleCounts() = new(nothing,nothing)
	MultipleCounts(counts, proba) = new(counts,proba)

end

function initialise_to_empty_vectors!(mc::MultipleCounts, type_proba, type_counts)

	mc.proba = Vector{type_proba}()
	mc.counts = Vector{type_counts}()

end

Base.show(io::IO, pb::MultipleCounts) = begin

	if pb.proba == nothing
		println(io, "Empty MultipleCounts")
	else

        #check if eltype of proba is an Int, and if so set variable `detector_counts` to true
        detector_counts = eltype(pb.proba) <: Int

		for i in 1:length(pb.proba)

			println(io, "output: \n")
			println(io, pb.counts[i])
            if detector_counts
                println(io, "counts = $(pb.proba[i])")
            else
                println(io, "p = $(pb.proba[i])")
            end
			
			println(io, "--------------------------------------")
		end
	end

end


"""
	to_threshold(mc::MultipleCounts)

Transforms a `MultipleCounts` into the equivalent for threshold detectors.
"""
function to_threshold(mc::MultipleCounts)

    count_proba = Dict()

    for (count, proba) in zip(mc.counts, mc.proba)
        new_count = to_threshold(count)

        if new_count in keys(count_proba)
            count_proba[new_count] += proba
        else
            count_proba[new_count] = proba
        end

    end

    # println("######")
    # @show count_proba

    counts = Vector{typeof(mc.counts[1])}()
    probas = Vector{typeof(mc.proba[1])}()

    for key in keys(count_proba)

        push!(counts, key)
        push!(probas, count_proba[key])

    end

    # @show counts
    # @show probas

    MultipleCounts(counts, probas)

end

function to_threshold!(mc::MultipleCounts)

    mc = to_threshold(mc)

end


"""
    BosonSamplingDistribution <: OutputMeasurementType

Container holding the entire boson sampling distribution for a given type of parameters, input, etc.
"""
mutable struct BosonSamplingDistribution <: OutputMeasurementType

    mc::Union{MultipleCounts, Nothing}
	BosonSamplingDistribution() = new(nothing)
	BosonSamplingDistribution(mc) = new(mc)

end

StateMeasurement(::Type{BosonSamplingDistribution}) = CompleteDistribution()
