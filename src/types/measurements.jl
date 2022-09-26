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

"""
    FockDetection(s::ModeOccupation)

Measuring the probability of getting the [`ModeOccupation`](@ref) `s` at the output.

    Fields:
        - s::ModeOccupation
"""
struct FockDetection <: OutputMeasurementType
    s::ModeOccupation
    FockDetection(s::ModeOccupation) = new(s) #at_most_one_photon_per_bin(s) ? new(s) : error("more than one detector per more")
end

StateMeasurement(::Type{FockDetection}) = FockStateMeasurement()

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


mutable struct TresholdDetection <: OutputMeasurementType
    s::Union{Vector{Int64}, Nothing}
    TresholdDetection() = new(nothing)
    TresholdDetection(s::Vector{Int64}) = new(s)
end

