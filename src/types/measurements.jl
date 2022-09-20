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
    PartitionSample <: OutputMeasurementType

Container holding a sample from `Partition` photon count.
"""
mutable struct PartitionSample <: OutputMeasurementType
    part_occ::Union{PartitionOccupancy, Nothing}

    PartitionSample() = new(nothing)
    PartitionSample(p::PartitionOccupancy) = new(p)
end

StateMeasurement(::Type{PartitionSample}) = PartitionMeasurement()
