
### measurements ###

abstract type OutputMeasurementType end

struct FockDetection <: OutputMeasurementType
end

struct PartitionCounts <: OutputMeasurementType
end

struct OutputMeasurement{T<:OutputMeasurementType}

    # as in most of boson sampling literature, this is for detectors blind to
    # internal degrees of freedom

    s::Union{ModeOccupation,Nothing} # position of the detectors for Fock measurement
    part_occupancy::Union{PartitionOccupancy,Nothing} # partition detectors for PartitionCounts

    # function OutputMeasurement{T}(s::ModeOccupation) where {T<:OutputMeasurementType}
    #     if T == FockDetection
    #         return at_most_one_photon_per_bin(s) ? new(s, nothing) : error("more than one detector per more")
    #     else
    #         return error(T, " not implemented")
    #     end
    # end

    OutputMeasurement{FockDetection}(s::ModeOccupation) = at_most_one_photon_per_bin(s) ? new(s, nothing) : error("more than one detector per more")
    OutputMeasurement(s::ModeOccupation) = OutputMeasurement{FockDetection}(s::ModeOccupation)

    OutputMeasurement{PartitionCounts}(part_occupancy::PartitionOccupancy) = new(nothing, part_occupancy)
    OutputMeasurement(part_occupancy::PartitionOccupancy) = OutputMeasurement{PartitionCounts}(part_occupancy::PartitionOccupancy)


end
