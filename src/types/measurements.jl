### measurements ###

abstract type OutputMeasurementType end

struct FockDetection <: OutputMeasurementType

    """detection as in the standard boson sampling"""

    s::ModeOccupation
    FockDetection(s::ModeOccupation) = at_most_one_photon_per_bin(s) ? new(s) : error("more than one detector per more")

end

struct PartitionCount <: OutputMeasurementType

    """measuring the probability of getting a specific
    count for a given partition"""

    part_occupancy::PartitionOccupancy
    PartitionCounts(part_occupancy::PartitionOccupancy) = new(part_occupancy)

end

struct PartitionCountsAll <: OutputMeasurementType

    """all possible counts probabilities in a partition"""

    part::Partition
    PartitionCountsAll(part::Partition) = new(part)

end

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

        at_most_one_photon_per_bin(s) ? new(s, nothing, nothing) : error("more than one detector per more")

    end
    OutputMeasurement(s::ModeOccupation) = OutputMeasurement{FockDetection}(s::ModeOccupation)



end
