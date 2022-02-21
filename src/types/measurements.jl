
### measurements ###

abstract type OutputMeasurementType end

struct FockDetection <: OutputMeasurementType
end

struct PartitionCounts <: OutputMeasurementType
end

struct OutputMeasurement{T<:OutputMeasurementType}

    # as in most of boson sampling literature, this is for detectors blind to
    # internal degrees of freedom

    s::ModeOccupation # position of the detectors for Fock measurement
    part::Partition # partition detectors for PartitionCounts

    function OutputMeasurement{T}(s::ModeOccupation) where {T<:OutputMeasurementType}
        if T == FockDetection
            return at_most_one_photon_per_bin(s) ? new(s, nothing) : error("more than one detector per more")
        else
            return error(T, " not implemented")
        end
    end

    function OutputMeasurement{T}(part::Partition) where {T<:OutputMeasurementType}
        if T == FockDetection
            return at_most_one_photon_per_bin(s) ? new(s, nothing) : error("more than one detector per more")
        else
            return error(T, " not implemented")
        end
    end

end
