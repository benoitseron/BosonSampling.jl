loss_amplitude_to_transmission_amplitude(loss::Real) = sqrt(1-loss^2)

"""
    virtual_interferometer_uniform_loss(real_interf::Matrix, η)
    virtual_interferometer_uniform_loss(real_interf::Interferometer, η)

Simulates a simple, uniformly lossy interferometer: take a 2m*2m interferometer and introduce beam splitters in front with transmittance `η`. Returns the corresponding virtual interferometer.
"""
function virtual_interferometer_uniform_loss(U::Matrix, η)

    # define a (2m) * (2m) virtual interferometer V

    # connect the remaining branches of the input beam splitters
    # (of course this could be made different but this is the simplest,
    # no thought case)

    m = size(U,1)
    V = Matrix{eltype(U)}(I, (2m, 2m))

    # loss beam splitters
    for i in 1:m
        V *= beam_splitter_modes(in_up = i, in_down = i+m, out_up = i, out_down = m+i, transmission_amplitude = η, n = 2m)
    end

    U_large = Matrix{eltype(U)}(I, (2m, 2m))
    U_large[1:m,1:m] = U
    V *= U_large

    #V = V[1:2m, 1:2m] # disregard virtual mode to connect second input branch of beam splitters
    ############ this must clearly not be good for unitarity

    V = copy(transpose(V))

    V

end

function virtual_interferometer_uniform_loss(real_interf::Interferometer, η)
    virtual_interferometer_uniform_loss(real_interf.U)
end


"""
    virtual_interferometer_general_loss(V::Matrix, W::Matrix, η::Vector{Real})

Generic lossy interferometer composed of two unitary matrices `V`, `W` forming the
physical matrix `U` = `V*W`. In between `V` and `W` are sandwiched a diagonal array
of beam splitters, with transmissivity `η` (`m`-dimensional vector corresponding
to the transmissivity of each layer) : `U_total` = `V*Diag(η)*W`. This generates `U_total`.
"""
function virtual_interferometer_general_loss(V::Matrix, W::Matrix, η::Vector{Real})

    # see GeneralLossInterferometer for info

    m = size(V,1)
    U_virtual = Matrix{eltype(V)}(I, (2m, 2m))
    U_virtual[1:m, 1:m] = W

    # loss beam splitters
    for i in 1:m
        U_virtual *= beam_splitter_modes(in_up = i, in_down = i+m, out_up = i, out_down = m+i, transmission_amplitude = η[i], n = 2m)
    end

    U_virtual *= V

    U_virtual = copy(transpose(U_virtual))

    U_virtual

end


"""
    to_lossy(s::Subset)

Transforms a subset of size `m` into a `2m` one with the last `m` modes being a empty (the environment modes).
"""
function to_lossy(s::Subset)
    subset = s.subset
    new_subset = zeros(Int, 2*length(subset))
    new_subset[1:length(subset)] .= subset
    Subset(new_subset)
end

"""
    to_lossy(part::Partition)

Transforms a partition of size `m` into a `2m` one with the last `m` modes being a subset (the environment modes).
"""
function to_lossy(part::Partition)

    n = part.subsets[1].n
    m = part.subsets[1].m

    environment = Subset(last_modes(m,2m))

    new_subsets =  Vector{Subset}()

    for subset in part.subsets
        push!(new_subsets, to_lossy(subset))
    end

    push!(new_subsets, environment)

    Partition(new_subsets)

end

"""
    to_lossy(mo::ModeOccupation)

Transforms a `ModeOccupation` into the same with extra padding to account for the environment modes.
"""
function to_lossy(mo::ModeOccupation)

    cat(mo,zeros(mo))

end

function to_lossy(state::Vector{Int})

    vcat(state,zeros(eltype(state), length(state)))

end

function to_lossy(interf::Interferometer)

    if isa(interf, LossyInterferometer)
        error("$interf is already a LossyInterferometer")
    else

        U = interf.U
        m = interf.m
        U_lossy = Matrix{eltype(U)}(I, 2 .* size(U))
        U_lossy[1:m,1:m] = U

        return UserDefinedLossyInterferometer(U_lossy)
    end

end

function to_lossy(i::Input{T}) where {T<:InputType}

    Input{T}(to_lossy(i.r), i.n, 2*i.m, i.G, i.distinguishability_param)

end

function to_lossy(o::OutputMeasurementType)

    # if StateMeasurement(typeof(output_measurement)) == FockStateMeasurement
    if typeof(o) == FockDetection
        FockDetection(to_lossy(o.s))

    else

        error("not implemented")

    end

end

"""
    lossy_target_modes(target_modes::Vector{Int})

Converts a vector of mode occupation into the same concatenated twice. This allows for modes occupied by circuit elements to have their loss mode attributed. For instance, a LossyLine targeting mode 1 with m=2 has

    target_modes = [1,0]
    lossy_target_modes(target_modes) = [1,0,1,0]

"""
function lossy_target_modes(target_modes::Vector{Int})

    vcat(target_modes, target_modes)

end




"""
    isa_transmissitivity(η::Real)
    isa_transmissitivity(η::Vector{Real})

Asserts that η is a valid transmissivity.
"""
function isa_transmissitivity(η::Real)
    (0<= η && η <= 1)
end

isa_transmissitivity(η::Vector{Real}) = isa_transmissitivity.(η)

"""
    UniformLossInterferometer <: LossyInterferometer
    UniformLossInterferometer(η::Real, U_physical::Matrix)
    UniformLossInterferometer(η::Real, U_physical::Interferometer)
    UniformLossInterferometer(m::Int, η::Real)

Simulates a simple, uniformly lossy interferometer: take a 2m*2m interferometer and introduce beam splitters in front with transmittance `η`. Returns the corresponding interferometer as a separate type.

The last form, `UniformLossInterferometer(m::Int, η::Real)` samples from a Haar random unitary.
"""
struct UniformLossInterferometer <: Interferometer
    m_real::Int
    m::Int
    η::Real #transmissivity of the upfront beamsplitters
    U_physical::Union{Matrix{Complex},Nothing} # physical matrix
    U::Union{Matrix{Complex},Nothing} # virtual 2m*2m interferometer

    function UniformLossInterferometer(η::Real, U_physical::Matrix)
        if isa_transmissitivity(η)
            U_virtual = virtual_interferometer_uniform_loss(U_physical, η)
            new(size(U_physical,1), size(U_virtual,1), η, U_physical, U_virtual)
        else
            error("incorrect η")
        end
    end

    UniformLossInterferometer(η::Real, U_physical::Interferometer) = UniformLossInterferometer(η, U_physical.U)

    UniformLossInterferometer(η::Real, m::Int) = UniformLossInterferometer(η, RandHaar(m))

end

"""
    GeneralLossInterferometer <: LossyInterferometer

Generic lossy interferometer composed of two unitary matrices `V`, `W` forming the
physical matrix `U` = `V*W`. In between `V` and `W` are sandwiched a diagonal array
of beam splitters, with transmissivity `η` (`m`-dimensional vector corresponding
to the transmissivity of each layer) : `U_total` = `V*Diag(η)*W`
"""
struct GeneralLossInterferometer <: Interferometer
    m_real::Int
    m::Int
    η::Vector{Real} #transmissivity of the upfront beamsplitters
    U_physical::Union{Matrix{Complex},Nothing} # physical matrix
    V::Union{Matrix{Complex},Nothing}
    W::Union{Matrix{Complex},Nothing}
    U::Union{Matrix{Complex},Nothing} # virtual 2m*2m interferometer

    function GeneralLossInterferometer(η::Vector{Real}, V::Matrix, W::Matrix)
        if isa_transmissitivity(η)
            U_virtual = virtual_interferometer_general_loss(V,W, η)
            new(size(U_physical,1), size(U_virtual,1), η, U_physical,V,W, U_virtual)
        else
            error("invalid η")
        end
    end

end


struct UserDefinedLossyInterferometer <: Interferometer
    m_real::Int
    m::Int
    η::Union{Real, Vector{Real}, Nothing} #transmissivity of the upfront beamsplitters
    U_physical::Union{Matrix{Complex},Nothing} # physical matrix
    U::Union{Matrix{Complex},Nothing} # virtual 2m*2m interferometer

    function UserDefinedLossyInterferometer(U::Matrix)

        m_real = Int(size(U,1)/2)
        new(m_real, 2*m_real, nothing, U[1:m_real,1:m_real], U)

    end

end


"""
    sort_by_lost_photons(pb::MultipleCounts)

Outputs a (n+1) sized array of MultipleCounts containing 0,...,n lost photons.
"""
function sort_by_lost_photons(pb::MultipleCounts)

    # number of lost photons is the number of photons in the last subset

    n = pb.counts[1].n
    sorted_array = [MultipleCounts() for i in 0:n]
    initialise_to_empty_vectors!.(sorted_array, Real, PartitionOccupancy)

    for (p, count) in zip(pb.proba, pb.counts)

        n_lost = count.counts.state[end]
        push!(sorted_array[n_lost+1].proba, p)
        push!(sorted_array[n_lost+1].counts, count)

    end

    sorted_array

end


"""
    tvd_k_lost_photons(k, pb_sorted, pd_sorted)

Gives the tvd between `pb_sorted` and `pd_sorted`, which should be an output of `sort_by_lost_photons(pb::MultipleCounts)`, for exactly `k` lost photons. Note that if you want the info regarding data considering up to `k` lost photons, you need to use [`tvd_less_than_k_lost_photons(k, pb_sorted, pd_sorted)`](@ref).
"""
function tvd_k_lost_photons(k, pb_sorted, pd_sorted)

    tvd(pb_sorted[k].proba,pd_sorted[k].proba)

end

"""
    tvd_less_than_k_lost_photons(k, pb_sorted, pd_sorted)

Gives the TVD obtained by considering 0,...,k photons lost. The tvd for each number of photons lost is summed using [`tvd_k_lost_photons(k, pb_sorted, pd_sorted)`](@ref).
"""
function tvd_less_than_k_lost_photons(k, pb_sorted, pd_sorted)

    sum(tvd(pb_sorted[j].proba,pd_sorted[j].proba) for j in 1:k)

end
