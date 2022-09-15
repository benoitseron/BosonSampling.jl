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

"""
    LossyInterferometer <: Interferometer

Interferometers with inclusion of loss: the real `Interferometer` has dimension `m_real * m_real` while we model it by a `2m_real * 2m_real` one where the last `m_real` modes are environment modes containing the lost photons.
"""
abstract type LossyInterferometer <: Interferometer end

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
struct UniformLossInterferometer <: LossyInterferometer
    m_real::Int
    m::Int
    η::Real #transmissivity of the upfront beamsplitters
    U_physical::Matrix{Complex} # physical matrix
    U::Matrix{Complex} # virtual 2m*2m interferometer

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
struct GeneralLossInterferometer <: LossyInterferometer
    m_real::Int
    m::Int
    η::Vector{Real} #transmissivity of the upfront beamsplitters
    U_physical::Matrix{Complex} # physical matrix
    V::Matrix{Complex}
    W::Matrix{Complex}
    U::Matrix{Complex} # virtual 2m*2m interferometer

    function GeneralLossInterferometer(η::Vector{Real}, V::Matrix, W::Matrix)
        if isa_transmissitivity(η)
            U_virtual = virtual_interferometer_general_loss(V,W, η)
            new(size(U_physical,1), size(U_virtual,1), η, U_physical,V,W, U_virtual)
        else
            error("invalid η")
        end
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

"""
    LossyBeamSplitter(transmission_amplitude, η_loss)

Creates a beam-splitter with tunable transmissivity and loss. Uniform model of loss: each input line i1,i2 has a beam splitter in front with transmission amplitude of `transmission_amplitude_loss` into an environment mode.
"""
struct LossyBeamSplitter <: Interferometer
    transmission_amplitude::Real
    η_loss::Real
    U::Matrix
    m::Int
    function LossyBeamSplitter(transmission_amplitude::Real, η_loss::Real)
        @argcheck between_one_and_zero(transmission_amplitude)
        @argcheck between_one_and_zero(η_loss)

        new(transmission_amplitude, η_loss, virtual_interferometer_uniform_loss(beam_splitter(transmission_amplitude),η_loss), 4)
    end
end


"""
    LossyLine(η_loss)

Optical line with some loss: represented by a `BeamSplitter` with a transmission amplitude of `transmission_amplitude_loss` into an environment mode.
"""
struct LossyLine <: Interferometer

    η_loss::Real
    U::Matrix
    m::Int

    function LossyLine(η_loss::Real)

        @argcheck between_one_and_zero(η_loss)

        new(η_loss, virtual_interferometer_uniform_loss(ones((1,1)), η_loss), 2)
    end
end

"""
    LossyCircuit(m_real::Int)

Lossy `Interferometer` constructed from `circuit_elements`.
"""
mutable struct LossyCircuit <: LossyInterferometer

    m_real::Int
    m::Int
    circuit_elements::Vector{Interferometer}
    U_physical::Union{Matrix{Complex}, Nothing} # physical matrix
    U::Union{Matrix{Complex}, Nothing} # virtual 2m*2m interferometer

    function LossyCircuit(m_real::Int)
        new(m_real, 2*m_real, [], nothing, nothing)
    end

end
