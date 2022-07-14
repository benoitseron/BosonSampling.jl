"""
    virtual_interferometer_uniform_loss(real_interf::Matrix, η)
    virtual_interferometer_uniform_loss(real_interf::Interferometer, η)

Simulates a simple, uniformly lossy interferometer: take a 2m*2m interferometer and introduce beam splitters in front with transmittance `η`. Returns the corresponding virtual interferometer.
"""
function virtual_interferometer_uniform_loss(real_interf::Matrix, η)

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

    m = size(U,1)
    U_virtual = Matrix{eltype(U)}(I, (2m, 2m))
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

    environment = Subset(last_modes(n,2m))

    new_subsets =  Vector{Subset}()

    for subset in part.subsets
        push!(new_subsets, to_lossy(subset))
    end

    push!(new_subsets, environment)

    Partition(new_subsets)

end

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
    m_virtual::Int
    η::Real #transmissivity of the upfront beamsplitters
    U_physical::Matrix{Complex} # physical matrix
    U_virtual::Matrix{Complex} # virtual 2m*2m interferometer

    function UniformLossInterferometer(η::Real, U_physical::Matrix)
        if isa_transmissitivity(η)
            U_virtual = virtual_interferometer_uniform_loss(U_physical, η)
            new(size(U_physical,1), size(U_virtual,1), η, U_physical, U_virtual)
        else
            error("incorrect η")
        end
    end

    UniformLossInterferometer(η::Real, U_physical::Interferometer) = UniformLossInterferometer(η, U_physical.U)

    UniformLossInterferometer(m::Int, η::Real) = UniformLossInterferometer(η, RandHaar(m))

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
    m_virtual::Int
    η::Vector{Real} #transmissivity of the upfront beamsplitters
    U_physical::Matrix{Complex} # physical matrix
    V::Matrix{Complex}
    W::Matrix{Complex}
    U_virtual::Matrix{Complex} # virtual 2m*2m interferometer

    function GeneralLossInterferometer(η::Vector{Real}, V::Matrix, W::Matrix)
        if isa_transmissitivity(η)
            U_virtual = virtual_interferometer_general_loss(V,W, η)
            new(size(U_physical,1), size(U_virtual,1), η, U_physical,V,W, U_virtual)
        else
            error("invalid η")
        end
    end

end
