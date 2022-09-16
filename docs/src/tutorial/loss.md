# Loss

Loss can be incorporated through [`BeamSplitter`](@ref)'s sending photons with some probability to extra environment modes. If a physical [`Interferometer`](@ref) has `m` modes, we create extra `m` modes representing lost photons. In reality, these would not be accessible, but we may still keep this information if necessary. This allows to post-select events upon a certain loss pattern, such as finding `l` (lost) photons in the environment modes.

## Conversions

In general, the function [`to_lossy`](@ref) converts physical `m`-mode objects into their `2m`-modes counterpart fitting the above model. For instance

    julia> n=3
    m=4

    first_modes(n,m)
    state = [1, 1, 1, 0]

    to_lossy(first_modes(n,m))
    state = [1, 1, 1, 0, 0, 0, 0, 0]

    # creating a Subset:
    Subset(first_modes(n,m)).m
    4

    # expanding it doesn't change the Subset
    to_lossy(Subset(first_modes(n,m)))
    subset = [1, 2, 3]

    # but it is now of the correct size
    to_lossy(Subset(first_modes(n,m))).m
    8

## Conventions

Each circuit element, such as [`BeamSplitter`](@ref) and [`PhaseShift`](@ref) can bear a certain amount of loss. We write it `η_loss`. It is the transmission amplitude of the beam splitter representing the loss process. Therefore the probability that a photon is not lost is `η_loss^2`.

## Lossy interferometers

The inclusion of loss creates bigger [`Interferometer`](@ref)'s, but half of their modes are not physical. For this reason, we use the subtype [`LossyInterferometer`](@ref).

The fields are named in such a way that all computations can be done without changes, as if we now used a `2m*2m` lossless interferometer. The physical quantities are labelled accordingly such as `m_real` and `U_physical`.

## Models implemented

Let us now discuss the various lossy elements available.
* [`UniformLossInterferometer`](@ref) : This simplest model is one where photons have an identical chance of being lost.
* [`GeneralLossInterferometer`](@ref) This is a generic model as described in ...
* Lossy circuit elements : When constructing a [`Circuit`](@ref) from elements, each element has its own loss characteristics. We also introduce lines, representing for instance optical fibers that have no interaction but can still be lossy.

## Circuits

When using `circuit_elements` to construct a lossy interferometer, the loss channel associated to mode `i` will always be mode `m+i`. Therefore, doing 
