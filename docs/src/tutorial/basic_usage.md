# Basic usage

This tutorial will introduce you to the definition of the building blocks for your
boson sampling experiment. The general workflow for a simple simulation is to define
an [`Input`](@ref) that enters into a [`Interferometer`](@ref) and ask what is the
probability to get a defined `OutputMeasurement`.

## Input

`BosonSampling.jl` provides three distinct types of input depending on the
distinguishability of the particles we want to make interfere: [`Bosonic`](@ref),
[`PartDist`](@ref) and [`Distinguishable`](@ref). In order to define the input, we first need to provide a [`ModeOccupation`](@ref) that describes the repartition of the particles among the modes.

    julia> n = 3; # photon number

    julia> m = 6; # mode number

    julia> my_mode_occupation = ModeOccupation(random_occupancy(n,m))
     state = [1, 0, 0, 1, 1, 0]

In the example above, `my_mode_occupation` has been created thanks to [`random_occupancy`](@ref) that randomly places `n` particles among `m` modes. Here we have one particle in the first, fourth and fifth modes.
Let's build an input made off indistinguishable photons by using the type [`Bosonic`](@ref)

    julia> my_input = Input{Bosonic}(my_mode_occupation)
     Type:Input{Bosonic}
     r:state = [1, 1, 0, 0, 1, 0]
     n:3
     m:6
     G:GramMatrix{Bosonic}(3, ComplexF64[1.0 + 0.0im 1.0 + 0.0im 1.0 + 0.0im; 1.0 + 0.0im 1.0 + 0.0im 1.0 + 0.0im; 1.0 + 0.0im 1.0 + 0.0im 1.0 + 0.0im], nothing, nothing, OrthonormalBasis(nothing))
     distinguishability_param:nothing

where `my_input` holds the information defined above and an additional field, the [`GramMatrix`](@ref):

    help?> GramMatrix

     Fields:
      - n::Int: photons number
      - S::Matrix: Gram matrix
      - rank::Union{Int, Nothing}
      - distinguishability_param::Union{Real, Nothing}
      - generating_vectors::OrthonormalBasis     

which contains everything about the distinguishability of the particles within `my_input`. The matrix itself can be accessed via the field `S`:

    julia> my_input.G.S
     3×3 Matrix{ComplexF64}:
     1.0+0.0im  1.0+0.0im  1.0+0.0im
     1.0+0.0im  1.0+0.0im  1.0+0.0im
     1.0+0.0im  1.0+0.0im  1.0+0.0im

One can do the same for [`Distinguishable`](@ref) particles placed in the [`first_modes`](@ref)

    julia> my_mode_occupation = first_modes(n,m);

    julia> my_input = Input{Distinguishable}(my_mode_occupation);

    julia> my_input.G.S
     3×3 Matrix{ComplexF64}:
     1.0+0.0im  0.0+0.0im  0.0+0.0im
     0.0+0.0im  1.0+0.0im  0.0+0.0im
     0.0+0.0im  0.0+0.0im  1.0+0.0im

We can move now to the [`PartDist`](@ref) case with a model of partially distinguishable particles defined by a [`RandomGramMatrix`](@ref)

    julia> my_input = Input{RandomGramMatrix}(first_modes(n,m));

where `my_input.G.S` is a randomly generated Gram matrix.

Finally, one can resort to a [`OneParameterInterpolation`](@ref) model taking a linear distinguishability
parameter as an additional argument in the definition of `my_input`:

    julia> my_mode_occupation = ModeOccupation(random_occupancy(n,m));

    julia> my_distinguishability_param = 0.7;

    julia> my_input = Input{OneParameterInterpolation}(my_mode_occupation, my_distinguishability_param)
     Type:Input{OneParameterInterpolation}
     r:state = [1, 1, 1, 0, 0, 0]
     n:3
     m:6
     G:GramMatrix{OneParameterInterpolation}(3, [1.0 0.7 0.7; 0.7 1.0 0.7; 0.7 0.7 1.0], nothing, 0.7, OrthonormalBasis(nothing))
     distinguishability_param:0.7

Notice that the [`Bosonic`](@ref) Gram matrix is recovered for `my_distinguishability_param = 1`
while we find the [`Distinguishable`](@ref) case for `my_distinguishability_param = 0`.  

## Interferometer

The second building block of our boson sampler is the interferometer
we want to apply on `my_input`. A common practice to study boson sampling is to
pick up at random a Haar distributed unitary matrix that will represent the interferometer.
This can be done as follow:

    julia> my_random_interf = RandHaar(m);

    julia> my_random_interf.U
     6×6 Matrix{ComplexF64}:
     -0.398793-0.162392im     0.500654+0.239935im   -0.0822224-0.189055im     0.361289+0.048139im   -0.0639807+0.521608im   -0.0608676-0.226158im
     0.331232+0.312918im    -0.362433-0.0585051im    0.172619+0.157846im       0.6224-0.0656408im   -0.186285+0.215539im    -0.233294-0.274952im
     -0.085377+0.00861048im  -0.427763+0.1271im      -0.140581-0.541889im   -0.0407117+0.219472im     0.499523-0.0486383im  -0.0711764-0.416309im
     -0.180636+0.294002im    -0.376353+0.0943096im   -0.489498-0.206616im   0.00169099-0.221399im    -0.260862+0.305118im     0.313454+0.373737im
     0.619915-0.0776736im     0.29879-0.323099im    -0.398401-0.371214im    -0.132369+0.0411683im   -0.242868+0.0913286im  -0.0282651-0.179273im
     0.179322+0.240927im    0.0369079+0.110875im     0.100647-0.0206654im   -0.153966+0.577663im     0.154379+0.372127im    -0.366647+0.481086im

where we have accessed to the matrix thanks to the field `.U`.
We may also need to use a specific interferometer such as a [`quantum Fourier transform`](https://en.wikipedia.org/wiki/Quantum_Fourier_transform) or the [`Hadamard transform`](https://en.wikipedia.org/wiki/Hadamard_transform):

    julia> my_fourier_interf = Fourier(m);

    julia> is_unitary(my_fourier_interf.U)
     true

    julia> my_hadamard_tf = Hadamard(2^m);

    julia> is_unitary(my_hadamard_tf.U)
     true

where we have checked the unitarity thanks to `is_unitary`.

The implemented interferometers are listed in [`Interferometer`](@ref) but it is still
possible to define our own unitary by resorting to the type [`UserDefinedInterferometer`](@ref):

    julia> sigma_y = [0 -1im; 1im 0];

    julia> my_interf = UserDefinedInterferometer(sigma_y)
     Type : UserDefinedInterferometer
     m : 2
     Unitary :
     ┌────────┬────────┐
     │ Col. 1 │ Col. 2 │
     ├────────┼────────┤
     │  0+0im │  0-1im │
     │  0+1im │  0+0im │
     └────────┴────────┘

## OutputMeasurement

Similary to the definition of the [`Input`](@ref), it is also possible to define an output configuration from a [`ModeOccupation`](@ref)

    julia> n = 3;

    julia> m=n^2;

    julia> my_mode_occupation = first_modes(n,m);

    julia> my_input = Input{Bosonic}(my_mode_occupation)
     Type:Input{Bosonic}
     r:state = [1, 1, 1, 0, 0, 0, 0, 0, 0]
     n:3
     m:9
     G:GramMatrix{Bosonic}(3, ComplexF64[1.0 + 0.0im 1.0 + 0.0im 1.0 + 0.0im; 1.0 + 0.0im 1.0 + 0.0im 1.0 + 0.0im; 1.0 + 0.0im 1.0 + 0.0im 1.0 + 0.0im], nothing, nothing, OrthonormalBasis(nothing))
     distinguishability_param:nothing

    julia> out = FockDetection(my_mode_occupation)
     FockDetection(state = [1, 1, 1, 0, 0, 0, 0, 0, 0])

using [`FockDetection`](@ref). Additionally, we can define an [`Event`](@ref) that stores our input-interferometer-output content

    julia> my_interf = Fourier(my_input.m)
     Type : Fourier
     m : 9

    julia> ev = Event(my_input, out, my_interf)
     Event{Bosonic, FockDetection}(Input:

     Type:Input{Bosonic}
     r:state = [1, 1, 1, 0, 0, 0, 0, 0, 0]
     n:3
     m:9
     G:GramMatrix{Bosonic}(3, ComplexF64[1.0 + 0.0im 1.0 + 0.0im 1.0 + 0.0im; 1.0 + 0.0im 1.0 + 0.0im 1.0 + 0.0im; 1.0 + 0.0im 1.0 + 0.0im 1.0 + 0.0im], nothing, nothing, OrthonormalBasis(nothing))
     distinguishability_param:nothing, FockDetection(state = [1, 1, 1, 0, 0, 0, 0, 0, 0]), EventProbability(nothing, nothing, nothing), Interferometer :

     Type : Fourier
     m : 9)

and then one can compute the probability that this event occurs

    julia> compute_probability!(ev)
     0.015964548319225575

Those steps can be repeated for different types of input. Let's say we want to
compute the probability that partially distinguishable photons pupulating the `n=3`
first modes of `m=9` modes end up in the `n` last output modes when interfering through
a random interferometer:

    julia> my_input = Input{RandomGramMatrix}(first_modes(n,m)); # input from a randomly generated Gram matrix

    julia> out = FockDetection(ModeOccupation([0,0,0,0,0,0,1,1,1]));

    julia> my_interf = RandHaar(m);

    julia> ev = Event(my_input, out, my_interf);

    julia> compute_probability!(ev)
     0.014458823860031098
