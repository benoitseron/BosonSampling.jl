# Circuit

`BosonSampling.jl` allows you to build circuits made of the available [`Interferometer`](@ref)s
and [`UserDefinedInterferometer`](@ref)s.

The first step is to define an empty circuit that will acts on `m=6` modes

    julia> my_circuit = Circuit(6);
    Interferometer :

    Type : Circuit
    m : 6

Next, we add one by one the components of the circuit, precising on which modes the
newly added component acts, by using [`add_element!`](@ref)

    julia> add_element!(circuit=my_circuit, interf=Fourier(4), target_modes=[1,2,3,4])
    6×6 Matrix{ComplexF64}:
    0.5+0.0im           0.5+0.0im           0.5+0.0im                   0.5+0.0im          0.0+0.0im  0.0+0.0im
    0.5+0.0im   3.06162e-17+0.5im          -0.5+6.12323e-17im  -9.18485e-17-0.5im          0.0+0.0im  0.0+0.0im
    0.5+0.0im          -0.5+6.12323e-17im   0.5-1.22465e-16im          -0.5+1.83697e-16im  0.0+0.0im  0.0+0.0im
    0.5+0.0im  -9.18485e-17-0.5im          -0.5+1.83697e-16im   2.75546e-16+0.5im          0.0+0.0im  0.0+0.0im
    0.0+0.0im           0.0+0.0im           0.0+0.0im                   0.0+0.0im          1.0+0.0im  0.0+0.0im
    0.0+0.0im           0.0+0.0im           0.0+0.0im                   0.0+0.0im          0.0+0.0im  1.0+0.0im

Here, we just have added a [`Fourier`](@ref) interferometer to our circuit that takes
the modes `[1,2,3,4]` at the input. The output matrix is the unitary representing
our circuit and will be updated at each call of [`add_element!`](@ref).

Let's add some more elements to our circuit:

    julia> add_element!(circuit=my_circuit, interf=BeamSplitter(1/sqrt(2)), target_modes=[5,6]);

    julia> add_element!(circuit=my_circuit, interf=RandHaar(6), target_modes=[1,2,3,4,5,6]);

The unitary representing our circuit can be accessed via the field `.U`,
as for any [`Interferometer`](@ref)

    julia> my_circuit.U
    6×6 Matrix{ComplexF64}:
    -0.483121-0.246661im    -0.232022-0.397813im   -0.027111-0.1335im     0.296595-0.471387im     -0.25528-0.282524im   0.0866359-0.111526im
    -0.372488+0.0893226im   -0.184263+0.0697938im    0.51829+0.200831im   0.315289+0.238577im   -0.0988814+0.298748im   0.0206645+0.499711im
    -0.504719-0.322371im   -0.0289979+0.458437im   -0.312735-0.156324im  -0.147983+0.354067im    0.0997703-0.0821812im   0.378552-0.0285867im
    -0.26704-0.191813im     0.174817-0.217575im     0.28131+0.345502im  -0.337596-0.230349im     0.505173+0.318283im     0.13136-0.273313im
    -0.227676-0.170793im     0.538223+0.116277im    0.180029-0.501201im  -0.116825-0.0909765im   -0.214777+0.228563im   -0.460068+0.0147747im
    0.0875484-0.0598966im   -0.258087-0.300614im   0.0235678-0.259943im   0.131754+0.42417im     -0.191588+0.497691im   0.0959991-0.522251im

Finally, the components of the circuit can also be retrieved via the field `.circuit_elements`

    julia> my_circuit.circuit_elements
    3-element Vector{Interferometer}:
    Interferometer :

    Type : Fourier
    m : 4
    Interferometer :

    Type : BeamSplitter
    m : 2
    Interferometer :

    Type : RandHaar
    m : 6

that are stored in a `Vector{Interferometer}`.     
