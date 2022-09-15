
"""
    Circuit(m::Int)

Creates an empty circuit with `m` input modes. The unitary representing the circuit
is accessed via the field `.U`.

    Fields:
        - m::Int
        - circuit_elements::Vector{Interferometer}
        - U::Union{Matrix{ComplexF64}, Nothing}
"""
mutable struct Circuit <: Interferometer

    m::Int
    circuit_elements::Vector{Interferometer}
    U::Union{Matrix, Nothing}

    function Circuit(m::Int)
        new(m, [], nothing)
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


"""
    add_element!(circuit::Circuit, interf::Interferometer; target_modes::Vector{Int})
    add_element!(;circuit::Circuit, interf::Interferometer, target_modes::Vector{Int})

Adds the circuit element `interf` that will be applied on `target_modes` to the `circuit`.
Will automatically update the unitary representing the circuit.
"""
function add_element!(;circuit::Circuit, interf::Interferometer, target_modes::Vector{Int})

    @argcheck interf.m == length(target_modes)

    push!(circuit.circuit_elements, interf)
    circuit.U == nothing ? circuit.U = Matrix{ComplexF64}(I, circuit.m, circuit.m) : nothing

    if interf.m == circuit.m
        circuit.U = interf.U * circuit.U
    else
        u = Matrix{ComplexF64}(I, circuit.m, circuit.m)
        for i in 1:size(interf.U)[1]
            for j in 1:size(interf.U)[2]
                u[target_modes[i], target_modes[j]] = interf.U[i,j]
            end
        end
        circuit.U = u * circuit.U
    end

end



add_element!(circuit::Circuit, interf::Interferometer; target_modes::Vector{Int}) = add_element!(circuit=circuit, interf=interf, target_modes=target_modes)

Base.show(io::IO, interf::Interferometer) = print(io, "Interferometer :\n\n", "Type : ", typeof(interf), "\n", "m : ", interf.m)

Base.show(io::IO, interf::UserDefinedInterferometer) = begin
     print(io, "Interferometer :\n\n", "Type : ", typeof(interf), "\n", "m : ", interf.m, "\nUnitary : \n")
     pretty_table(io, interf.U)
end
