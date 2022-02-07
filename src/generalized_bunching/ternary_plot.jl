include("counter_example_circuit.jl")

H_circuit = H

n = 7
S_bosonic = ones(ComplexF64, n,n)
S_dist = Matrix{ComplexF64}(I,n,n)

p_bosonic = abs(permanent_ryser(H_circuit .* S_bosonic))
p_counter = abs(permanent_ryser(H_circuit .* S_circuit))

function ternary_point(x,y; H_circuit, S_circuit)

    """returns the relative bunching probability P(S(x,y))/P(S_bosonic)
    for
        S(x,y) = x * S_bosonic + y * S_dist + (1-x-y) * S_circuit
    in the n=7 counter example circuit"""

    S(x,y) = x * S_bosonic + y * S_dist + (1-x-y) * S_circuit

    abs(permanent_ryser(H_circuit .* (S(x,y)))) / p_bosonic

end

ternary_point(1,0; H_circuit = H_circuit, S_circuit = S_circuit)

n_x = 301
n_y = 301
x_array = range(0,1, length = n_x)
y_array = range(0,1, length = n_y)

ternary_data = Matrix(undef,n_x*n_y, 3)

i = 1

for x in x_array
    for y in y_array

        if x+y<=1-10eps() #so that it is a ternary plot

            ternary_data[i,1] = x
            ternary_data[i,2] = y
            ternary_data[i,3] = ternary_point(x,y; H_circuit = H_circuit, S_circuit = S_circuit)

            i+=1
        end
    end
end

# now we need to remove undefined points

last_defined_index = 1

for i in 1:n_x*n_y
    try ternary_data[i,1]

    catch
        last_defined_index = i-1
        break
    end
end

ternary_data = ternary_data[1:last_defined_index, :]

CSV.write("counter_examples/ternary_data.csv",  DataFrame(ternary_data,:auto), header=false)

### sub ternary ###

ternary_point(1,0; H_circuit = H_circuit, S_circuit = S_circuit)

n_x = 301
n_y = 301
x_array = range(0,1, length = n_x)
y_array = range(0,0.2, length = n_y)

ternary_data = Matrix(undef,n_x*n_y, 3)

i = 1

for x in x_array
    for y in y_array

        if x+y<=1-10eps() #so that it is a ternary plot

            ternary_data[i,1] = x
            ternary_data[i,2] = y
            ternary_data[i,3] = ternary_point(x,y; H_circuit = H_circuit, S_circuit = S_circuit)

            i+=1
        end
    end
end

# now we need to remove undefined points

last_defined_index = 1

for i in 1:n_x*n_y
    try ternary_data[i,1]

    catch
        last_defined_index = i-1
        break
    end
end

ternary_data = ternary_data[1:last_defined_index, :]

CSV.write("counter_examples/sub_ternary_data.csv",  DataFrame(ternary_data,:auto), header=false)
