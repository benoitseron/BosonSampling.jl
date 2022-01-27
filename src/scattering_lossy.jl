n = 4
r = 2
m = n^2
mu = 0.8

complete_mode_occupation!(v) = append!(v, [0 for i in 1:m-r])

for arrangement in combinations(collect(1:n),r)
    println(complete_mode_occupation!(arrangement))
end


# scattering probability for a uniformly lossy
# interferometer with mu transmissivity



mu^(2*(r)) * (1-mu^2)^(n-r)

append!([1,3], [0 for i in 1:m-r])
