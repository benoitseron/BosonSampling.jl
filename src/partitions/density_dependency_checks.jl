max_density = 1
min_density = 0.03
steps = 30
n_iter = 100

invert_densities = [max_density * (max_density/min_density)^((i-1)/(steps-1)) for i in 1:steps]

function power_law_with_n(n,k)

    partition_sizes = k:k
    m_array = Int.(floor.(n * invert_densities))


    tvd_array = zeros((length(partition_sizes), length(m_array)))
    var_array = copy(tvd_array)

    for (k,n_subsets) in enumerate(partition_sizes)

        #@show n_subsets

        for i in 1:length(m_array)

            this_tvd = tvd_equilibrated_partition_real_average(m_array[i], n_subsets, n, niter = n_iter)

            tvd_array[k,i] = (n_subsets <= m_array[i] ? this_tvd[1] : missing)
            var_array[k,i] = (n_subsets <= m_array[i] ? this_tvd[2] : missing)
        end

    end

    x_data = reverse(1 ./ invert_densities)
    y_data = reverse(tvd_array[1,:])

    get_power_law_log_log(x_data,y_data)

end

for n in 5:2:13
    @show n
    power_law_with_n(n,2)
end

# n = 7
# power law: y = 0.44038585499823646 * x^0.982801094275387
# n = 9
# power law: y = 0.4232947463279576 * x^0.9788828718055166
# n = 11
# power law: y = 0.4123148441313412 * x^0.9564544056604489
# n = 13
# power law: y = 0.4052999461220922 * x^0.9403720479501786

# getting the constants

for k in 2:3

    println("c($k) = $(power_law_with_n(5,k)[3])")

end
