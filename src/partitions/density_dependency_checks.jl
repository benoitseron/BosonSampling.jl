max_density = 1
min_density = 0.03
steps = 10
n_iter = 100

invert_densities = [max_density * (max_density/min_density)^((i-1)/(steps-1)) for i in 1:steps]

function power_law_with_n(n)

    m_array = Int.(floor.(n * invert_densities))
    partition_sizes = 2:2

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
    power_law_with_n(n)
end

# n = 7
# power law: y = 0.4333180401690799 * x^0.9796220494428659
# n = 9
# power law: y = 0.4140980068718625 * x^0.9579794008667261
# n = 11
# power law: y = 0.3953762282562418 * x^0.9261607896204195
# n = 13
# power law: y = 0.3883203879554887 * x^0.9041350377811694
# n = 15
# power law: y = 0.37997840718797066 * x^0.8905876525664023
