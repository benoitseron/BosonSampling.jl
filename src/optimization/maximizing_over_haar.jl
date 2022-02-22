# brute force tests to see which matrices U maximize tvd(P^D - P^B)


function brute_force_maximize_over_haar_partition_pdf(; n_modes = 4, n_trials = 10,n_photons = 2, partition_size = 2, dist = tvd, saveimg = true, start_with_fourier_special_partition = false)

    if start_with_fourier_special_partition
        # apply the inverse transformation to the fourier matrix so that if the output for U = fourier_matrix
        # would be (1 0 1 0 1 0) it is now (111000)
        # note that this special partition only really make sense for partition size = n_modes / 2
        # thus we will skip all other partition sizes

        if 2*partition_size != n_modes
            throw(ArgumentError("special partition only make sense when the partition is half the modes"))
        elseif n_modes%2 != 0
            throw(ArgumentError("n_modes must be even to make sense"))
        end

        U_best = inv(permutation_matrix_special_partition_fourier(n_modes)) * fourier_matrix(n_modes)
    else
        U_best = fourier_matrix(n_modes)
    end

    is_the_fourier_matrix = true
    dist_best = distance_partition_pdfs(U = U_best, n_photons = n_photons, partition_size = partition_size, dist = dist)

    for trial in 1:n_trials

        U_this_iter = rand_haar(n_modes)
        dist_this_iter = distance_partition_pdfs(U = U_this_iter, n_photons = n_photons, partition_size = partition_size, dist = dist)

        if dist_this_iter > dist_best

            U_best = U_this_iter
            dist_best = dist_this_iter
            is_the_fourier_matrix = false

        end

    end

    if saveimg
        # plotting the best result

        cd(starting_directory)
        make_directory("images", delete_contents_if_existing = false)
    	make_directory("most sensitive interferometer", delete_contents_if_existing = false)
        make_directory("$dist", delete_contents_if_existing = false)
        make_directory("fourier_special_partition : $start_with_fourier_special_partition", delete_contents_if_existing = false)

        input_state = zeros(Int, n_modes)
        input_state[1:n_photons] = ones(Int,n_photons)
        partition_occupancy_vector = zeros(Int, n_modes)
        partition_occupancy_vector[1:partition_size] = ones(Int, partition_size)

        part = occupancy_vector_to_partition(partition_occupancy_vector)

        pdf_bosonic = proba_partition(U_best, partition_occupancy_vector, input_state = input_state)
        pdf_dist = partition_probability_distribution_distinguishable(part, U_best, number_photons = n_photons)

        x = 0:n_photons
        y = [pdf_bosonic, pdf_dist]

        plt = Plots.scatter(x,y, label = ["bosonic" "dist"], ylims = (0,1.05))

        matrix_plot_real = heatmap(real.(U_best))
        matrix_plot_img = heatmap(imag.(U_best))
        #display(plt)
        savefig(plt, "most sensitive matrix to dist (distributions), n_modes = $n_modes, n_photons = $n_photons, partition_size = $partition_size, n_trials = $n_trials, is_the_fourier_matrix = $is_the_fourier_matrix")
        savefig(matrix_plot_real, "most sensitive matrix to dist (real), n_modes = $n_modes, n_photons = $n_photons, partition_size = $partition_size, n_trials = $n_trials, is_the_fourier_matrix = $is_the_fourier_matrix")
        savefig(matrix_plot_img, "most sensitive matrix to dist (img), n_modes = $n_modes, n_photons = $n_photons, partition_size = $partition_size, n_trials = $n_trials, is_the_fourier_matrix = $is_the_fourier_matrix")

        save_matrix(U_best, "most sensitive matrix to dist (matrix), n_modes = $n_modes, n_photons = $n_photons, partition_size = $partition_size, n_trials = $n_trials, is_the_fourier_matrix = $is_the_fourier_matrix")

        cd(starting_directory)

    end

    (dist_best, U_best)

end

for n_modes = 2:16
    for n_photons = 1:n_modes
        for partition_size = 1:n_modes
            try
                brute_force_maximize_over_haar_partition_pdf(;n_modes = n_modes, n_trials = 10000, n_photons = n_photons, partition_size = partition_size, dist = sqr, start_with_fourier_special_partition = true)
            catch
            end
        end
    end
end

brute_force_maximize_over_haar_partition_pdf(n_modes = 100, n_trials = 10000, n_photons = 10, partition_size = 50, dist = sqr, start_with_fourier_special_partition = true))
