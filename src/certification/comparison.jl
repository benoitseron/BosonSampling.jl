### comparison with subset certifcation ###

plotly()

n_array = collect(2:20)

m_array = 5 .* n_array

const_density = [choose_best_average_subset(n = n_array[i],m = m_array[i])[2] for i in 1:length(n_array)]
const_density_full_bunching = [choose_subset_size(n = n_array[i],m = m_array[i])[3]  for i in 1:length(n_array)]

scatter(n_array, const_density, label = "partitions");
scatter!(const_density_full_bunching, label = "full bunching");
title!("tvd const density m = 5n")
savefig("src/certification/images/tvd_const_density")

m_array = n_array .^2

no_collision = [choose_best_average_subset(n = n_array[i],m = m_array[i])[2] for i in 1:length(n_array)]
no_collision_full_bunching = [choose_subset_size(n = n_array[i],m = m_array[i])[3] for i in 1:length(n_array)]

scatter(n_array, no_collision, label = "partitions");
scatter!(no_collision_full_bunching, label = "full bunching");
title!("tvd no collision density m = n^2")
savefig("src/certification/images/tvd_no_collision")
