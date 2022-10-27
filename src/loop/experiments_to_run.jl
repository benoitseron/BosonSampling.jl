###### first experiment ######

n = 8
sparsity = 2
m = sparsity * n
n_subsets = 2
n_subsets > 3 && m > 12 ? (@warn "may be slow") : nothing

x = 0.9

equilibrated_input(sparsity, m) = ModeOccupation([((i-1) % sparsity) == 0 ? 1 : 0 for i in 1:m])

equilibrated_input(sparsity, m)

(-tvd_reflectivities(η_thermalization(m)))

part = equilibrated_partition(m, n_subsets)

###### second experiment ######

n = 3
sparsity = 3
m = sparsity * n
n_subsets = 2
n_subsets > 3 && m > 12 ? (@warn "may be slow") : nothing

x = 0.9

equilibrated_input(sparsity, m) = ModeOccupation([((i-1) % sparsity) == 0 ? 1 : 0 for i in 1:m])

equilibrated_input(sparsity, m)

sol.minimizer = [3.440486091233136e-5, 2.5123712802059534e-5, 0.7071242076767412, 2.8385514629156598e-5, 1.026309739414248e-5, 0.9920363137406634, 0.9885666673304567, 0.7392252054990487]

tvd(pdf_bos, pdf_dist) = 0.8099999950860921

###### third ######

n = 2
sparsity = 3
m = sparsity * n
n_subsets = 3
n_subsets > 3 && m > 12 ? (@warn "may be slow") : nothing

x = 0.9

equilibrated_input(sparsity, m) = ModeOccupation([((i-1) % sparsity) == 0 ? 1 : 0 for i in 1:m])

equilibrated_input(sparsity, m)

sol.minimizer = [0.00010685158506912827, 4.484886054731138e-5, 0.7071202195365436, 3.7341888184701493e-5, 0.9553391960360447]

tvd(pdf_bos, pdf_dist) = 0.8099999856935723

###### pseudo_photon_number_resolution ######



n = 10
steps_pnr = 2 # number of bins for pseudo pnr
steps_pnr > 2 && n > 6 ? (@warn "may be slow") : nothing
m = n + steps_pnr

x = 0.8

d = Uniform(0,2pi)
ϕ = nothing # rand(d,m)
η_loss_lines = 0.9 * ones(m)
η_loss_bs = 1. * ones(m-1)

η = vcat(η_thermalization(n), η_pnr(steps_pnr))

η = [0.75, 0.5555555555555556, 0.4375, 0.3599999999999999, 0.30555555555555547, 0.26530612244897966, 0.234375, 0.2098765432098766, 0.18999999999999995, 0.3333333333333333, 0.6666666666666666]

# in this case the input is r → state = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0]

# you can run the thermalization file to see the distribution, and look at the output results of mc_b, mc_d for which patterns are seen with what probability
