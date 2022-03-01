# this file checks the conventions

# see conventions.md

using BosonSampling
using Test

U = matrix_test(2)

interferometer_convention = "tichy"
@test scattering_matrix(U, [1,1],[2,0], interferometer_convention) == [[1 1]; [3 3]]

@test scattering_matrix(U, [1,1],[2,0], "tichy") == scattering_matrix(U, [1,1],[2,0], "shchesnovich")
