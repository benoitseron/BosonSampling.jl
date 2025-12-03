"""
Generate boson sampling data and save to JLD2.

Usage:
    julia --project=. examples/generate_samples.jl

Output:
    Creates a .jld2 file with all parameters and samples.
"""

using BosonSampling
using JLD2
using Dates
using ProgressMeter

# =============================================================================
# Parameters
# =============================================================================

n = 10             # number of photons
m = 100             # number of modes
n_samples = 1000    # number of samples to generate

# Gram matrix (partial distinguishability)
r = 10               # rank of Gram matrix
S = rand_gram_matrix_from_orthonormal_basis(n, r)

# Alternatively, use your own Gram matrix:
# S = your_gram_matrix

# =============================================================================
# Setup
# =============================================================================

input = Input{UserDefinedGramMatrix}(first_modes(n, m), S)
interf = RandHaar(m)
U = interf.U  # store the unitary matrix

# =============================================================================
# Generate samples
# =============================================================================

println("Generating $n_samples samples (n=$n, m=$m, rank=$r)")

samples = sample_multiple(input, interf, n_samples, show_progress=true)

# =============================================================================
# Save to file
# =============================================================================

outdir = "data/generated_samples"
mkpath(outdir)

filename = joinpath(outdir, "boson_samples_n$(n)_m$(m)_$(Dates.format(now(), "yyyymmdd_HHMMSS")).jld2")

jldsave(filename;
    # Samples
    samples,

    # Physical parameters
    n,
    m,
    n_samples,

    # Gram matrix (distinguishability)
    S,
    gram_rank = r,

    # Interferometer
    U,

    # Metadata
    generated_at = now(),
    package_version = "BosonSampling.jl"
)

println("Saved to: $filename")

# =============================================================================
# How to load the data
# =============================================================================
#=
using JLD2

data = load("data/generated_samples/your_file.jld2")

samples = data["samples"]       # Vector of ModeOccupation
S = data["S"]                   # Gram matrix
U = data["U"]                   # Unitary matrix
n = data["n"]                   # photon number
m = data["m"]                   # mode number
=#
