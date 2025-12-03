using BosonSampling
using BenchmarkTools

println("="^70)
println("Benchmarking sample_multiple vs repeated sample!()")
println("="^70)

# Large system parameters
n = 10  # photons
m = 100  # modes
n_samples = 100

println("\nParameters:")
println("  n = $n photons")
println("  m = $m modes")
println("  n_samples = $n_samples")

# Create input with partial distinguishability
r = 5  # rank of Gram matrix
S = rand_gram_matrix_from_orthonormal_basis(n, r)
input = Input{UserDefinedGramMatrix}(first_modes(n, m), S)
interf = RandHaar(m)

println("\n" * "="^70)
println("Method 1: Repeated sample!() calls (slow)")
println("="^70)

time_repeated = @belapsed begin
    for i in 1:$n_samples
        ev = Event($input, FockSample(), $interf)
        sample!(ev)
    end
end

println("\nTotal time for $n_samples samples: $(round(time_repeated, digits=3)) s")
println("Time per sample: $(round(time_repeated*1000/n_samples, digits=2)) ms")

println("\n" * "="^70)
println("Method 2: sample_multiple() (fast - pre-computes transformations)")
println("="^70)

time_multiple = @belapsed sample_multiple($input, $interf, $n_samples)

println("\nTotal time for $n_samples samples: $(round(time_multiple, digits=3)) s")
println("Time per sample: $(round(time_multiple*1000/n_samples, digits=2)) ms")

println("\n" * "="^70)
println("Performance Comparison")
println("="^70)

speedup = time_repeated / time_multiple
time_saved = time_repeated - time_multiple

println("\n🚀 Speedup: $(round(speedup, digits=2))x faster!")
println("⏱️  Time saved: $(round(time_saved, digits=3)) s ($(round(time_saved/time_repeated*100, digits=1))% reduction)")
println("\nFor $n_samples samples:")
println("  Old method: $(round(time_repeated, digits=3)) s")
println("  New method: $(round(time_multiple, digits=3)) s")

# Breakdown of what's pre-computed
println("\n" * "="^70)
println("What sample_multiple() pre-computes (done once):")
println("="^70)
println("  ✓ Gram matrix decomposition")
println("  ✓ $n Householder matrices ($(n)×$(r)×$(r) = $(n*r*r) elements)")
println("  ✓ Shuffle permutation matrix ($(n*m)×$(n*m) = $(n*m*n*m) elements)")
println("  ✓ Block diagonal interferometers")
println("  ✓ Full interferometer construction (matrix multiplications)")
println("\nRepeated for each sample:")
println("  • Clifford sampling in expanded space")
println("  • Binning from $(n*m) modes to $m physical modes")

println("\n" * "="^70)
println("✓ Benchmark complete!")
println("="^70)
