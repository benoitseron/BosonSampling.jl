using BosonSampling
using BenchmarkTools

println("="^70)
println("Benchmarking Parallel vs Sequential sample_multiple()")
println("="^70)

# Check available threads
n_threads = Threads.nthreads()
println("\nAvailable threads: $n_threads")
if n_threads == 1
    println("⚠️  Warning: Running with only 1 thread!")
    println("   Set JULIA_NUM_THREADS environment variable for parallel speedup")
    println("   Example: JULIA_NUM_THREADS=8 julia benchmark_threaded.jl")
end

# Parameters
n = 10  # photons
m = 100  # modes
n_samples = 1000

println("\nParameters:")
println("  n = $n photons")
println("  m = $m modes")
println("  n_samples = $n_samples")

# Create input
r = 5
S = rand_gram_matrix_from_orthonormal_basis(n, r)
input = Input{UserDefinedGramMatrix}(first_modes(n, m), S)
interf = RandHaar(m)

println("\n" * "="^70)
println("Method 1: Sequential sample_multiple()")
println("="^70)

time_sequential = @belapsed sample_multiple($input, $interf, $n_samples, threaded=false)

println("\nTotal time: $(round(time_sequential, digits=3)) s")
println("Per sample: $(round(time_sequential*1000/n_samples, digits=3)) ms")

println("\n" * "="^70)
println("Method 2: Parallel sample_multiple() with $n_threads threads")
println("="^70)

time_parallel = @belapsed sample_multiple($input, $interf, $n_samples, threaded=true)

println("\nTotal time: $(round(time_parallel, digits=3)) s")
println("Per sample: $(round(time_parallel*1000/n_samples, digits=3)) ms")

println("\n" * "="^70)
println("Performance Comparison")
println("="^70)

if n_threads > 1
    speedup = time_sequential / time_parallel
    efficiency = speedup / n_threads * 100

    println("\n🚀 Parallel speedup: $(round(speedup, digits=2))x faster!")
    println("⚡ Parallel efficiency: $(round(efficiency, digits=1))% (ideal: 100%)")
    println("\nTime comparison:")
    println("  Sequential: $(round(time_sequential, digits=3)) s")
    println("  Parallel:   $(round(time_parallel, digits=3)) s")
    println("  Time saved: $(round(time_sequential - time_parallel, digits=3)) s")

    # Combined speedup vs original repeated sample!()
    # Estimate: original ~0.16s per sample based on previous benchmark
    time_original_estimate = 0.16 * n_samples
    combined_speedup = time_original_estimate / time_parallel
    println("\n" * "="^70)
    println("Combined Speedup (vs original repeated sample!)")
    println("="^70)
    println("  Estimated original time: $(round(time_original_estimate, digits=1)) s")
    println("  New parallel time: $(round(time_parallel, digits=3)) s")
    println("  🎯 Total speedup: $(round(combined_speedup, digits=1))x faster!")
else
    println("\n⚠️  No parallel speedup with 1 thread")
    println("   Run Julia with multiple threads to see parallel benefits:")
    println("   JULIA_NUM_THREADS=8 julia benchmark_threaded.jl")
end

println("\n" * "="^70)
println("✓ Benchmark complete!")
println("="^70)
