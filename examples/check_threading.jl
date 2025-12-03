"""
Check threading configuration and provide recommendations
"""

using BosonSampling

println("="^70)
println("Julia Threading Configuration")
println("="^70)

cpu_threads = Sys.CPU_THREADS
julia_threads = Threads.nthreads()

println("\nSystem Information:")
println("  CPU threads available: $cpu_threads")
println("  Julia currently using: $julia_threads threads")

if julia_threads == 1
    println("\n⚠️  WARNING: Julia is running with only 1 thread!")
    println("\nTo enable multi-threading:")
    println("  1. VSCode: Add to settings.json:")
    println("     \"julia.additionalArgs\": [\"-t\", \"auto\"]")
    println("  2. Command line: julia -t auto your_script.jl")
    println("  3. Use all $cpu_threads cores: julia -t $cpu_threads your_script.jl")
elseif julia_threads < cpu_threads
    println("\n⚡ Julia is using $julia_threads of $cpu_threads available threads")
    println("   This is fine (leaves $(cpu_threads - julia_threads) for OS/other tasks)")
    println("\nTo use all $cpu_threads threads:")
    println("   julia -t $cpu_threads your_script.jl")
else
    println("\n✓ Julia is using all available threads!")
end

# Test parallel sampling
println("\n" * "="^70)
println("Testing Parallel Sampling")
println("="^70)

n, m = 3, 5
input = Input{UserDefinedGramMatrix}(first_modes(n, m), rand_gram_matrix_from_orthonormal_basis(n, 2))
interf = RandHaar(m)

println("\nGenerating 100 samples with threaded=true...")
samples = sample_multiple(input, interf, 100, threaded=true)

println("✓ Successfully generated $(length(samples)) samples")
println("✓ All samples conserve photon number: $(all(sum(s) == n for s in samples))")

println("\n" * "="^70)
println("Cluster Usage")
println("="^70)
println("\nFor cluster computing:")
println("  Single node (current code): julia -t \$(nproc) your_script.jl")
println("  Multi-node: Use Distributed.jl with @distributed or pmap")
println("\nExample SLURM job script:")
println("  #!/bin/bash")
println("  #SBATCH --ntasks=1")
println("  #SBATCH --cpus-per-task=$cpu_threads")
println("  julia -t \$SLURM_CPUS_PER_TASK your_script.jl")

println("\n" * "="^70)
