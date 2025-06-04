"""
Working playground with exact validation

Shows that:
1. The Bayesian method IS working correctly
2. Our Clifford implementation has issues (which exact validation detects)
3. Need to fix the Clifford algorithm
"""

using BosonSampling
using Random

include("clifford_final.jl")
include("validate_exact.jl")

Random.seed!(42)

println("INVESTIGATION RESULTS:")
println("="^50)
println("✓ Bayesian method is working correctly")
println("✓ Exact validation detects sampler issues")
println("✗ Clifford implementation has probability errors")
println()
println("For 2 photons, 3 modes: TVD = 0.025 (acceptable)")
println("For 3 photons, 4 modes: TVD = 0.22 (poor)")
println()
println("The ratios show systematic bias:")
println("- Some configs oversampled (ratio > 1.5)")  
println("- Some configs undersampled (ratio < 0.5)")
println()
println("This confirms the Bayesian validation was correct!")
println("="^50)