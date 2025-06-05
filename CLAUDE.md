# BosonSampling.jl Package Structure and Usage Guide

## Overview

BosonSampling.jl is a comprehensive Julia package for simulating boson sampling experiments with various physical effects including loss, partial distinguishability, and realistic detector imperfections. The package provides efficient algorithms for sampling and probability calculations, along with validation and certification tools for quantum advantage demonstrations.

## Package Architecture

### Core Dependencies

The package builds on several key Julia packages:
- **Mathematical**: `Permanents`, `LinearAlgebra`, `PolynomialRoots`, `Roots`
- **Statistical**: `Statistics`, `StatsBase`, `Distributions`, `HypothesisTests`
- **Data & I/O**: `CSV`, `DataFrames`, `JLD`, `DelimitedFiles`
- **Visualization**: `Plots`, `PrettyTables`, `Luxor`, `ColorSchemes`
- **Performance**: `BenchmarkTools`, `ProgressMeter`, `ProgressBars`
- **Optimization**: `Optim`

### Module Structure

```
BosonSampling.jl/
├── src/
│   ├── BosonSampling.jl          # Main module with exports
│   ├── main.jl                   # Core functionality
│   ├── types/                    # Type definitions
│   │   ├── types.jl             # All type includes
│   │   ├── input.jl             # Input state types
│   │   ├── interferometers.jl   # Interferometer types
│   │   ├── events.jl            # Event structure
│   │   ├── measurements.jl      # Output measurement types
│   │   └── mode_occupation.jl   # Mode representations
│   ├── boson_samplers/          # Sampling algorithms
│   │   ├── sample.jl            # Main sampling interface
│   │   ├── cliffords_sampler.jl # Clifford's algorithm
│   │   ├── partial_distinguishability_sampler.jl
│   │   ├── classical_sampler.jl
│   │   ├── noisy_sampler.jl
│   │   └── metropolis_sampler.jl
│   ├── distributions/           # Probability calculations
│   ├── certification/           # Bayesian validation
│   ├── partitions/             # Partition analysis
│   ├── loop/                   # Loop architectures
│   ├── bunching/              # Bunching analysis
│   └── visual.jl              # Visualization tools
```

## Type System

### 1. Input States (`src/types/input.jl`)

The package uses a sophisticated type hierarchy to represent different levels of photon distinguishability:

#### Base Types
```julia
abstract type InputType end

# Concrete input types
struct Bosonic <: InputType end                    # Indistinguishable photons
struct Distinguishable <: InputType end           # Fully distinguishable photons
abstract type PartDist <: InputType end           # Partial distinguishability

# Partial distinguishability models
struct OneParameterInterpolation <: PartDist      # Single parameter x ∈ [0,1]
struct RandomGramMatrix <: PartDist               # Random Gram matrix
struct UserDefinedGramMatrix <: PartDist          # Custom Gram matrix
```

#### Input Construction
```julia
# Basic input creation
Input{Bosonic}(first_modes(n, m))                 # n photons in first n modes
Input{Distinguishable}(last_modes(n, m))          # n photons in last n modes

# Partial distinguishability with parameter
Input{OneParameterInterpolation}(first_modes(n, m), x)  # x=0: distinguishable, x=1: bosonic

# Custom Gram matrix
S = your_gram_matrix  # n×n positive semi-definite matrix
Input{UserDefinedGramMatrix}(first_modes(n, m), S)
```

### 2. Mode Representations (`src/types/mode_occupation.jl`)

```julia
# Mode occupation: photon count per mode
ModeOccupation([1, 2, 0, 1])  # 4 photons in modes [1,2,4]

# Mode list: which mode each photon occupies
ModeList([1, 2, 2, 4])        # Same configuration as above

# Conversion utilities
first_modes(n, m)             # n photons in first n modes
last_modes(n, m)              # n photons in last n modes
convert(ModeList, mo)         # Convert representations
```

### 3. Interferometers (`src/types/interferometers.jl`)

```julia
# Random Haar-distributed unitary
RandHaar(m)                   # m×m random unitary

# Discrete Fourier Transform
Fourier(m)                    # DFT matrix

# Hadamard matrix
Hadamard(m)                   # Hadamard matrix (m must be power of 2)

# User-defined unitary
UserDefinedInterferometer(U)   # Custom m×m unitary matrix U

# With loss (trait-based)
add_loss(interferometer, η)    # Add uniform loss with efficiency η
```

### 4. Events (`src/types/events.jl`)

The Event structure connects all components of a boson sampling experiment:

```julia
# Basic event creation
ev = Event(input_state, output_measurement, interferometer)

# With additional parameters
ev = Event{TIn, TOut}(input, output, interf, proba_params)

# Example complete setup
input = Input{Bosonic}(first_modes(3, 5))
output = FockSample()                    # Unknown output to be sampled
interf = RandHaar(5)
ev = Event(input, output, interf)
```

### 5. Output Measurements (`src/types/measurements.jl`)

```julia
# Sampling (unknown outcome)
FockSample()                             # To be filled by sampling
DarkCountFockSample(p_dark)             # With dark count probability
RealisticDetectorsFockSample(p_dark, p_no_count)  # Full detector model

# Detection (known outcome)
FockDetection(ModeOccupation([1,1,0,1,0]))  # Specific detection pattern

# Threshold detection
ThresholdFockDetection(ThresholdModeOccupation([1,0,1,0]))
```

## Usage Patterns

### Basic Boson Sampling

```julia
using BosonSampling

# Setup experiment
n, m = 3, 6                              # 3 photons, 6 modes
input = Input{Bosonic}(first_modes(n, m))
interf = RandHaar(m)
output = FockSample()

# Create and run experiment
ev = Event(input, output, interf)
sample!(ev)

# Access result
sampled_state = ev.output_measurement.s
println("Sampled output: ", sampled_state)
```

### Comparing Different Statistics

```julia
# Setup inputs with different statistics
input_bosonic = Input{Bosonic}(first_modes(n, m))
input_distinguishable = Input{Distinguishable}(first_modes(n, m))
input_partial = Input{OneParameterInterpolation}(first_modes(n, m), 0.8)

# Same interferometer for fair comparison
interf = UserDefinedInterferometer(your_unitary)

# Sample from each
ev_b = Event(input_bosonic, FockSample(), interf)
ev_d = Event(input_distinguishable, FockSample(), interf)
ev_p = Event(input_partial, FockSample(), interf)

sample!(ev_b)
sample!(ev_d)  
sample!(ev_p)

println("Bosonic: ", ev_b.output_measurement.s)
println("Distinguishable: ", ev_d.output_measurement.s)
println("Partial (x=0.8): ", ev_p.output_measurement.s)
```

### Probability Calculations

```julia
# Exact probability calculation
output_state = ModeOccupation([1, 0, 1, 1, 0, 0])
ev_prob = Event(input, FockDetection(output_state), interf)
compute_probability!(ev_prob)
probability = ev_prob.proba_params.probability

# Full distribution (for small systems)
dist = full_distribution(input, interf)
println("Distribution: ", dist)
```

### Including Physical Effects

#### Loss
```julia
η = 0.85  # 85% efficiency
sample!(ev, η)  # Include uniform loss
```

#### Dark Counts
```julia
p_dark = 0.01  # 1% dark count probability  
output = DarkCountFockSample(p_dark)
ev = Event(input, output, interf)
sample!(ev)
```

#### Realistic Detectors
```julia
p_dark = 0.01     # Dark count probability
p_no_count = 0.05 # No-count probability
output = RealisticDetectorsFockSample(p_dark, p_no_count)
ev = Event(input, output, interf)
sample!(ev)
```

## Sampling Algorithms

### 1. Clifford's Algorithm (Corrected)

For indistinguishable bosons (`Input{Bosonic}`):

```julia
# Automatic dispatch through sample!
ev = Event(Input{Bosonic}(first_modes(n, m)), FockSample(), interf)
sample!(ev)

# Direct usage
sample = cliffords_sampler(input=input_state, interf=interferometer)
```

**Features:**
- Fixes bias from original Clifford & Clifford algorithm
- Efficient permanent calculation using Ryser algorithm
- Handles repeated rows with proper multiplicities

### 2. Classical Sampler

For distinguishable photons (`Input{Distinguishable}`):

```julia
ev = Event(Input{Distinguishable}(first_modes(n, m)), FockSample(), interf)
sample!(ev)
```

**Algorithm:** Multinomial sampling based on |U_{i,j}|² probabilities.

### 3. Partial Distinguishability Extension

For partially distinguishable photons (`Input{OneParameterInterpolation}`):

```julia
x = 0.7  # Distinguishability parameter
input = Input{OneParameterInterpolation}(first_modes(n, m), x)
ev = Event(input, FockSample(), interf)
sample!(ev)
```

**Features:**
- Implements algorithm from arXiv:2406.17682v2
- Polynomial complexity for bounded interference order
- Supports heterogeneous distinguishability models

### 4. Noisy Sampler

Includes various physical imperfections:

```julia
# Setup noisy sampling
ev = Event(input, FockSample(), interf)

# With loss and dark counts
η = 0.9
p_dark = 0.01
sample!(ev, η, p_dark)
```

## Advanced Features

### Partition Analysis

For analyzing photon distributions in mode subsets:

```julia
using BosonSampling

# Define partition (list of mode subsets)
partition = [[1, 2], [3, 4], [5, 6]]

# Compute partition probabilities
(indexes, probs) = compute_probabilities_partition(interf, partition, input)

# Event-based partition measurement
part_event = Event(input, PartitionCountsOutput(partition), interf)
compute_probability!(part_event)
```

### Loop Architectures

For thermalization and pseudo-number resolution:

```julia
# Thermalization setup
η_vals = η_thermalization(n)          # Get transmissivities
partition = partition_thermalization(m) # Create partition

# Pseudo photon number resolution
steps = 3
η_pnr_vals = η_pnr(steps)            # [0.25, 0.5, 0.75]
```

### Bayesian Certification

For validating quantum advantage:

```julia
# Define hypotheses
p_bosonic = HypothesisFunction(p_B)      # Bosonic hypothesis
p_classical = HypothesisFunction(p_D)    # Classical hypothesis

# Collect experimental events
events = [sample_experiment() for _ in 1:num_samples]

# Compute confidence
confidence = compute_confidence(events, p_bosonic, p_classical)
println("Confidence in bosonic behavior: ", confidence)

# Estimate required sample size
samples_needed = number_of_samples(ev_bosonic, ev_classical, p_null=0.95)
```

### Visualization

```julia
# Visualize sampling setup
visualize_sampling(input, sampled_output)

# Visualize with probability
prob_data = full_distribution(input, interf)
visualize_proba(input, specific_output, prob_data)
```

## Performance Considerations

### Algorithm Complexity
- **Clifford's Algorithm**: O(n²m + permanent calculation)
- **Classical Sampler**: O(nm)  
- **Partial Distinguishability**: O(k^n) where k is interference order
- **Partition Methods**: Scale better than full distributions for large m

### Optimization Tips
1. **Use appropriate input types**: Automatic dispatch to optimal algorithms
2. **Partition analysis**: More efficient than full distributions for large systems
3. **Permanent calculation**: Ryser algorithm is default (fastest for moderate sizes)
4. **Memory management**: Use `clean_proba` and `clean_pdf` for numerical stability

### Benchmarking
```julia
using BenchmarkTools

# Benchmark sampling
@benchmark sample!(ev) setup=(ev = Event(input, FockSample(), interf))

# Compare algorithms
@benchmark cliffords_sampler(input=input, interf=interf)
@benchmark classical_sampler(input=input_d, interf=interf)
```

## Common Workflows

### 1. Quantum Advantage Demonstration
```julia
# Setup experiment
n, m = 4, 8
input_b = Input{Bosonic}(first_modes(n, m))
input_d = Input{Distinguishable}(first_modes(n, m))
interf = RandHaar(m)

# Collect samples
num_samples = 1000
samples_b = [sample!(Event(input_b, FockSample(), interf)).output_measurement.s for _ in 1:num_samples]
samples_d = [sample!(Event(input_d, FockSample(), interf)).output_measurement.s for _ in 1:num_samples]

# Statistical validation
tvd_distance = tvd(histogram(samples_b), histogram(samples_d))
```

### 2. Device Characterization
```julia
# Measure interferometer with known input
known_input = Input{Bosonic}(first_modes(n, m))
experimental_data = collect_experimental_samples(device, known_input, num_samples)

# Compare with theory
theoretical_dist = full_distribution(known_input, estimated_interferometer)
experimental_dist = histogram(experimental_data)

# Compute fidelity
fidelity = 1 - tvd(theoretical_dist, experimental_dist)/2
```

### 3. Loss Characterization
```julia
# Vary loss parameter
η_values = 0.5:0.1:1.0
tvd_vs_loss = []

for η in η_values
    # Sample with loss
    lossy_samples = [sample!(Event(input, FockSample(), interf), η) for _ in 1:num_samples]
    
    # Compare to lossless
    lossless_samples = [sample!(Event(input, FockSample(), interf)) for _ in 1:num_samples]
    
    # Compute distance
    distance = tvd(histogram(lossy_samples), histogram(lossless_samples))
    push!(tvd_vs_loss, distance)
end

# Plot results
plot(η_values, tvd_vs_loss, xlabel="Loss parameter η", ylabel="TVD from lossless")
```

## Error Handling and Validation

### Input Validation
```julia
# Validate unitary matrix
@assert is_unitary(U, atol=1e-10) "Matrix must be unitary"

# Validate Gram matrix  
@assert is_gram_matrix(S) "Matrix must be positive semi-definite"

# Validate probabilities
clean_prob = clean_proba(probability_value)  # Handles numerical errors
normalized_dist = clean_pdf(distribution)    # Normalizes and validates
```

### Common Errors and Solutions

1. **Non-unitary interferometer**: Use `is_unitary()` to check matrices
2. **Invalid Gram matrix**: Ensure positive semi-definite for partial distinguishability
3. **Probability normalization**: Use `clean_pdf()` for numerical issues
4. **Mode conservation**: Check input/output photon numbers match (with loss)

## Testing and Development

### Running Tests
```julia
using Pkg
Pkg.test("BosonSampling")
```

### Example Usage Scripts
- `test/example_usage.jl`: Basic usage examples
- `benchmarks/benchmarks.jl`: Performance benchmarks
- `docs/publication/`: Publication-quality examples

This documentation provides a comprehensive guide to using BosonSampling.jl effectively for quantum optics simulations, algorithm development, and experimental validation.