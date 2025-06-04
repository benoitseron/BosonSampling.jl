# Clifford Sampler Development

## Files

- `clifford_final.jl` - Corrected implementation based on the existing codebase algorithm
- `playground.jl` - Simple test workspace with Bayesian validation

## Status

The Clifford implementation follows the same algorithm as the built-in `clifford_sampler_unoptimised` but both seem to have issues with Bayesian validation (getting 0.0 confidence). This may indicate:

1. Issues with the validation setup
2. Problems with the interferometer configuration  
3. Systematic issues in event creation

## Usage

```julia
include("clifford_final.jl")

# Test the implementation
input = Input{Bosonic}(first_modes(3, 6))
interf = RandHaar(6) 
sample = clifford_sampler_input_final(input, interf)
```

## Note

The algorithm implementation appears correct but requires further investigation into why Bayesian validation fails.