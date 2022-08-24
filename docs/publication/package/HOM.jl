using BosonSampling
using Plots

# Set the model of partial distinguishability
T = OneParameterInterpolation
# Define the unbalanced beams-plitter
B = BeamSplitter(1/sqrt(2))
# Set each particle in a different mode
r = ModeOccupation([4,4])

# Will store the coincidence probability
P_coinc = []
p_ = []

for t in -3:0.01:3
    # distinguishability
    dist = exp(-t^2)
    i = Input{T}(r,dist)
    # Compute the full output mode occupation
    # statistics
    p_ =
    theoretical_distribution(input=i,interf=B)
    # Store the probabiltiy to observe [1,1]
    push!(P_coinc, p_[2]+p_[3])
end

plot!(P_coinc)
