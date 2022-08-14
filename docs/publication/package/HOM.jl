using BosonSampling

# Set the model of partial distinguishability
T = OneParameterInterpolation
# Define the unbalanced beams-plitter
B = BeamSplitter(1/sqrt(2))
# Set each particle in a different mode
r = ModeOccupation([1,1])

# Will store the coincidence probability
P_coinc = Vector{Float64}(undef,0)

for t in -4:0.01:4
    # distinguishability
    dist = exp(-t^2)
    i = Input{T}(r,dist)
    # Compute the full output mode occupation
    # statistics
    p_ =
    theoretical_distribution(input=i,interf=B)
    # Store the probabiltiy to observe [1,1]
    push!(proba_bunching, p_[2]+p_[3])
end
