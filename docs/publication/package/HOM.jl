using BosonSampling

# Set experimental parameters
Δω = 1
# Set the model of partial distinguishability
T = OneParameterInterpolation
# Define the unbalanced beams-plitter
B = BeamSplitter(1/sqrt(2))
# Set each particle in a different mode
r_i = ModeOccupation([1,1])

# Define the output as detecting a coincidence
r_f = ModeOccupation([1,1])
o = FockDetection(r_f)

# Will store the events probability
events = []

for Δt in -4:0.01:4
    # distinguishability
    dist = exp(-(Δω * Δt)^2)
    i = Input{T}(r_i,dist)

    # Create the event
    ev = Event(i,o,B)
    # Compute its probability to occur
    compute_probability!(ev)

    # Store the event and its probability
    push!(events, ev)
end
