"""
Test Householder sampler integration with UserDefinedGramMatrix
"""

@testset "Householder sampler with UserDefinedGramMatrix" begin

    # Parameters
    n = 3  # number of photons
    m = 5  # number of modes
    r = 2  # rank of Gram matrix

    # Generate random Gram matrix
    S = rand_gram_matrix_from_orthonormal_basis(n, r)

    # Create input state with partial distinguishability
    input = Input{UserDefinedGramMatrix}(first_modes(n, m), S)

    # Create interferometer
    interf = RandHaar(m)

    # Create event and sample
    output = FockSample()
    ev = Event(input, output, interf)
    sample!(ev)

    # Test photon number conservation
    sampled_output = ev.output_measurement.s
    @test sum(sampled_output.state) == n

    # Test multiple samples for consistency
    for i in 1:10
        ev_test = Event(input, FockSample(), interf)
        sample!(ev_test)
        @test sum(ev_test.output_measurement.s.state) == n
        @test length(ev_test.output_measurement.s.state) == m
    end

    # Test with different Gram matrix ranks
    for test_r in 1:min(n, 3)
        S_test = rand_gram_matrix_from_orthonormal_basis(n, test_r)
        input_test = Input{UserDefinedGramMatrix}(first_modes(n, m), S_test)
        ev_test = Event(input_test, FockSample(), interf)
        sample!(ev_test)
        @test sum(ev_test.output_measurement.s.state) == n
    end

end
