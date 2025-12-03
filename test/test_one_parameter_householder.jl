"""
Test OneParameterInterpolation with Householder sampler
"""

@testset "OneParameterInterpolation with Householder sampler" begin

    # Parameters
    n = 3  # number of photons
    m = 5  # number of modes

    # Create interferometer (same for all tests)
    interf = RandHaar(m)

    # Test different distinguishability parameters
    x_values = [0.0, 0.3, 0.7, 1.0]  # 0=distinguishable, 1=indistinguishable

    for x in x_values
        # Create input with OneParameterInterpolation
        input = Input{OneParameterInterpolation}(first_modes(n, m), x)

        # Verify Gram matrix properties
        @test size(input.G.S) == (n, n)
        @test all(input.G.S[i,i] ≈ 1.0 for i in 1:n)  # Diagonal should be 1
        @test input.G.distinguishability_param == x

        # Sample and verify photon conservation
        for i in 1:5
            ev = Event(input, FockSample(), interf)
            sample!(ev)
            output = ev.output_measurement.s
            @test sum(output.state) == n
            @test length(output.state) == m
        end
    end

    # Test edge cases
    @testset "Edge case: x=0 (fully distinguishable)" begin
        input_dist = Input{OneParameterInterpolation}(first_modes(n, m), 0.0)
        # Should be identity matrix
        @test all(input_dist.G.S[i,j] ≈ (i==j ? 1.0 : 0.0) for i in 1:n, j in 1:n)

        ev = Event(input_dist, FockSample(), interf)
        sample!(ev)
        @test sum(ev.output_measurement.s.state) == n
    end

    @testset "Edge case: x=1 (fully indistinguishable)" begin
        input_indist = Input{OneParameterInterpolation}(first_modes(n, m), 1.0)
        @test all(input_indist.G.S .≈ 1.0)  # Should be all ones

        ev = Event(input_indist, FockSample(), interf)
        sample!(ev)
        @test sum(ev.output_measurement.s.state) == n
    end

    # Compare with UserDefinedGramMatrix (should give equivalent results)
    @testset "Equivalence with UserDefinedGramMatrix" begin
        x = 0.5
        input_oneparam = Input{OneParameterInterpolation}(first_modes(n, m), x)

        # Create equivalent UserDefinedGramMatrix
        S_manual = fill(x, n, n)
        for i in 1:n
            S_manual[i, i] = 1.0
        end
        input_userdefined = Input{UserDefinedGramMatrix}(first_modes(n, m), S_manual)

        # Both should have same Gram matrix
        @test input_oneparam.G.S ≈ input_userdefined.G.S

        # Both should sample successfully
        ev1 = Event(input_oneparam, FockSample(), interf)
        ev2 = Event(input_userdefined, FockSample(), interf)
        sample!(ev1)
        sample!(ev2)
        @test sum(ev1.output_measurement.s.state) == n
        @test sum(ev2.output_measurement.s.state) == n
    end

end
