begin
    using Revise

    using BosonSampling
    using Plots
    using ProgressMeter
    using Distributions
    using Random
    using Test
    using ArgCheck
    using StatsBase
    using ColorSchemes
    using Interpolations
    using Dierckx
    using LinearAlgebra
    using PrettyTables
    using LaTeXStrings
    using JLD
    using AutoHashEquals
    using LinearRegression
    
    using DataStructures
    using FFTW
    # include package using @with_kw
    using Parameters

end

### test case with one mode ###

@testset "gaussian partition" begin
    
    @testset "single mode" begin
            
            
        """

        pdf_one_mode(n,r)

        Returns the probability of having n photons in one mode, given the displacement r. Input of 1 mode squeezed sttate.
        """
        pdf_one_mode(n,r::Real) = begin
                if n % 2 == 1
                        return 0.0
                else

                        μ = cosh(r)
                        ν = sinh(r)

                        return convert(Float64,1/μ * (ν/μ)^(n) * factorial(big(n)) / (2^(div(n,2)) * factorial(big(div(n,2))))^2)
                end
        end

        m = 1
        input_state = GeneralGaussian(m = m, r = 0.3 * ones(m))
        interferometer = Fourier(m)
        part = equilibrated_partition(m, 1)

        # part = Partition(Subset(first_modes(1, m)))
        # part.subsets[1].subset

        n_max = 51
        mc = compute_probabilities_partition_gaussian(interferometer, part, input_state, n_max)

        pdf = mc.proba

        pdf_one_mode_array = [pdf_one_mode(i,input_state.r[1]) for i in 0:n_max]

        # bar(pdf, alpha = 0.5)
        # bar!(pdf_one_mode_array, alpha = 0.5)

        # bar(pdf .- pdf_one_mode_array)

        if length(pdf) != length(pdf_one_mode_array)
            error("lengths are not equal")
        else

            #bar(pdf .- pdf_one_mode_array)

            @test pdf ≈ pdf_one_mode_array atol = 1e-10
        end



    end

end

@testset "gaussian partition - chicago algorithm" begin
    
    @testset "single mode" begin
            
            
        """

        pdf_one_mode(n,r)

        Returns the probability of having n photons in one mode, given the displacement r. Input of 1 mode squeezed sttate.
        """
        pdf_one_mode(n,r::Real) = begin
                if n % 2 == 1
                        return 0.0
                else

                        μ = cosh(r)
                        ν = sinh(r)

                        return convert(Float64,1/μ * (ν/μ)^(n) * factorial(big(n)) / (2^(div(n,2)) * factorial(big(div(n,2))))^2)
                end
        end

        m = 1
        input_state = GeneralGaussian(m = m, r = 0.3 * ones(m))
        interferometer = Fourier(m)
        part = equilibrated_partition(m, 1)

        # part = Partition(Subset(first_modes(1, m)))
        # part.subsets[1].subset

        n_max = 51
        mc = compute_probabilities_partition_gaussian_chicago(interferometer, part, input_state, n_max)

        pdf = mc.proba

        bar(pdf, alpha = 0.5)

        pdf_one_mode_array = [pdf_one_mode(i,input_state.r[1]) for i in 0:n_max]

        # bar(pdf, alpha = 0.5)
        # bar!(pdf_one_mode_array, alpha = 0.5)

        if length(pdf) != length(pdf_one_mode_array)
            error("lengths are not equal")
        else

            #bar(pdf .- pdf_one_mode_array)

            @test pdf ≈ pdf_one_mode_array atol = 1e-10
        end

    end

end



