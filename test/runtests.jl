using ShockWaves
using Test

@testset "Normal Shock Wave" begin
    # Write your tests here.
    
    ## Tests for the NormalShockWave functions

    # Test the PostShockMachNumber
    function TestNormalShockPostShockMachNumber()
        nsw = NormalShockWave()

        test_cases = [
                      ( 1.40, 1.2, 0.84217 ),
                      ( 1.40, 1.5, 0.70109 ),
                      ( 1.40, 2.5, 0.51299 ),
                      ( 1.67, 1.2, 0.84629 ),
                      ( 1.67, 1.5, 0.71583 ),
                      ( 1.67, 2.5, 0.55339 ),
                     ]

        absolute_tolerance = 0.00001
        @testset "Post Shock Mach Number Tests" begin
            for ( γ, M₁, ExpectedValue ) in test_cases
                ActualValue = PostShockMachNumber( nsw, γ, M₁ )
                @test isapprox( ExpectedValue, ActualValue, atol=absolute_tolerance )
            end
        end
    end
    TestNormalShockPostShockMachNumber()

    # Test the DensityRatio function
    function TestNormalShockDensityRatio()
        nsw = NormalShockWave()

        test_cases = [
                      ( 1.40, 1.2, 1.34161 ),
                      ( 1.40, 1.5, 1.86207 ),
                      ( 1.40, 2.5, 3.33333 ),
                      ( 1.67, 1.2, 1.29682 ),
                      ( 1.67, 1.5, 1.71276 ),
                      ( 1.67, 2.5, 2.69697 ),
                     ]
        absolute_tolerance = 0.00001
        @testset "Density Ratio Tests" begin
            for ( γ, M₁, ExpectedValue ) in test_cases
                ActualValue = DensityRatio( nsw, γ, M₁ )
                @test isapprox( ExpectedValue, ActualValue, atol=absolute_tolerance )
            end
        end
    end
    TestNormalShockDensityRatio()

    # Test the PressureRatio function
    function TestNormalShockPressureRatio()
        nsw = NormalShockWave()

        test_cases = [
                      ( 1.40, 1.2, 1.51333 ),
                      ( 1.40, 1.5, 2.45833 ),
                      ( 1.40, 2.5, 7.12500 ),
                      ( 1.67, 1.2, 1.55041 ),
                      ( 1.67, 1.5, 2.56367 ),
                      ( 1.67, 2.5, 7.56742 ),
                     ]
        absolute_tolerance = 0.00001
        @testset "Pressure Ratio Tests" begin
            for ( γ, M₁, ExpectedValue ) in test_cases
                ActualValue = PressureRatio( nsw, γ, M₁ )
                @test isapprox( ExpectedValue, ActualValue, atol=absolute_tolerance )
            end
        end
    end
    TestNormalShockPressureRatio()

    # Test the TemperatureRatio function
    function TestNormalShockTemperatureRatio()
        nsw = NormalShockWave()

        test_cases = [
                      ( 1.40, 1.2, 1.12800 ),
                      ( 1.40, 1.5, 1.32021 ),
                      ( 1.40, 2.5, 2.13750 ),
                      ( 1.67, 1.2, 1.19555 ),
                      ( 1.67, 1.5, 1.49681 ),
                      ( 1.67, 2.5, 2.80590 ),
                     ]
        absolute_tolerance = 0.00001
        @testset "Pressure Ratio Tests" begin
            for ( γ, M₁, ExpectedValue ) in test_cases
                ActualValue = TemperatureRatio( nsw, γ, M₁ )
                @test isapprox( ExpectedValue, ActualValue, atol=absolute_tolerance )
            end
        end
    end
    TestNormalShockTemperatureRatio()
end
