using ShockWaves
using Test

@testset "Helper functions" begin
    # Test the speed of sound function
    function TestHelperSpeedOfSound()
        test_cases = [
                      ( 1.40, 1.0, 1.0, 1.18322 ),
                      ( 1.40, 2.0, 1.0, 1.67332 ),
                      ( 1.40, 2.0, 2.0, 1.18322 ),
                      ( 1.40, 1.0, 2.0, 0.83666 ),
                      ( 1.67, 1.0, 1.0, 1.29228 ),
                      ( 1.67, 2.0, 1.0, 1.82757 ),
                      ( 1.67, 2.0, 2.0, 1.29228 ),
                      ( 1.67, 1.0, 2.0, 0.91378 ),
                     ]

        absolute_tolerance = 0.00001
        @testset "Speed of Sound" begin
            for ( γ, P, ρ, ExpectedValue ) in test_cases
                ActualValue = SpeedOfSound( γ, P, ρ )
                @test isapprox( ExpectedValue, ActualValue, atol=absolute_tolerance )
            end
        end
    end
    TestHelperSpeedOfSound()

end

# Tests for the NormalShockWave functions
@testset "Normal Shock Wave" begin
    # Test the PostShockMachNumber function
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
        @testset "Post Shock Mach Number" begin
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
        @testset "Density Ratio" begin
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
        @testset "Pressure Ratio" begin
            for ( γ, M₁, ExpectedValue ) in test_cases
                ActualValue = PressureRatio( nsw, γ, M₁ )
                @test isapprox( ExpectedValue, ActualValue, atol=absolute_tolerance )
            end
        end
    end
    TestNormalShockPressureRatio()

    # Test the calculation of Mach number from the pressure ratio
    function TestNormalShockMachNumberFromPressureRatio()
        nsw = NormalShockWave()

        # Inputs are γ, P₂P₁, ExpectedValue
        test_cases = [
                      ( 1.40, 1.2, 1.08233 ),
                      ( 1.40, 1.5, 1.19523 ),
                      ( 1.40, 2.5, 1.51186 ),
                      ( 1.67, 1.2, 1.07698 ),
                      ( 1.67, 1.5, 1.18309 ),
                      ( 1.67, 2.5, 1.48294 ),
                     ]
        
        absolute_tolerance = 0.00001
        @testset "Shock Wave Mach Number from Pressure Ratio" begin
            for ( γ, P₂P₁, ExpectedValue ) in test_cases
                # Test the function
                ActualValue = MachNumberFromPressureRatio( nsw, γ, P₂P₁ )
                @test isapprox( ExpectedValue, ActualValue, atol=absolute_tolerance )
            end
        end
    end
    TestNormalShockMachNumberFromPressureRatio()

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
        @testset "Temperature Ratio" begin
            for ( γ, M₁, ExpectedValue ) in test_cases
                ActualValue = TemperatureRatio( nsw, γ, M₁ )
                @test isapprox( ExpectedValue, ActualValue, atol=absolute_tolerance )
            end
        end
    end
    TestNormalShockTemperatureRatio()
end

# Tests for the MovingShockWave functions
@testset "Moving Shock Wave" begin
    # Test the DensityRatio function
    function TestMovingShockDensityRatio()
        msw = MovingShockWave()

        test_cases = [
                      ( 1.40, 1.1, 1.07042 ),
                      ( 1.40, 1.5, 1.33333 ),
                      ( 1.40, 2.5, 1.88235 ),
                      ( 1.67, 1.1, 1.05870 ),
                      ( 1.67, 1.5, 1.27211 ),
                      ( 1.67, 2.5, 1.69045 ),
                     ]
        absolute_tolerance = 0.00001
        @testset "Density Ratio" begin
            for ( γ, P₂P₁, ExpectedValue ) in test_cases
                ActualValue = DensityRatio( msw, γ, P₂P₁ )
                @test isapprox( ExpectedValue, ActualValue, atol=absolute_tolerance )
            end
        end
    end
    TestMovingShockDensityRatio()

    # Test the DensityRatioFromMachNumber function
    function TestMovingShockDensityRatioFromMachNumber()
        msw = MovingShockWave()

        test_cases = [
                      ( 1.40, 1.2, 1.34161 ),
                      ( 1.40, 1.5, 1.86207 ),
                      ( 1.40, 2.5, 3.33333 ),
                      ( 1.67, 1.2, 1.29682 ),
                      ( 1.67, 1.5, 1.71276 ),
                      ( 1.67, 2.5, 2.69697 ),
                     ]
        absolute_tolerance = 0.00001
        @testset "Density Ratio From Mach Number" begin
            for ( γ, M₁, ExpectedValue ) in test_cases
                ActualValue = DensityRatioFromMachNumber( msw, γ, M₁ )
                @test isapprox( ExpectedValue, ActualValue, atol=absolute_tolerance )
            end
        end
    end
    TestMovingShockDensityRatioFromMachNumber()

    # Test the TemperatureRatio function
    function TestMovingShockTemperatureRatio()
        msw = MovingShockWave()

        test_cases = [
                      ( 1.40, 1.1, 1.02763 ),
                      ( 1.40, 1.5, 1.12500 ),
                      ( 1.40, 2.5, 1.32813 ),
                      ( 1.67, 1.1, 1.03901 ),
                      ( 1.67, 1.5, 1.17914 ),
                      ( 1.67, 2.5, 1.47890 ),
                     ]
        absolute_tolerance = 0.00001
        @testset "Temperature Ratio" begin
            for ( γ, P₂P₁, ExpectedValue ) in test_cases
                ActualValue = TemperatureRatio( msw, γ, P₂P₁ )
                @test isapprox( ExpectedValue, ActualValue, atol=absolute_tolerance )
            end
        end
    end
    TestMovingShockTemperatureRatio()
    
    # Test the TemperatureRatioFromMachNumber function
    function TestMovingShockTemperatureRatioFromMachNumber()
        msw = MovingShockWave()

        test_cases = [
                      ( 1.40, 1.2, 1.12800 ),
                      ( 1.40, 1.5, 1.32021 ),
                      ( 1.40, 2.5, 2.13750 ),
                      ( 1.67, 1.2, 1.19555 ),
                      ( 1.67, 1.5, 1.49681 ),
                      ( 1.67, 2.5, 2.80590 ),
                     ]
        absolute_tolerance = 0.00001
        @testset "Temperature Ratio From Mach Number" begin
            for ( γ, M₁, ExpectedValue ) in test_cases
                ActualValue = TemperatureRatioFromMachNumber( msw, γ, M₁ )
                @test isapprox( ExpectedValue, ActualValue, atol=absolute_tolerance )
            end
        end
    end
    TestMovingShockTemperatureRatioFromMachNumber()

    # Test the ShockWaveVelocity function
    # This function tests both the version that accepts the upstream parameters as well as one that accepts the upstream speed of sound
    function TestMovingShockShockWaveVelocity()
        msw = MovingShockWave()


        # Inputs are γ, P₂P₁, P₁, ρ₁, ExpectedValue
        test_cases = [
                      ( 1.40, 1.1, 1.0, 1.0, 1.23288 ),
                      ( 1.40, 1.5, 1.0, 1.0, 1.41421 ),
                      ( 1.40, 1.1, 2.0, 1.0, 1.74356 ),
                      ( 1.40, 1.5, 2.0, 1.0, 2.00000 ),
                      ( 1.40, 1.1, 1.0, 2.0, 0.87178 ),
                      ( 1.40, 1.5, 1.0, 2.0, 1.00000 ),
                      ( 1.40, 1.1, 2.0, 2.0, 1.23288 ),
                      ( 1.40, 1.5, 2.0, 2.0, 1.41421 ),
                      ( 1.67, 1.1, 1.0, 1.0, 1.34294 ),
                      ( 1.67, 1.5, 1.0, 1.0, 1.52889 ),
                      ( 1.67, 1.1, 2.0, 1.0, 1.89921 ),
                      ( 1.67, 1.5, 2.0, 1.0, 2.16217 ),
                      ( 1.67, 1.1, 1.0, 2.0, 0.94961 ),
                      ( 1.67, 1.5, 1.0, 2.0, 1.08109 ),
                      ( 1.67, 1.1, 2.0, 2.0, 1.34294 ),
                      ( 1.67, 1.5, 2.0, 2.0, 1.52889 ),
                     ]
        
        absolute_tolerance = 0.00001
        @testset "Shock Wave Velocity From Pressure Ratio" begin
            for ( γ, P₂P₁, P₁, ρ₁, ExpectedValue ) in test_cases
                # First, check that the function itself is returning correct values
                # Since we tested SpeedOfSound already further up, we can assume that if those tests past then it's okay to use here too
                a₁ = SpeedOfSound( γ, P₁, ρ₁ )
                DirectValue = ShockWaveVelocity( msw, a₁, γ, P₂P₁ ) # The result of the function calling with a pre-calculated speed of sound. This ensures the function is returning the right value
                @test isapprox( ExpectedValue, DirectValue, atol=absolute_tolerance )
                # Test the version that takes in each property individually. This ensures that the speed of sound is being calculated/passed correctly in the function
                ActualValue = ShockWaveVelocity( msw, γ, P₁, ρ₁, P₂P₁ )
                @test isapprox( ExpectedValue, ActualValue, atol=absolute_tolerance )
            end
        end
    end
    TestMovingShockShockWaveVelocity()

    # Test the calculation of shock wave velocity from Mach number
    function TestMovingShockShockWaveVelocityFromMachNumber()
        msw = MovingShockWave()


        # Inputs are γ, M, P₁, ρ₁, ExpectedValue
        test_cases = [
                      ( 1.40, 1.04198, 1.0, 1.0, 1.23288 ),
                      ( 1.40, 1.19523, 1.0, 1.0, 1.41421 ),
                      ( 1.40, 1.04198, 2.0, 1.0, 1.74356 ),
                      ( 1.40, 1.19523, 2.0, 1.0, 2.00000 ),
                      ( 1.40, 1.04198, 1.0, 2.0, 0.87178 ),
                      ( 1.40, 1.19523, 1.0, 2.0, 1.00000 ),
                      ( 1.40, 1.04198, 2.0, 2.0, 1.23288 ),
                      ( 1.40, 1.19523, 2.0, 2.0, 1.41421 ),
                      ( 1.67, 1.03920, 1.0, 1.0, 1.34294 ),
                      ( 1.67, 1.18309, 1.0, 1.0, 1.52889 ),
                      ( 1.67, 1.03920, 2.0, 1.0, 1.89921 ),
                      ( 1.67, 1.18309, 2.0, 1.0, 2.16217 ),
                      ( 1.67, 1.03920, 1.0, 2.0, 0.94961 ),
                      ( 1.67, 1.18309, 1.0, 2.0, 1.08109 ),
                      ( 1.67, 1.03920, 2.0, 2.0, 1.34294 ),
                      ( 1.67, 1.18309, 2.0, 2.0, 1.52889 ),
                     ]
        
        absolute_tolerance = 0.00001
        @testset "Shock Wave Velocity From Mach Number" begin
            for ( γ, M, P₁, ρ₁, ExpectedValue ) in test_cases
                # First, check that the function itself is returning correct values
                # Since we tested SpeedOfSound already further up, we can assume that if those tests past then it's okay to use here too
                a₁ = SpeedOfSound( γ, P₁, ρ₁ )
                DirectValue = ShockWaveVelocityFromMachNumber( msw, a₁, γ, M ) # The result of the function calling with a pre-calculated speed of sound. This ensures the function is returning the right value
                @test isapprox( ExpectedValue, DirectValue, atol=absolute_tolerance )
                # Test the version that takes in each property individually. This ensures that the speed of sound is being calculated/passed correctly in the function
                ActualValue = ShockWaveVelocityFromMachNumber( msw, γ, P₁, ρ₁, M )
                @test isapprox( ExpectedValue, ActualValue, atol=absolute_tolerance )
            end
        end
    end
    TestMovingShockShockWaveVelocityFromMachNumber()

    # Test the calculation of shock wave Mach number from pressure ratio
    function TestMovingShockMachNumberFromPressureRatio()
        msw = MovingShockWave()

        # Inputs are γ, P₂P₁, ExpectedValue
        test_cases = [
                      ( 1.40, 1.2, 1.08233 ),
                      ( 1.40, 1.5, 1.19523 ),
                      ( 1.40, 2.5, 1.51186 ),
                      ( 1.67, 1.2, 1.07698 ),
                      ( 1.67, 1.5, 1.18309 ),
                      ( 1.67, 2.5, 1.48294 ),
                     ]
        
        absolute_tolerance = 0.00001
        @testset "Shock Wave Mach Number from Pressure Ratio" begin
            for ( γ, P₂P₁, ExpectedValue ) in test_cases
                # Test the function
                ActualValue = MachNumberFromPressureRatio( msw, γ, P₂P₁ )
                @test isapprox( ExpectedValue, ActualValue, atol=absolute_tolerance )
            end
        end
    end
    TestMovingShockMachNumberFromPressureRatio()
end
