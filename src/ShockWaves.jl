module ShockWaves

# An abstract type defining a shock wave
# Each of the files below will define a concrete type as a descendant of this type to specialize their functions
abstract type AbstractShockWave end

# Define a few convenience functions

"""
    SpeedOfSound( γ::T, P::T, ρ::T ) where { T <: AbstractFloat }

Computes the speed of sound given the supplied values of the ratio of specific heats, pressure, and density

# Arguments
    γ: Ratio of specific heats
    P: Pressure
    ρ: Density

# Returns
    A value of type T representing the speed of sound

# Equation
    This function calculates the speed of sound according to
          a = √( γ P / ρ )
    where
        ⋅ a is the speed of sound
        ⋅ γ is the ratio of specific heats
        ⋅ P is the pressure
        ⋅ ρ is the density

# Note
    This function is simply a rearrangement of the expression
          a = √( γ R T )
    where
        ⋅ R is the individual gas constant, equal to the universal gas constant divided by the molecular weight, R₀/μ
        ⋅ T is the temperature
"""
function SpeedOfSound( γ::T, P::T, ρ::T ) where { T <: AbstractFloat }
    return sqrt( γ * P / ρ )
end

include("NormalShockWave.jl")       # Functions to compute properties across a normal shock wave
include("MovingShockWave.jl")       # Functions to compute properties across a moving shock wave

export SpeedOfSound

end
