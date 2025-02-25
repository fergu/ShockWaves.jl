# These are functions to compute properties across a normal shock wave
# These relationships are taken from
# Anderson, John D. "Modern Compressible Flow with Historical Perspective" 3rd ed. (2003) McGraw-Hill. ISBN: O-01-242443-5
# Specifically, these equations are from chapter 3

"""
    NormalShockWave <: AbstractShockWave

A type describing a normal shock wave.
This type is primarily meant to help specialize function calls for normal shock waves (as opposed to other types of shocks)
"""
struct NormalShockWave <: AbstractShockWave end

"""
    PostShockMachNumber( nsw::NormalShockWave, γ::T, M₁::T ) where { T <: AbstractFloat }

Computes the post-shock Mach number following a normal shock wave

# Arguments
    nsw: A NormalShockWave() type
    γ: The ratio of specific heats for the gas the shock wave is in
    M₁: The Mach number of the shock wave, must be greater than 1.0

# Returns
    A value of type T representing the post-shock Mach number

# Equation
    The post-shock Mach number is calculated using Eq. 3.51 as:
          M₂² = (1 + [(γ - 1)/2]M₁²)/(γM₁² - (γ-1)/2)
    where
        ⋅ M₁ is the pre-shock Mach number
        ⋅ M₂ is the post-shock Mach number
        ⋅ γ  is the ratio of specific heats for the gas
    
# Note
    For convenience, this function returns the Mach number directly rather than its square. That is, the returned value is the square root of Eq. 3.51.
"""
function PostShockMachNumber( nsw::NormalShockWave, γ::T, M₁::T ) where { T <: AbstractFloat }
    return sqrt( ( 1 + ( ( γ - 1.0 ) / 2.0 ) * M₁^2 ) / ( γ * M₁^2 - ( γ - 1.0 ) / 2.0 ) )
end

"""
    DensityRatio( nsw::NormalShockWave, γ::T, M₁::T ) where { T <: AbstractFloat }

Computes the density ratio ρ₂/ρ₁ across a normal shock wave, where state 1 is before the shock and state 2 is after the shock.

# Arguments
    nsw: A NormalShockWave() type
    γ: The ratio of specific heats of the gas the shock wave is in
    M₁: The Mach number of the shock wave, must be greater than 1.0

# Returns
    A value of type T representing the density ratio across the shock wave, ρ₂/ρ₁

# Equation
    The density ratio is computed using Eq. 3.53 as:
          ρ₂/ρ₁ = ((γ + 1)M₁²) / (2 + (γ - 1)M₁²)
    where
        ⋅ ρ₁ is the density before the shock wave
        ⋅ ρ₂ is the density after the shock wave
        ⋅ γ  is the ratio of specific heats for the gas
        ⋅ M₁ is the Mach number of the shock wave
"""
function DensityRatio( nsw::NormalShockWave, γ::T, M₁::T ) where { T <: AbstractFloat }
    return ( ( γ + 1.0 ) * M₁^2 ) / ( 2.0 + ( γ - 1.0 ) * M₁^2 )
end

"""
    PressureRatio( nsw::NormalShockWave, γ::T, M₁::T ) where { T <: AbstractFloat }

Computes the pressure ratio, P₂/P₁ across a normal shock wave, where state 1 is before the shock and state 2 is after the shock.

# Arguments
    nsw: A NormalShockWave type
    γ: The ratio of specific heats of the gas the shock wave is in
    M₁: The Mach number of the shock wave, must be greater than 1.0

# Returns
    A value of type T representing the pressure ratio across the shock wave, P₂/P₁

# Equation
    The pressure ratio is computed using Eq. 3.57 as:
          P₂/P₁ = 1 + (2γ)/(γ+1) (M₁² - 1)
    where
        ⋅ P₁ is the pressure before the shock wave
        ⋅ P₂ is the pressure after the shock wave
        ⋅ γ  is the ratio of specific heats for the gas
        ⋅ M₁ is the Mach number of the shock wave
"""
function PressureRatio( nsw::NormalShockWave, γ::T, M₁::T ) where { T <: AbstractFloat }
    return 1.0 + ( 2 * γ ) / ( γ + 1 ) * ( M₁^2 - 1 )
end

"""
    TemperatureRatio( nsw::NormalShockWave, γ::T, M₁::T ) where { T <: AbstractFloat }

Computes the temperature ratio, T₂/T₁, across a normal shock wave, where state 1 is before the shock and state 2 is after the shock.

# Arguments
    ::NormalShockWave: A NormalShockWave type
    γ: The ratio of specific heats of the gas the shock wave is in
    M₁: The Mach number of the shock wave, must be greater than 1.0

# Returns
    A value of type T representing the temperature ratio across the shock wave, T₂/T₁

# Equation
    The pressure ratio is computed using Eq. 3.58 as:
          T₂/T₁ = (P₂/P₁)(ρ₁/ρ₂)
    where
        ⋅ T₁ is the temperature before the shock wave
        ⋅ T₂ is the temperature after the shock wave
        ⋅ P₁ is the pressure before the shock wave
        ⋅ P₂ is the pressure after the shock wave
        ⋅ ρ₁ is the density before the shock wave
        ⋅ ρ₂ is the density after the shock wave

# Note
    Anderson substitutes the expressions for ρ₂/ρ₁ (Eq. 3.53) and P₂/P₁ (Eq. 3.57) into Eq. 3.58 to obtain an analytic expression for T₂/T₁ (Eq. 3.59) in terms of the Mach number and γ.
    However, as an implementation detail, this function calls the DensityRatio and PressureRatio functions and computes Eq. 3.58 from those results rather than implementing Eq. 3.59.
    This is done to ensure consistency with Eqs. 3.53 and 3.57.
"""
function TemperatureRatio( nsw::NormalShockWave, γ::T, M₁::T ) where { T <: AbstractFloat }
    return PressureRatio( nsw, γ, M₁ ) / DensityRatio( nsw, γ, M₁ )
end

export NormalShockWave, PostShockMachNumber, DensityRatio, PressureRatio, TemperatureRatio
