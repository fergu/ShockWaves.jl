# These are relationships for moving (normal) shock waves, from Anderson's Modern Compressible Flow Ch 7

"""
    MovingShockWave <: AbstractShockWave

A type describing a moving shock wave.
This type is primarily meant to help specialize function calls for moving shock waves (as opposed to other types of shocks)
These shocks are treated as normal shock waves for most thermodynamic properties. 
The main difference is that the equations are expressed in terms of shock wave pressure ratios rather than a Mach number
"""
struct MovingShockWave <: AbstractShockWave end

"""
    DensityRatio( msw::MovingShockWave, γ::T, P₂P₁::T ) where { T <: AbstractFloat }

Computes the density ratio ρ₂/ρ₁ across a moving shock wave, where state 1 is before the shock and state 2 is after the shock.

# Arguments
    msw: A MovingShockWave() type
    γ: The ratio of specific heats of the gas the shock wave is propagating into
    P₂P₁: The pressure ratio of the shock wave, must be greater than 1.0

# Returns
    A value of type T representing the density ratio across the shock wave, ρ₂/ρ₁

# Equation
    The density ratio is computed using Eq. 7.11 as:
          ρ₂/ρ₁ = (1 + (γ+1)/(γ-1) + P₂/P₁) / (1 + (γ+1)/(γ-1)*P₂/P₁)
    where
        ⋅ ρ₁ is the density before the shock wave
        ⋅ ρ₂ is the density after the shock wave
        ⋅ γ  is the ratio of specific heats for the gas
        ⋅ P₂/P₁ is the pressure ratio across the shock wave
"""
function DensityRatio( msw::MovingShockWave, γ::T, P₂P₁::T ) where { T <: AbstractFloat }
    return (1.0 + ((γ+1.0) / (γ-1.0)) * P₂P₁)/( ((γ+1.0)/(γ-1.0)) + P₂P₁)
end

"""
    DensityRatioFromMachNumber( msw::MovingShockWave, γ::T, M::T ) where { T <: AbstractFloat }

Computes the density ratio ρ₂/ρ₁ across a moving shock wave based on the shock wave Mach number, where state 1 is before the shock and state 2 is after the shock.

# Arguments
    msw: A MovingShockWave() type
    γ: The ratio of specific heats of the gas the shock wave is propagating into
    M: The Mach number of the shock wave, must be greater than 1.0

# Returns
    A value of type T representing the density ratio across the shock wave, ρ₂/ρ₁

# Equation
    The density ratio is computed using Eq. 7.11 as:
          ρ₂/ρ₁ = (1 + (γ+1)/(γ-1) + P₂/P₁) / (1 + (γ+1)/(γ-1)*P₂/P₁)
    where
        ⋅ ρ₁ is the density before the shock wave
        ⋅ ρ₂ is the density after the shock wave
        ⋅ γ  is the ratio of specific heats for the gas
        ⋅ P₂/P₁ is the pressure ratio across the shock wave

# Note
    This function is a convenience wrapper around DensityRatio( ::MovingShockWave, γ, P₂P₁ ) to compute the density ratio in terms of a Mach number instead of a pressure ratio.
    The function PressureRatio( NormalShockWave(), γ, M ) is called to determine the pressure ratio from the Mach number, and this result is passed to DensityRatio( msw, γ, P₂P₁ )
"""
function DensityRatioFromMachNumber( msw::MovingShockWave, γ::T, M::T ) where { T <: AbstractFloat }
    P₂P₁ = PressureRatio( NormalShockWave(), γ, M )
    return DensityRatio( msw, γ, P₂P₁ )
end

"""
    TemperatureRatio( msw::MovingShockWave, γ::T, P₂P₁::T ) where { T <: AbstractFloat }

Computes the temperature ratio T₂/T₁ across a moving shock wave, where state 1 is before the shock and state 2 is after the shock.
The shock wave strength is based on the supplied pressure ratio across the shock wave, P₂P₁

# Arguments
    msw: A MovingShockWave() type
    γ: The ratio of specific heats of the gas the shock wave is propagating into
    P₂P₁: The pressure ratio of the shock wave, must be greater than 1.0

# Returns
    A value of type T representing the temperature ratio across the shock wave, T₂/T₁

# Equation
    The temperature ratio is computed using Eq. 7.10 as:
          T₂/T₁ = P₂/P₁ ( ((γ+1)/(γ-1)) + P₂/P₁ ) / ( 1 + ((γ+1)/(γ-1))*P₂/P₁ )
    where
        ⋅ T₁ is the temperature before the shock wave
        ⋅ T₂ is the temperature after the shock wave
        ⋅ γ  is the ratio of specific heats for the gas
        ⋅ P₂/P₁ is the pressure ratio across the shock wave
"""
function TemperatureRatio( msw::MovingShockWave, γ::T, P₂P₁::T ) where { T <: AbstractFloat }
    return P₂P₁*( ( (γ+1.0)/(γ-1.0) ) + P₂P₁ ) / ( 1.0 + ( (γ+1.0)/(γ-1.0) ) * P₂P₁ )
end

"""
    TemperatureRatioFromMachNumber( msw::MovingShockWave, γ::T, M::T ) where { T <: AbstractFloat }

Computes the temperature ratio T₂/T₁ across a moving shock wave, where state 1 is before the shock and state 2 is after the shock.

# Arguments
    msw: A MovingShockWave() type
    γ: The ratio of specific heats of the gas the shock wave is propagating into
    P₂P₁: The pressure ratio of the shock wave, must be greater than 1.0

# Returns
    A value of type T representing the temperature ratio across the shock wave, T₂/T₁

# Equation
    The temperature ratio is computed using Eq. 7.10 as:
          T₂/T₁ = P₂/P₁ ( ((γ+1)/(γ-1)) + P₂/P₁ ) / ( 1 + ((γ+1)/(γ-1))*P₂/P₁ )
    where
        ⋅ T₁ is the temperature before the shock wave
        ⋅ T₂ is the temperature after the shock wave
        ⋅ γ  is the ratio of specific heats for the gas
        ⋅ P₂/P₁ is the pressure ratio across the shock wave

# Note
    This function is a convenience wrapper around TemperatureRatio( ::MovingShockWave, γ, P₂P₁ ) to compute the temperature ratio in terms of a Mach number instead of a pressure ratio.
    The function PressureRatio( NormalShockWave(), γ, M ) is called to determine the pressure ratio from the Mach number, and this result is passed to TemperatureRatio( msw, γ, P₂P₁ )
"""
function TemperatureRatioFromMachNumber( msw::MovingShockWave, γ::T, M::T ) where { T <: AbstractFloat }
    P₂P₁ = PressureRatio( NormalShockWave(), γ, M )
    return TemperatureRatio( msw, γ, P₂P₁ )
end

"""
    ShockWaveVelocity( msw::MovingShockWave, γ::T, P₁::T, ρ₁::T, P₂P₁::T ) where { T <: AbstractFloat }

Computes the velocity of a moving shock wave

# Arguments
    msw: A MovingShockWave() type
    a₁: The speed of sound of the gas the shock is propagating into
    γ: The ratio of specific heats of the gas the shock wave is propagating into
    P₂P₁: The pressure ratio across the shock wave, must be greater than 1.0

# Returns
    A value of type T representing the velocity of the shock wave

# Equation
    The equation for the velocity of a moving shock wave is given by Eq. 7.14 as
          W = a₁ √( ( (γ+1)/2γ ) ( P₂/P₁ - 1 ) + 1 )
    where
        ⋅ a₁ = √(γ P₁/ρ₁) is the speed of sound of the gas the shock is propagating in to
        ⋅ γ is the ratio of specific heats of the gas the shock is propagating in to
        ⋅ P₁ is the pressure before the shock wave
        ⋅ P₂ is the pressure after the shock wave
"""
function ShockWaveVelocity( msw::MovingShockWave, a₁::T, γ::T, P₂P₁::T ) where { T <: AbstractFloat }
    return a₁ * sqrt( (γ + 1)/(2.0 * γ) * (P₂P₁ - 1.0) + 1.0 )
end

"""
    ShockWaveVelocity( msw::MovingShockWave, γ::T, P₁::T, ρ₁::T, P₂P₁::T ) where { T <: AbstractFloat }

Computes the velocity of a moving shock wave

# Arguments
    msw: A MovingShockWave() type
    γ: The ratio of specific heats of the gas the shock wave is propagating into
    P₁: The pressure in front of the shock wave
    ρ₁: The density in front of the shock wave
    P₂P₁: The pressure ratio across the shock wave, must be greater than 1.0

# Returns
    A value of type T representing the velocity of the shock wave

# Equation
    The equation for the velocity of a moving shock wave is given by Eq. 7.14 as
          W = a₁ √( ( (γ+1)/2γ ) ( P₂/P₁ - 1 ) + 1 )
    where
        ⋅ a₁ = √(γ P₁/ρ₁) is the speed of sound of the gas the shock is propagating in to
        ⋅ γ is the ratio of specific heats of the gas the shock is propagating in to
        ⋅ P₁ is the pressure before the shock wave
        ⋅ P₂ is the pressure after the shock wave

# Note
    This is intended as a convenience function to supply the upstream gas properties assuming the speed of sound is not already known.
    If the speed of sound of the upstream gas is already known, consider using the overloaded version of this function, ShockWaveVelocity( ::MovingShockWave, a₁, γ, P₂P₁ ) as it is what this function ultimately calls
"""
function ShockWaveVelocity( msw::MovingShockWave, γ::T, P₁::T, ρ₁::T, P₂P₁::T ) where { T <: AbstractFloat }
    a₁ = SpeedOfSound( γ, P₁, ρ₁ )
    return ShockWaveVelocity( msw, a₁, γ, P₂P₁ )
end

"""
    ShockWaveVelocity( msw::MovingShockWave, γ::T, P₁::T, ρ₁::T, P₂P₁::T ) where { T <: AbstractFloat }

Computes the velocity of a moving shock wave

# Arguments
    msw: A MovingShockWave() type
    γ: The ratio of specific heats of the gas the shock wave is propagating into
    P₁: The pressure in front of the shock wave
    ρ₁: The density in front of the shock wave
    M: The Mach number of the shock wave, must be greater than 1.0

# Returns
    A value of type T representing the velocity of the shock wave

# Equation
    The equation for the velocity of a moving shock wave is given by Eq. 7.14 as
          W = a₁ √( ( (γ+1)/2γ ) ( P₂/P₁ - 1 ) + 1 )
    where
        ⋅ a₁ = √(γ P₁/ρ₁) is the speed of sound of the gas the shock is propagating in to
        ⋅ γ is the ratio of specific heats of the gas the shock is propagating in to
        ⋅ P₁ is the pressure before the shock wave
        ⋅ P₂ is the pressure after the shock wave

# Note
    This function is a convenience wrapper around ShockWaveVelocity( ::MovingShockWave, γ, P₂P₁ ) to compute the shock wave velocity in terms of a Mach number instead of a pressure ratio.
    The function PressureRatio( NormalShockWave(), γ, M ) is called to determine the pressure ratio from the Mach number, and this result is passed to ShockWaveVelocity( msw, a₁, γ, P₂P₁ )
"""
function ShockWaveVelocityFromMachNumber( msw::MovingShockWave, a₁::T, γ::T, M::T ) where { T <: AbstractFloat }
    P₂P₁ = PressureRatio( NormalShockWave(), γ, M )
    return ShockWaveVelocity( msw, a₁, γ, P₂P₁ )
end

"""
    ShockWaveVelocity( msw::MovingShockWave, γ::T, P₁::T, ρ₁::T, P₂P₁::T ) where { T <: AbstractFloat }

Computes the velocity of a moving shock wave

# Arguments
    msw: A MovingShockWave() type
    γ: The ratio of specific heats of the gas the shock wave is propagating into
    P₁: The pressure in front of the shock wave
    ρ₁: The density in front of the shock wave
    M: The Mach number of the shock wave, must be greater than 1.0

# Returns
    A value of type T representing the velocity of the shock wave

# Equation
    The equation for the velocity of a moving shock wave is given by Eq. 7.14 as
          W = a₁ √( ( (γ+1)/2γ ) ( P₂/P₁ - 1 ) + 1 )
    where
        ⋅ a₁ = √(γ P₁/ρ₁) is the speed of sound of the gas the shock is propagating in to
        ⋅ γ is the ratio of specific heats of the gas the shock is propagating in to
        ⋅ P₁ is the pressure before the shock wave
        ⋅ P₂ is the pressure after the shock wave

# Note
    This function is a convenience wrapper around ShockWaveVelocity( ::MovingShockWave, γ, P₂P₁ ) to compute the shock wave velocity in terms of a Mach number instead of a pressure ratio.
    This function computes the speed of sound of the upstream gas, a₁, from the supplied values of γ, P₁, and ρ₁, and passes these to ShockWaveVelocityFromMachNumber( msw, a₁, γ, M )
"""
function ShockWaveVelocityFromMachNumber( msw::MovingShockWave, γ::T, P₁::T, ρ₁::T, M::T ) where { T <: AbstractFloat }
    a₁ = SpeedOfSound( γ, P₁, ρ₁ )
    return ShockWaveVelocityFromMachNumber( msw, a₁, γ, M )
end

"""
    MachNumberFromPressureRatio( msw::MovingShockWave, γ::T, P₂P₁::T ) where { T <: AbstractFloat }

Computes the Mach number of a moving shock wave based on the pressure ratio across the shock wave

# Arguments
    msw: A MovingShockWave() type
    γ: The ratio of specific heats of the gas the shock wave is propagating into
    P₂P₁: The pressure ratio across the shock wave, must be greater than 1.0

# Returns
    A value of type T representing the Mach number of the shock wave

# Note
    This is simply a convenience function around MachNumberFromPressureRatio( ::NormalShockWave, γ, P₂P₁ )
    For more information on how this function works, see the documentation for that function.
"""
function MachNumberFromPressureRatio( msw::MovingShockWave, γ::T, P₂P₁::T ) where { T <: AbstractFloat }
    return MachNumberFromPressureRatio( NormalShockWave(), γ, P₂P₁ )
end

"""
    PostShockVelocity( msw::MovingShockWave, γ::T, P₁::T, ρ₁::T, P₂P₁::T ) where { T <: AbstractFloat }

Computes the velocity of the fluid following the passage of a moving shock wave

# Arguments
    msw: A MovingShockWave() type
    γ: The ratio of specific heats of the gas the shock wave is propagating into
    P₁: The pressure in front of the shock wave
    ρ₁: The density in front of the shock wave
    P₂P₁: The pressure ratio across the shock wave, must be greater than 1.0

# Returns
    A value of type T representing the Mach number of the shock wave

# Equation
    The equation for the post-shock velocity is given by Eq. 7.15 as
          uₚ = W ( 1 - ρ₁/ρ₂ )
    where
        ⋅ uₚ is the post-shock velocity
        ⋅ W  is the velocity of the shock wave
        ⋅ ρ₁ is the density before the shock wave
        ⋅ ρ₂ is the density after the shock wave

# Note
    As an implementation detail, this function uses the supplied parameters to call DensityRatio( msw, ... ) and ShockWaveVelocity, and uses these results to compute uₚ.
    Anderson shows in Eq. 7.16 that the analytic expression for ρ₁/ρ₂ (Eq. 7.11) and wave velocity (Eq. 7.14) can be used to derive an analytic expression if this approach is preferred
"""
function PostShockVelocity( msw::MovingShockWave, γ::T, P₁::T, ρ₁::T, P₂P₁::T ) where { T <: AbstractFloat }
    ρ₂ρ₁ = DensityRatio( msw, γ, P₂P₁ )
    W = ShockWaveVelocity( msw, γ, P₁, ρ₁, P₂P₁ )
    return W * ( 1.0 - 1.0/ρ₂ρ₁ )
end

"""
    PostShockVelocity( msw::MovingShockWave, γ::T, P₁::T, ρ₁::T, P₂P₁::T ) where { T <: AbstractFloat }

Computes the velocity of the fluid following the passage of a moving shock wave

# Arguments
    msw: A MovingShockWave() type
    γ: The ratio of specific heats of the gas the shock wave is propagating into
    P₁: The pressure in front of the shock wave
    ρ₁: The density in front of the shock wave
    P₂P₁: The pressure ratio across the shock wave, must be greater than 1.0

# Returns
    A value of type T representing the Mach number of the shock wave

# Equation
    The equation for the post-shock velocity is given by Eq. 7.15 as
          uₚ = W ( 1 - ρ₁/ρ₂ )
    where
        ⋅ uₚ is the post-shock velocity
        ⋅ W  is the velocity of the shock wave
        ⋅ ρ₁ is the density before the shock wave
        ⋅ ρ₂ is the density after the shock wave

# Note
    As an implementation detail, this function uses the supplied parameters to call DensityRatio( msw, ... ) and ShockWaveVelocity, and uses these results to compute uₚ.
    Anderson shows in Eq. 7.16 that the analytic expression for ρ₁/ρ₂ (Eq. 7.11) and wave velocity (Eq. 7.14) can be used to derive an analytic expression if this approach is preferred

    This function uses the supplied Mach number to compute a pressure ratio from the normal shock relationships (PressureRatio( NormalShockWave(), ... )) and passes that value to the PostShockVelocity() function
"""
function PostShockVelocityFromMachNumber( msw::MovingShockWave, γ::T, P₁::T, ρ₁::T, M::T ) where { T <: AbstractFloat }
    P₂P₁ = PressureRatio( NormalShockWave(), γ, M )
    return PostShockVelocity( msw, γ, P₁, ρ₁, P₂P₁ )
end

export MovingShockWave
export DensityRatio, DensityRatioFromMachNumber
export TemperatureRatio, TemperatureRatioFromMachNumber
export ShockWaveVelocity, ShockWaveVelocityFromMachNumber
export MachNumberFromPressureRatio
export PostShockVelocity, PostShockVelocityFromMachNumber
