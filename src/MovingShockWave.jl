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

"""
    DiaphragmPressureRatioFromShockPressureRatio( msw::MovingShockWave, γ₁::T, P₁::T, ρ₁::T, γ₄::T, P₄::T, ρ₄::T, P₂P₁::T ) where { T <: AbstractFloat }

Returns the diaphragm pressure ratio required to form a shock wave with a desired shock wave pressure ratio given the driver and driven state

# Arguments
    γ: The ratio of specific heats of the gas
    P: The pressure
    ρ: The density
    P₂P₁: The desired shock wave pressure ratio

    Subscripts on γ, P, and ρ indicate whether the variable refers to the driven (1) or driver (4) sections

# Returns
    A value of type T representing the diaphragm pressure ratio required to create a shock with the desired pressure ratio

# Equation
    This function computes the right hand side of Eq. 7.94, which is
        P₄/P₁ = P₂/P₁ ( 1 - ( (γ₄-1) * (a₁/a₄) * (P₂/P₁-1) ) / √( 2γ₁ * ( 2γ₁ + (γ₁+1) (P₂/P₁-1) ) ) )^(-2γ₄/(γ₄-1))
    where
        ⋅ P is the pressure
        ⋅ γ is the ratio of specific heats
        ⋅ a = sqrt(γ*P/ρ) is the speed of sound

        and the subscripts refer to the driven section (1), behind the shock wave (2), or the driver section (4)
"""
function DiaphragmPressureRatioFromShockPressureRatio( msw::MovingShockWave, γ₁::T, P₁::T, ρ₁::T, γ₄::T, P₄::T, ρ₄::T, P₂P₁::T ) where { T <: AbstractFloat }
    a₁ = SpeedOfSound( γ₁, P₁, ρ₁ ) # Speed of sound in the driven section
    a₄ = SpeedOfSound( γ₄, P₄, ρ₄ ) # Speed of sound in the driver section
    a₁a₄ = a₁ / a₄ # Ratio of speed of sound
    return P₂P₁ * ( 1.0 - ( ( γ₄ - 1.0 ) * ( a₁a₄ ) * ( P₂P₁ - 1.0 ) ) / sqrt( 2.0 * γ₁ * ( 2.0 * γ₁ + ( γ₁ + 1.0 ) * ( P₂P₁ - 1.0 ) ) ) )^(-2.0 * γ₄ / ( γ₄-1.0 ) )
end

"""
    DiaphragmPressureRatioFromMachNumber( msw::MovingShockWave, γ₁::T, P₁::T, ρ₁::T, γ₄::T, P₄::T, ρ₄::T, M::T ) where { T <: AbstractFloat }

Returns the diaphragm pressure ratio required to form a shock wave with a desired Mach number given the driver and driven state

# Arguments
    γ: The ratio of specific heats of the gas
    P: The pressure
    ρ: The density
    M: The desired shock wave Mach number

    Subscripts on γ, P, and ρ indicate whether the variable refers to the driven (1) or driver (4) sections

# Returns
    A value of type T representing the diaphragm pressure ratio required to create a shock with the desired pressure ratio

# Notes
    This function is a wrapper around DiaphragmPressureRatioFromShockPressureRatio(). It computes the pressure ratio of the shock, P₂/P₁, from the supplied Mach number using the PressureRatio() function for normal shock waves and passes the result to DiaphragmPressureRatioFromShockPressureRatio()
"""
function DiaphragmPressureRatioFromMachNumber( msw::MovingShockWave, γ₁::T, P₁::T, ρ₁::T, γ₄::T, P₄::T, ρ₄::T, M::T ) where { T <: AbstractFloat }
    P₂P₁ = PressureRatio( NormalShockWave(), γ₁, M )
    return DiaphragmPressureRatioFromShockPressureRatio( msw, γ₁, P₁, ρ₁, γ₄, P₄, ρ₄, P₂P₁ )
end

"""
    ShockPressureRatioFromDiaphragmPressureRatio( msw::MovingShockWave, γ₁::T, P₁::T, ρ₁::T, γ₄::T, P₄::T, ρ₄::T; P₂P₁::T=2.0, ΔP::T=1e-3 ) where { T <: AbstractFloat }

Returns the shock wave pressure ratio that is generated from the rupture of a diaphragm with gases on either side at a supplied state

# Arguments
    γ: The ratio of specific heats of the gas
    P: The pressure
    ρ: The density

    Subscripts on γ, P, and ρ indicate whether the variable refers to the driven (1) or driver (4) sections

# Optional Arguments
    P₂P₁: An initial guess for the shock wave pressure ratio. Default is 2.0
    ΔP: The step size used to compute the numerical derivative. See Notes below for how this parameter is used

# Returns
    A value of type T representing the shock wave pressure ratio generated by the rupture of a diaphragm with a given pressure ratio

# Equation
    This function finds the solution of Eq. 7.94, as implemented by DiaphragmPressureRatioFromShockPressureRatio(). See the documentation of that function for more details

# Notes
    As an implementation detail, the solution to Eq. 7.94 is found using an iterative Newton-Raphson method with an initial guess of P₂P₁.
    In order to compute the derivative, two function evaluations are performed at each iteration - the first at the current guess for P₂P₁, and the second at P₂P₁ + ΔP, where ΔP is a small increment.
"""
function ShockPressureRatioFromDiaphragmPressureRatio( msw::MovingShockWave, γ₁::T, P₁::T, ρ₁::T, γ₄::T, P₄::T, ρ₄::T; P₂P₁::T=2.0, ΔP::T=1e-3 ) where { T <: AbstractFloat }
    P₄P₁ = P₄ / P₁ # Pressure ratio across the diaphragm

    # Now need to find the root of the equation
    # We'll use a bisection method for this as it is fairly reliable
    ϵ = 100 # Error
    ϵₘ = 1e-6 # Maximum error
    i = 0 # Iteration number
    iₘ = 100 # Maximum number of iterations

    while ( ϵ > ϵₘ )
        if ( i >= iₘ )
            @warn "Exceeded maximum number of iterations. Returning current estimate of P₂/P₁=$P₂P₁ (ϵ=$ϵ)."
            break
        end
        # We'll compute the function at P/₂P₁, and at P₂/P₁ + ΔP, to determine a derivative numerically
        f₀ = DiaphragmPressureRatioFromShockPressureRatio( msw, γ₁, P₁, ρ₁, γ₄, P₄, ρ₄, P₂P₁ ) - P₄P₁ # Evaluation of our function at P₂P₁ 
        f₁ = DiaphragmPressureRatioFromShockPressureRatio( msw, γ₁, P₁, ρ₁, γ₄, P₄, ρ₄, P₂P₁ + ΔP ) - P₄P₁ # Evaluation of our function at P₂P₁ + ΔP
        fₚ = ( f₁ - f₀ ) / ΔP # Evaluation of the derivative
        P₂P₁ = P₂P₁ - f₀ / fₚ # Compute our next guess for pressure
        ϵ = f₀ # Reset our error estimate to the value of the function evaluation
    end

    return P₂P₁
end

"""
    MachNumberFromDiaphragmPressureRatio( msw::MovingShockWave, γ₁::T, P₁::T, ρ₁::T, γ₄::T, P₄::T, ρ₄::T; P₂P₁::T=2.0, ΔP::T=1e-3 ) where { T <: AbstractFloat }

Returns the Mach number of the shock wave that is generated from the rupture of a diaphragm with gases on either side at a supplied state

# Arguments
    γ: The ratio of specific heats of the gas
    P: The pressure
    ρ: The density

    Subscripts on γ, P, and ρ indicate whether the variable refers to the driven (1) or driver (4) sections

# Optional Arguments
    P₂P₁: An initial guess for the shock wave pressure ratio. Default is 2.0
    ΔP: The step size used to compute the numerical derivative. See Notes below for how this parameter is used

# Returns
    A value of type T representing the shock wave pressure ratio generated by the rupture of a diaphragm with a given pressure ratio

# Equation
    This function finds the solution of Eq. 7.94, as implemented by DiaphragmPressureRatioFromShockPressureRatio(). See the documentation of that function for more details.
    Additionally, this function is simply a wrapper around ShockPressureRatioFromDiaphragmPressureRatio(), and computes a Mach number from the pressure ratio returned from that function. For more information, see the documentation of that function.
"""
function MachNumberFromDiaphragmPressureRatio( msw::MovingShockWave, γ₁::T, P₁::T, ρ₁::T, γ₄::T, P₄::T, ρ₄::T; P₂P₁::T=2.0, ΔP::T=1e-3 ) where { T <: AbstractFloat }
    P₂P₁ = ShockPressureRatioFromDiaphragmPressureRatio( msw, γ₁, P₁, ρ₁, γ₄, P₄, ρ₄; P₂P₁, ΔP )
    MachNumberFromPressureRatio( msw, γ₁, P₂P₁ )
end

export MovingShockWave
export DensityRatio, DensityRatioFromMachNumber
export TemperatureRatio, TemperatureRatioFromMachNumber
export ShockWaveVelocity, ShockWaveVelocityFromMachNumber
export MachNumberFromPressureRatio
export PostShockVelocity, PostShockVelocityFromMachNumber
export DiaphragmPressureRatioFromShockPressureRatio, DiaphragmPressureRatioFromMachNumber 
export ShockPressureRatioFromDiaphragmPressureRatio, MachNumberFromDiaphragmPressureRatio
