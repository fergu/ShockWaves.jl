module ShockWaves

# An abstract type defining a shock wave
# Each of the files below will define a concrete type as a descendant of this type to specialize their functions
abstract type AbstractShockWave end

include("NormalShockWave.jl")      # Functions to compute properties across a normal shock wave

end
