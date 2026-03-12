# =============================================================================
# S-Parameter Types and Data Structures
# =============================================================================
#
# This file defines the data structures for S-parameter calculation results.
#
# =============================================================================

"""
    SParameterResult{FT<:Real}

Container for S-parameter calculation results at single or multiple frequencies.

# Fields
- `S::Array{Complex{FT},3}`: S-parameter matrix (N×N×F where F = number of frequencies)
- `Z_port::Array{Complex{FT},3}`: Port impedance matrix (N×N×F)
- `frequencies::Vector{FT}`: Frequency points in Hz
- `Z0::FT`: Reference impedance in Ohms
- `portIDs::Vector{Int}`: Port identifier numbers

# Access Patterns
```julia
# Single frequency result
sparams.S[:,:,1]          # Full S-matrix at first frequency
sparams.S[1,1,:]          # S11 across all frequencies
sparams.S[2,1,5]          # S21 at 5th frequency point
```
"""
struct SParameterResult{FT<:Real}
    S::Array{Complex{FT},3}        # N×N×F array
    Z_port::Array{Complex{FT},3}   # N×N×F array
    frequencies::Vector{FT}
    Z0::FT
    portIDs::Vector{Int}
    
    function SParameterResult(
        S::AbstractArray{Complex{FT}},
        Z_port::AbstractArray{Complex{FT}},
        frequencies::Vector{FT},
        Z0::FT,
        portIDs::Vector{Int}
    ) where {FT<:Real}
        # Validate dimensions
        if size(S, 3) != length(frequencies)
            error("Frequency dimension mismatch: S has $(size(S, 3)) frequency points, " *
                  "but frequencies vector has $(length(frequencies)) points")
        end
        if size(S, 1) != size(S, 2)
            error("S-parameter matrix must be square, got $(size(S, 1))×$(size(S, 2))")
        end
        if size(S, 1) != length(portIDs)
            error("Port dimension mismatch: S-matrix is $(size(S, 1))×$(size(S, 1)), " *
                  "but there are $(length(portIDs)) port IDs")
        end
        if size(S) != size(Z_port)
            error("S and Z_port must have same dimensions: $(size(S)) vs $(size(Z_port))")
        end
        
        new{FT}(S, Z_port, frequencies, Z0, portIDs)
    end
end

"""
    SParameterResult(S::Matrix{CT}, Z_port::Matrix{CT}, freq::FT, Z0::FT, portIDs::Vector{Int}) where {CT,FT}

Convenience constructor for single-frequency results.
"""
function SParameterResult(
    S::Matrix{CT},
    Z_port::Matrix{CT},
    freq::FT,
    Z0::FT,
    portIDs::Vector{Int}
) where {CT<:Complex, FT<:Real}
    # Reshape to 3D arrays with single frequency dimension
    S_3d = reshape(S, size(S, 1), size(S, 2), 1)
    Z_port_3d = reshape(Z_port, size(Z_port, 1), size(Z_port, 2), 1)
    return SParameterResult(S_3d, Z_port_3d, [freq], Z0, portIDs)
end

"""
    Base.getproperty(result::SParameterResult, sym::Symbol)

Convenience accessors for common S-parameters.
"""
function Base.getproperty(result::SParameterResult{FT}, sym::Symbol) where FT
    if sym == :num_ports
        return length(result.portIDs)
    elseif sym == :num_frequencies
        return length(result.frequencies)
    else
        return getfield(result, sym)
    end
end

"""
    SParameterSweepConfig{FT<:Real}

Configuration for frequency sweep S-parameter analysis.

# Fields
- `frequencies::Vector{FT}`: Frequency points
- `Z0::FT`: Reference impedance
- `solverT::Symbol`: Default solver type
- `solverKwargs::Dict{Symbol,Any}`: Default solver parameters
"""
struct SParameterSweepConfig{FT<:Real}
    frequencies::Vector{FT}
    Z0::FT
    solverT::Symbol
    solverKwargs::Dict{Symbol,Any}
    
    function SParameterSweepConfig(
        frequencies::Vector{FT},
        Z0::FT = FT(50.0),
        solverT::Symbol = :direct;
        kwargs...
    ) where {FT<:Real}
        return new{FT}(frequencies, Z0, solverT, Dict(kwargs))
    end
end

# =============================================================================
# Exports
# =============================================================================

export SParameterResult, SParameterSweepConfig
