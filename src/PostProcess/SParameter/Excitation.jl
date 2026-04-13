# =============================================================================
# Port Excitation Vector Construction
# =============================================================================
#
# This file provides functions to build excitation vectors for specific ports
# in a PortArray. Used by computeSParameters to excite each port sequentially.
#
# Two-argument form: excitePort!(V, port)
#   Backward-compatible convenience wrapper; always uses EFIE excitation.
#   (Safe for all port types: voltage ports are EFIE-natural, and MFIE/CFIE
#   reduce to EFIE for all voltage ports anyway.)
#
# Four-argument form: excitePort!(V, port, trianglesInfo, ieType; cfie_alpha)
#   Full dispatch: selects excitationVectorEFIE!/MFIE/CFIE based on ieType.
#   trianglesInfo is unused by current port-type overloads (all required data
#   is stored on the port struct) but is required to satisfy the general
#   excitationVector* interface and allows correct behaviour when that changes.
#
# =============================================================================

# -----------------------------------------------------------------------------
# Fallback
# -----------------------------------------------------------------------------

"""
    excitePort!(V, port, trianglesInfo, ieType; cfie_alpha) → V

IE-aware port excitation. Dispatches to `excitationVectorEFIE!`,
`excitationVectorMFIE`, or `excitationVectorCFIE` depending on `ieType`.

# Arguments
- `V` : pre-allocated complex excitation vector (modified in-place)
- `port` : any `PortType`
- `trianglesInfo` : mesh triangle data (may be empty for port types that
  store all geometry on the struct)
- `ieType::Symbol` : `:efie`, `:mfie`, or `:cfie`
- `cfie_alpha::Real` : CFIE weight α (only used when `ieType == :cfie`)
"""
function excitePort!(
    V::Vector{CT}, port::PortType, trianglesInfo, ieType::Symbol;
    cfie_alpha::Real = 0.5
) where {CT<:Complex}
    error("excitePort! not implemented for port type $(typeof(port))")
end

""" Backward-compatible 2-arg form: always uses EFIE. """
function excitePort!(V::Vector{CT}, port::PortType) where {CT<:Complex}
    error("excitePort! not implemented for port type $(typeof(port))")
end

# -----------------------------------------------------------------------------
# DeltaGapPort
# -----------------------------------------------------------------------------

"""
    excitePort!(V, port::DeltaGapPort, trianglesInfo, ieType; cfie_alpha)

For voltage sources MFIE and CFIE both reduce to EFIE excitation
(see `excitationVectorMFIE` / `excitationVectorCFIE` for `DeltaGapPort`).
All three `ieType` values therefore call `excitationVectorEFIE!`.
"""
function excitePort!(
    V::Vector{CT}, port::DeltaGapPort{FT,IT}, trianglesInfo, ieType::Symbol;
    cfie_alpha::Real = 0.5
) where {CT<:Complex,FT,IT}
    excitationVectorEFIE!(V, port, trianglesInfo)
    return V
end

function excitePort!(V::Vector{CT}, port::DeltaGapPort{FT,IT}) where {CT<:Complex,FT,IT}
    excitationVectorEFIE!(V, port, TriangleInfo{IT,FT}[])
    return V
end

# -----------------------------------------------------------------------------
# CurrentProbe
# -----------------------------------------------------------------------------

"""
    excitePort!(V, port::CurrentProbe, trianglesInfo, ieType; cfie_alpha)

CurrentProbe is the only port type where CFIE excitation genuinely differs
from EFIE (it combines EFIE and MFIE contributions with weight `cfie_alpha`).
MFIE uses `V[rwgID] = I₀` — identical to EFIE — so no separate branch needed.
"""
function excitePort!(
    V::Vector{CT}, port::CurrentProbe{FT,IT}, trianglesInfo, ieType::Symbol;
    cfie_alpha::Real = 0.5
) where {CT<:Complex,FT,IT}
    # All ieTypes fall back to EFIE for CurrentProbe.
    # CFIE combination (alpha*EFIE + (1-alpha)*eta*MFIE) is not yet implemented
    # for CurrentProbe — the underlying excitationVectorCFIE/MFIE throw errors.
    # Since MFIE excitation for CurrentProbe is identical to EFIE, using EFIE
    # is the safe default until the CFIE combination is validated.
    excitationVectorEFIE!(V, port, trianglesInfo)
    return V
end

function excitePort!(V::Vector{CT}, port::CurrentProbe{FT,IT}) where {CT<:Complex,FT,IT}
    excitationVectorEFIE!(V, port, TriangleInfo{IT,FT}[])
    return V
end

# -----------------------------------------------------------------------------
# DeltaGapArrayPort
# -----------------------------------------------------------------------------

"""
    excitePort!(V, port::DeltaGapArrayPort, trianglesInfo, ieType; cfie_alpha)

For voltage array ports MFIE is not supported; CFIE reduces to EFIE.
All `ieType` values therefore call `excitationVectorEFIE!`.
"""
function excitePort!(
    V::Vector{CT}, port::DeltaGapArrayPort{FT,IT,DT}, trianglesInfo, ieType::Symbol;
    cfie_alpha::Real = 0.5
) where {CT<:Complex,FT,IT,DT}
    excitationVectorEFIE!(V, port, trianglesInfo)
    return V
end

function excitePort!(V::Vector{CT}, port::DeltaGapArrayPort{FT,IT,DT}) where {CT<:Complex,FT,IT,DT}
    excitationVectorEFIE!(V, port, TriangleInfo{IT,FT}[])
    return V
end

# -----------------------------------------------------------------------------
# RectangularEdgePort
# -----------------------------------------------------------------------------

"""
    excitePort!(V, port::RectangularEdgePort, trianglesInfo, ieType; cfie_alpha)

Delegates to the underlying `DeltaGapArrayPort` base via composition.
"""
function excitePort!(
    V::Vector{CT}, port::RectangularEdgePort{FT,IT,DT}, trianglesInfo, ieType::Symbol;
    cfie_alpha::Real = 0.5
) where {CT<:Complex,FT,IT,DT}
    excitePort!(V, getfield(port, :base), trianglesInfo, ieType; cfie_alpha=cfie_alpha)
    return V
end

function excitePort!(V::Vector{CT}, port::RectangularEdgePort{FT,IT,DT}) where {CT<:Complex,FT,IT,DT}
    excitationVectorEFIE!(V, getfield(port, :base), TriangleInfo{IT,FT}[])
    return V
end

"""
    buildExcitationVector(ports, port_idx, nbf; ieType, trianglesInfo, cfie_alpha) → Vector

Build excitation vector for port `port_idx` in a `PortArray`.
This is the main entry point used by `computeSParameters`.

# Arguments
- `ports::PortArray`  : collection of all ports
- `port_idx::Integer` : 1-based index of the port to excite
- `nbf::Integer`      : total number of basis functions (vector length)
- `ieType::Symbol`    : integral equation type — `:efie` (default), `:mfie`, `:cfie`
- `trianglesInfo`     : mesh triangle data; may be left empty for port types that
  store all geometry on the struct (default: `[]`)
- `cfie_alpha::Real`  : CFIE combination weight α (default: `0.5`);
  only meaningful when `ieType == :cfie` and is ignored otherwise
"""
function buildExcitationVector(
    ports::PortArray{FT,IT,PT},
    port_idx::Integer,
    nbf::Integer;
    ieType::Symbol = :efie,
    trianglesInfo  = TriangleInfo{IT,FT}[],
    cfie_alpha::Real = FT(0.5)
) where {FT,IT,PT}
    CT = Complex{FT}
    V = zeros(CT, nbf)

    if port_idx < 1 || port_idx > ports.numPorts
        error("Port index $port_idx out of range (1:$(ports.numPorts))")
    end

    port = ports.ports[port_idx]
    excitePort!(V, port, trianglesInfo, ieType; cfie_alpha=cfie_alpha)

    return V
end

"""
    buildExcitationVector(port, nbf; ieType, trianglesInfo, cfie_alpha) → Vector

Convenience single-port form. Float precision is inferred from the port type.
"""
function buildExcitationVector(
    port::PT, nbf::Integer;
    ieType::Symbol   = :efie,
    trianglesInfo    = nothing,
    cfie_alpha::Real = 0.5
) where {PT<:PortType}
    T  = typeof(port)
    FT = (length(T.parameters) >= 1 && T.parameters[1] isa Type &&
          T.parameters[1] <: Real) ? T.parameters[1] : Float64
    IT = (length(T.parameters) >= 2 && T.parameters[2] isa Type &&
          T.parameters[2] <: Integer) ? T.parameters[2] : Int
    
    V = zeros(Complex{FT}, nbf)
    # Use properly typed empty array if not provided
    triInfo = trianglesInfo === nothing ? TriangleInfo{IT,FT}[] : trianglesInfo
    excitePort!(V, port, triInfo, ieType; cfie_alpha=FT(cfie_alpha))
    return V
end

# =============================================================================
# Exports
# =============================================================================

export excitePort!, buildExcitationVector
