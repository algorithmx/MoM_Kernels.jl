# =============================================================================
# Port Impedance Measurement
# =============================================================================
#
# This file provides functions to measure port currents and compute port 
# impedance for given excitation.
#
# Uses multiple dispatch to handle different port types.
#
# Note: extractPortVoltage functions have been removed as they computed
# (Z·I_j)[rwgID] which gives zero for off-diagonal entries in multi-port
# analysis. The computeSParameters function now uses the Y-matrix inversion
# approach which correctly computes all Z_ij terms.
#
# =============================================================================

"""
    extractPortCurrent(port::DeltaGapPort{FT,IT}, I_coeff::Vector{CT}) where {CT,FT,IT}

Extract port current for a DeltaGapPort (single-edge voltage source).

Applies the same length factor as `excitePort!` (power-balance duality,
A0_EDGE_PORT.md §6.4):

    I_port = f × I_coeff[rwgID]

where f is the dual of the excitation factor:
  - Full RWG (triID_neg > 0, two adjacent triangles): f = edgel / 2
  - Half/boundary RWG (triID_neg == 0, one triangle):  f = edgel

This matches the `length_factor` used in `excitePort!(V, port::DeltaGapPort)`.
"""
function extractPortCurrent(port::DeltaGapPort{FT,IT}, I_coeff::Vector{CT}) where {CT,FT,IT}
    rwgID = port.rwgID

    if rwgID <= 0 || rwgID > length(I_coeff)
        error("Invalid port rwgID: $(rwgID)")
    end

    # Dual of excitePort! length_factor: l/2 for full RWG, l for half (boundary) RWG.
    length_factor = port.triID_neg > 0 ? port.edgel / 2 : port.edgel

    return length_factor * I_coeff[rwgID]
end

"""
    extractPortCurrent(port::CurrentProbe{FT,IT}, I_coeff::Vector{CT}) where {CT,FT,IT}

Extract port current for a CurrentProbe (single-edge current source).

For a CurrentProbe the excitation entry is `V[rwgID] = port.I` with no
length factor (the testing integral collapses to \$V_m = I_0\$ directly,
see CurrentProbe docstring). By power-balance duality the length factor
is therefore 1, so the raw RWG coefficient is the port current:

    I_port = I_coeff[rwgID]
"""
function extractPortCurrent(port::CurrentProbe{FT,IT}, I_coeff::Vector{CT}) where {CT,FT,IT}
    rwgID = port.rwgID
    
    if rwgID <= 0 || rwgID > length(I_coeff)
        error("Invalid port rwgID: $(rwgID)")
    end
    
    return I_coeff[rwgID]
end

"""
    extractPortCurrent(port::DeltaGapArrayPort{FT,IT,DT}, I_coeff::Vector{CT}) where {CT,FT,IT,DT}

Extract total port current for DeltaGapArrayPort.

The length factor per edge must be the *dual* of the factor used in `excitePort!`
so that the power balance P = ½ V_port* I_port holds:

    P = ½ Σ_k V_k* I_k = ½ V_port* Σ_k w_k f_k I_k

The factor f_k mirrors the excitation convention (A0_EDGE_PORT.md §5.1–5.2):
  - Full RWG (triID_neg > 0, two adjacent triangles): f_k = l_k / 2
  - Half/boundary RWG (triID_neg == 0, one triangle):  f_k = l_k

This matches exactly the length_factor used in `excitePort!` and `_excitation_array!`.

    I_total = Σ_edge  weight_edge × I_coeff[rwgID] × f_k
"""
function extractPortCurrent(port::DeltaGapArrayPort{FT,IT,DT}, I_coeff::Vector{CT}) where {CT<:Complex, FT<:Real, IT<:Integer, DT}
    if !port.isBound
        error("Port must be bound to mesh")
    end
    
    I_total = zero(CT)
    
    for i in eachindex(port.rwgIDs)
        rwgID = port.rwgIDs[i]
        weight = port.edgeWeights[i]
        
        if rwgID <= 0 || rwgID > length(I_coeff)
            continue
        end
        
        # Dual of excitePort! length_factor: l/2 for full RWG, l for half (boundary) RWG.
        length_factor = port.triID_neg[i] > 0 ?
                        port.edgeLengths[i] / 2 :
                        port.edgeLengths[i]
        
        I_total += weight * I_coeff[rwgID] * length_factor
    end
    
    return I_total
end

"""
    extractPortCurrent(port::RectangularEdgePort{FT,IT,DT}, I_coeff::Vector{CT}) where {CT,FT,IT,DT}

Extract total port current for a RectangularEdgePort.

Delegates to the underlying `DeltaGapArrayPort` base via composition.
"""
function extractPortCurrent(port::RectangularEdgePort{FT,IT,DT}, I_coeff::Vector{CT}) where {CT<:Complex, FT<:Real, IT<:Integer, DT}
    return extractPortCurrent(getfield(port, :base), I_coeff)
end

"""
    measurePortImpedance(port_i::PortType, Z, I_j::Vector{CT}, port_j::PortType) where CT

Measure port impedance Z_ij = V_i / I_j.

Computes the impedance seen at port i when port j is excited.

# Arguments
- `port_i::PortType`: Measurement port (port i)
- `Z`: Impedance matrix or operator (unused in current implementations)
- `I_j::Vector`: Solution current vector from exciting port j
- `port_j::PortType`: Excitation port (port j)

# Returns
- `Complex{FT}`: Port impedance Z_ij

# Note on Multi-Port Analysis
For multi-port Z_ij (i≠j), the correct approach uses the Y-matrix method
implemented in `computeSParameters`. The single-port versions of this
function (DeltaGapPort, CurrentProbe) are not suitable for multi-port
analysis as they compute (Z·I_j)[rwgID_i] which equals zero for passive ports.

The array port versions (DeltaGapArrayPort, RectangularEdgePort) compute
driving-point impedance Z_drive = V_j^src / I_j^port, which is only valid
for single-port analysis (port_i == port_j).
"""
function measurePortImpedance(
    port_i::PortType,
    Z,
    I_j::Vector{CT},
    port_j::PortType
) where {CT<:Complex}
    error("measurePortImpedance not implemented for port types $(typeof(port_i)), $(typeof(port_j))")
end

# Note: measurePortImpedance for DeltaGapPort and CurrentProbe have been removed
# as they used an incorrect algorithm (computed (Z·I_j)[rwgID] which gives zero
# for off-diagonal entries). Use computeSParameters for correct multi-port analysis.

"""
    measurePortImpedance(port_i::DeltaGapArrayPort{FT,IT,DT}, Z, I_j::Vector{CT}, port_j::PortType) where {CT,FT,IT,DT}

Driving-point impedance for a DeltaGapArrayPort:

    Z_drive = V_j^src / I_j^port

where `V_j^src = port_j.V` is the source voltage that generated `I_j`, and
`I_j^port = extractPortCurrent(port_j, I_j)` is the resulting port current.

# Single-port vs multi-port
This formula is exact for single-port analysis (`port_i == port_j`).
For multi-port Z_ij (i≠j) the correct value requires inverting the full
short-circuit admittance matrix; use `computeSParameters` for that.
"""
function measurePortImpedance(
    port_i::DeltaGapArrayPort{FT,IT,DT},
    Z,
    I_j::Vector{CT},
    port_j::PortType
) where {CT<:Complex, FT<:Real, IT<:Integer, DT}
    I_exc = extractPortCurrent(port_j, I_j)
    
    if abs(I_exc) < eps(FT)
        return Complex{FT}(Inf, 0)
    end
    
    return port_j.V / I_exc
end

"""
    measurePortImpedance(port_i::RectangularEdgePort{FT,IT,DT}, Z, I_j::Vector{CT}, port_j::PortType) where {CT,FT,IT,DT}

Driving-point impedance for a RectangularEdgePort. Delegates to the
underlying `DeltaGapArrayPort` base for current extraction.

    Z_drive = V_j^src / I_j^port

See `measurePortImpedance(::DeltaGapArrayPort, ...)` for the single-port vs
multi-port caveat.
"""
function measurePortImpedance(
    port_i::RectangularEdgePort{FT,IT,DT},
    Z,
    I_j::Vector{CT},
    port_j::PortType
) where {CT<:Complex, FT<:Real, IT<:Integer, DT}
    I_exc = extractPortCurrent(port_j, I_j)

    if abs(I_exc) < eps(FT)
        return Complex{FT}(Inf, 0)
    end

    return port_j.V / I_exc
end

# =============================================================================
# Exports
# =============================================================================

export extractPortCurrent, measurePortImpedance
