# =============================================================================
# Port Excitation Vector Construction
# =============================================================================
#
# This file provides functions to build excitation vectors for specific ports
# in a PortArray. Used by computeSParameters to excite each port sequentially.
#
# =============================================================================

"""
    excitePort!(V::Vector{CT}, port::PortType) where CT

Set excitation vector for a single port.

For voltage ports: V[port.rwgID] = V_port × length_factor
For current ports: V[port.rwgID] = I_port

# Arguments
- `V::Vector{Complex{FT}}`: Pre-allocated excitation vector (modified in-place)
- `port::PortType`: Port to excite

# Returns
- `V`: The modified excitation vector
"""
function excitePort!(V::Vector{CT}, port::PortType) where {CT<:Complex}
    error("excitePort! not implemented for port type $(typeof(port))")
end

"""
    excitePort!(V::Vector{CT}, port::DeltaGapPort{FT,IT}) where {CT,FT,IT}

Excite a DeltaGapPort (voltage source).

The excitation is: V[rwgID] = V_port × (edge_length / 2) for full RWG
                    V[rwgID] = V_port × edge_length for half RWG (boundary)
"""
function excitePort!(V::Vector{CT}, port::DeltaGapPort{FT,IT}) where {CT<:Complex,FT,IT}
    rwgID = port.rwgID
    
    if rwgID <= 0 || rwgID > length(V)
        error("Invalid port rwgID: $(rwgID), vector length: $(length(V))")
    end
    
    # Full RWG: V × l/2, Half RWG (boundary): V × l
    length_factor = port.triID_neg > 0 ? port.edgel / 2 : port.edgel
    
    V[rwgID] = port.V * length_factor
    
    return V
end

"""
    excitePort!(V::Vector{CT}, port::CurrentProbe{FT,IT}) where {CT,FT,IT}

Excite a CurrentProbe (current source).

The excitation is: V[rwgID] = I_probe
"""
function excitePort!(V::Vector{CT}, port::CurrentProbe{FT,IT}) where {CT<:Complex,FT,IT}
    rwgID = port.rwgID
    
    if rwgID <= 0 || rwgID > length(V)
        error("Invalid port rwgID: $(rwgID), vector length: $(length(V))")
    end
    
    # Direct current injection
    V[rwgID] = port.I
    
    return V
end

"""
    excitePort!(V::Vector{CT}, port::DeltaGapArrayPort{FT,IT,DT}) where {CT,FT,IT,DT}

Excite a DeltaGapArrayPort (multi-edge voltage port).

Applies excitation to all boundary edges with their respective weights:
V[rwgID] += V_port × weight × length_factor

Requires the port to be bound to mesh.
"""
function excitePort!(V::Vector{CT}, port::DeltaGapArrayPort{FT,IT,DT}) where {CT,FT,IT,DT}
    if !port.isBound
        error("DeltaGapArrayPort must be bound to mesh before excitation. " *
              "Call bind_to_mesh!() first.")
    end
    
    for i in eachindex(port.rwgIDs)
        rwgID = port.rwgIDs[i]
        
        if rwgID <= 0 || rwgID > length(V)
            continue  # Skip invalid indices
        end
        
        # Weight for this edge
        weight = port.edgeWeights[i]
        
        # Length factor: l/2 for full RWG, l for half RWG
        length_factor = port.triID_neg[i] > 0 ? 
                        port.edgeLengths[i] / 2 : 
                        port.edgeLengths[i]
        
        # Apply excitation
        V[rwgID] += port.V * weight * length_factor
    end
    
    return V
end

"""
    excitePort!(V::Vector{CT}, port::RectangularEdgePort{FT,IT,DT}) where {CT,FT,IT,DT}

Excite a RectangularEdgePort.

Delegates to DeltaGapArrayPort implementation.
"""
function excitePort!(V::Vector{CT}, port::RectangularEdgePort{FT,IT,DT}) where {CT,FT,IT,DT}
    # Convert to DeltaGapArrayPort and excite
    # This is a view conversion - no data is copied
    gap_port = _to_delta_gap_array_port(port)
    return excitePort!(V, gap_port)
end

"""
    buildExcitationVector(ports::PortArray{FT,IT,PT}, port_idx::Integer, nbf::Integer) where {FT,IT,PT}

Build excitation vector for a specific port in a PortArray.

This is the main entry point used by computeSParameters.

# Arguments
- `ports::PortArray`: Port array containing all ports
- `port_idx::Integer`: Index of port to excite (1-based)
- `nbf::Integer`: Total number of basis functions (vector length)

# Returns
- `Vector{Complex{FT}}`: Excitation vector with single port excited
"""
function buildExcitationVector(
    ports::PortArray{FT,IT,PT},
    port_idx::Integer,
    nbf::Integer
) where {FT,IT,PT}
    CT = Complex{FT}
    V = zeros(CT, nbf)
    
    if port_idx < 1 || port_idx > ports.numPorts
        error("Port index $port_idx out of range (1:$(ports.numPorts))")
    end
    
    port = ports.ports[port_idx]
    excitePort!(V, port)
    
    return V
end

"""
    buildExcitationVector(port::PT, nbf::Integer) where {PT<:PortType}

Build excitation vector for a single port.

Convenience method for single-port analysis.
"""
function buildExcitationVector(port::PT, nbf::Integer) where {PT<:PortType}
    # Infer element type from port (try to get from port's value type)
    V = zeros(ComplexF64, nbf)
    excitePort!(V, port)
    return V
end

# =============================================================================
# Exports
# =============================================================================

export excitePort!, buildExcitationVector
