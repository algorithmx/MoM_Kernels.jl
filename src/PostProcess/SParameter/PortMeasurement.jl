# =============================================================================
# Port Impedance Measurement
# =============================================================================
#
# This file provides functions to measure port voltages and currents,
# and compute port impedance Z_ij = V_i / I_j for given excitation.
#
# Uses multiple dispatch to handle different port types.
#
# =============================================================================

"""
    extractPortVoltage(port::PortType, Z, I_coeff::Vector{CT}) where CT

Extract port voltage from solution vector.

For voltage ports: Returns the source voltage (known)
For current ports: Computes V = Z[port.rwgID, :] ⋅ I_coeff
"""
function extractPortVoltage(port::PortType, Z, I_coeff::Vector{CT}) where {CT<:Complex}
    error("extractPortVoltage not implemented for port type $(typeof(port))")
end

"""
    extractPortVoltage(port::DeltaGapPort{FT,IT}, Z, I_coeff::Vector{CT}) where {CT,FT,IT}

Extract voltage for a DeltaGapPort (voltage source).

For voltage sources, we compute the actual voltage response from the solution,
not the source voltage. This gives the true port impedance.
"""
function extractPortVoltage(
    port::DeltaGapPort{FT,IT},
    Z,
    I_coeff::Vector{CT}
) where {CT<:Complex, FT<:Real, IT<:Integer}
    rwgID = port.rwgID
    
    if rwgID <= 0 || rwgID > length(I_coeff)
        error("Invalid port rwgID: $(rwgID)")
    end
    
    # Compute voltage from impedance matrix row
    # V_port = Σ_k Z[port, k] * I_coeff[k]
    V_port = zero(CT)
    
    # Handle both dense matrices and MLFMA operators
    if Z isa AbstractMatrix
        # Dense matrix: direct row access
        for k in 1:length(I_coeff)
            V_port += Z[rwgID, k] * I_coeff[k]
        end
    else
        # For MLFMA or other operators, use matrix-vector product
        # Z_row = Z * e_rwgID (unit vector at port location)
        e_port = zeros(CT, length(I_coeff))
        e_port[rwgID] = one(CT)
        Z_row = similar(e_port)
        
        # Apply Z to unit vector to get the row
        # This works for any LinearMap-like operator
        mul!(Z_row, Z, e_port)
        
        V_port = dot(Z_row, I_coeff)
    end
    
    return V_port
end

"""
    extractPortVoltage(port::CurrentProbe{FT,IT}, Z, I_coeff::Vector{CT}) where {CT,FT,IT}

Extract voltage for a CurrentProbe (current source).

For current sources, we must compute the voltage from the impedance matrix:
V_port = Z[port.rwgID, :] ⋅ I_coeff
"""
function extractPortVoltage(
    port::CurrentProbe{FT,IT},
    Z,
    I_coeff::Vector{CT}
) where {CT<:Complex, FT<:Real, IT<:Integer}
    rwgID = port.rwgID
    
    if rwgID <= 0 || rwgID > length(I_coeff)
        error("Invalid port rwgID: $(rwgID)")
    end
    
    # Compute voltage from impedance matrix
    V_port = zero(CT)
    
    if Z isa AbstractMatrix
        for k in 1:length(I_coeff)
            V_port += Z[rwgID, k] * I_coeff[k]
        end
    else
        # MLFMA or operator
        e_port = zeros(CT, length(I_coeff))
        e_port[rwgID] = one(CT)
        Z_row = similar(e_port)
        mul!(Z_row, Z, e_port)
        V_port = dot(Z_row, I_coeff)
    end
    
    return V_port
end

"""
    extractPortVoltage(port::DeltaGapArrayPort{FT,IT,DT}, Z, I_coeff::Vector{CT}) where {CT,FT,IT,DT}

Extract voltage for a DeltaGapArrayPort (multi-edge port).

Computes weighted sum of voltages on all boundary edges:
V_total = Σ_edge weight_edge × V_edge
"""
function extractPortVoltage(
    port::DeltaGapArrayPort{FT,IT,DT},
    Z,
    I_coeff::Vector{CT}
) where {CT<:Complex, FT<:Real, IT<:Integer, DT}
    if !port.isBound
        error("Port must be bound to mesh")
    end
    
    V_total = zero(CT)
    
    for i in eachindex(port.rwgIDs)
        rwgID = port.rwgIDs[i]
        weight = port.edgeWeights[i]
        
        if rwgID <= 0 || rwgID > length(I_coeff)
            continue
        end
        
        # Compute voltage at this edge
        V_edge = zero(CT)
        if Z isa AbstractMatrix
            for k in 1:length(I_coeff)
                V_edge += Z[rwgID, k] * I_coeff[k]
            end
        else
            # MLFMA
            e_edge = zeros(CT, length(I_coeff))
            e_edge[rwgID] = one(CT)
            Z_row = similar(e_edge)
            mul!(Z_row, Z, e_edge)
            V_edge = dot(Z_row, I_coeff)
        end
        
        V_total += weight * V_edge
    end
    
    return V_total
end

"""
    extractPortCurrent(port::DeltaGapPort{FT,IT}, I_coeff::Vector{CT}) where {CT,FT,IT}

Extract port current from solution vector.

Returns I_coeff[port.rwgID] for single-edge ports.
"""
function extractPortCurrent(port::DeltaGapPort{FT,IT}, I_coeff::Vector{CT}) where {CT,FT,IT}
    rwgID = port.rwgID
    
    if rwgID <= 0 || rwgID > length(I_coeff)
        error("Invalid port rwgID: $(rwgID)")
    end
    
    return I_coeff[rwgID]
end

"""
    extractPortCurrent(port::CurrentProbe{FT,IT}, I_coeff::Vector{CT}) where {CT,FT,IT}

Extract port current from solution vector for CurrentProbe.

Returns I_coeff[port.rwgID] for single-edge ports.
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

Sums currents on all boundary edges with their weights and length factors:
I_total = Σ_edge weight_edge × I_edge × length_factor
"""
function extractPortCurrent(
    port::DeltaGapArrayPort{FT,IT,DT},
    I_coeff::Vector{CT}
) where {CT<:Complex, FT<:Real, IT<:Integer, DT}
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
        
        # Length factor: l/2 for full RWG, l for half RWG
        length_factor = port.triID_neg[i] > 0 ? 
                        port.edgeLengths[i] / 2 : 
                        port.edgeLengths[i]
        
        I_total += weight * I_coeff[rwgID] * length_factor
    end
    
    return I_total
end

"""
    measurePortImpedance(port_i::PortType, Z, I_j::Vector{CT}, port_j::PortType) where CT

Measure port impedance Z_ij = V_i / I_j.

Computes the impedance seen at port i when port j is excited.

# Arguments
- `port_i::PortType`: Measurement port (port i)
- `Z`: Impedance matrix or operator
- `I_j::Vector`: Solution current vector from exciting port j
- `port_j::PortType`: Excitation port (port j)

# Returns
- `Complex{FT}`: Port impedance Z_ij
"""
function measurePortImpedance(
    port_i::PortType,
    Z,
    I_j::Vector{CT},
    port_j::PortType
) where {CT<:Complex}
    error("measurePortImpedance not implemented for port types $(typeof(port_i)), $(typeof(port_j))")
end

"""
    measurePortImpedance(port_i::DeltaGapPort{FT,IT}, Z, I_j::Vector{CT}, port_j::PortType) where {CT,FT,IT}

Measure impedance for DeltaGapPort (single-edge voltage port).
"""
function measurePortImpedance(
    port_i::DeltaGapPort{FT,IT},
    Z,
    I_j::Vector{CT},
    port_j::PortType
) where {CT<:Complex, FT<:Real, IT<:Integer}
    V_i = extractPortVoltage(port_i, Z, I_j)
    I_i = extractPortCurrent(port_i, I_j)
    
    if abs(I_i) < eps(FT)
        return Complex{FT}(Inf, 0)
    end
    
    return V_i / I_i
end

"""
    measurePortImpedance(port_i::CurrentProbe{FT,IT}, Z, I_j::Vector{CT}, port_j::PortType) where {CT,FT,IT}

Measure impedance for CurrentProbe (single-edge current port).
"""
function measurePortImpedance(
    port_i::CurrentProbe{FT,IT},
    Z,
    I_j::Vector{CT},
    port_j::PortType
) where {CT<:Complex, FT<:Real, IT<:Integer}
    V_i = extractPortVoltage(port_i, Z, I_j)
    I_i = extractPortCurrent(port_i, I_j)
    
    if abs(I_i) < eps(FT)
        return Complex{FT}(Inf, 0)
    end
    
    return V_i / I_i
end

"""
    measurePortImpedance(port_i::DeltaGapArrayPort{FT,IT,DT}, Z, I_j::Vector{CT}, port_j::PortType) where {CT,FT,IT,DT}

Measure impedance for DeltaGapArrayPort (multi-edge port).
"""
function measurePortImpedance(
    port_i::DeltaGapArrayPort{FT,IT,DT},
    Z,
    I_j::Vector{CT},
    port_j::PortType
) where {CT<:Complex, FT<:Real, IT<:Integer, DT}
    V_i = extractPortVoltage(port_i, Z, I_j)
    I_i = extractPortCurrent(port_i, I_j)
    
    if abs(I_i) < eps(FT)
        return Complex{FT}(Inf, 0)
    end
    
    return V_i / I_i
end

"""
    measurePortImpedance(port_i::RectangularEdgePort{FT,IT,DT}, Z, I_j::Vector{CT}, port_j::PortType) where {CT,FT,IT,DT}

Measure impedance for RectangularEdgePort.

Delegates to DeltaGapArrayPort implementation via conversion.
"""
function measurePortImpedance(
    port_i::RectangularEdgePort{FT,IT,DT},
    Z,
    I_j::Vector{CT},
    port_j::PortType
) where {CT<:Complex, FT<:Real, IT<:Integer, DT}
    # Convert to DeltaGapArrayPort and measure
    gap_port = _to_delta_gap_array_port(port_i)
    return measurePortImpedance(gap_port, Z, I_j, port_j)
end

# =============================================================================
# Exports
# =============================================================================

export extractPortVoltage, extractPortCurrent, measurePortImpedance
