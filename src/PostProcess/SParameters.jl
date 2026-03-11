# =============================================================================
# S-Parameter Calculation Module
# =============================================================================
#
# This module provides S-parameter extraction for MoM simulations.
#
# Key Formulas:
#   Input Impedance:     Zin = V_port / I_port
#   Reflection Coeff:    S11 = (Zin - Z0) / (Zin + Z0)
#   S-matrix:            S = (Z_port - Z0*I) * (Z_port + Z0*I)^(-1)
#
# =============================================================================

# -----------------------------------------------------------------------------
# DeltaGapPort (Voltage Source)
# -----------------------------------------------------------------------------

"""
    computeInputImpedance(port::DeltaGapPort, Z_matrix, V_excitation)

Compute input impedance for a DeltaGapPort.

Formula: Zin = V_port / I_port
"""
function computeInputImpedance(
    port::DeltaGapPort{FT, IT},
    Z_matrix::Matrix{Complex{FT}},
    V_excitation::Vector{Complex{FT}}
) where {FT<:Real, IT<:Integer}

    nbf = size(Z_matrix, 1)

    if port.rwgID <= 0 || port.rwgID > nbf
        error("Invalid port rwgID: $(port.rwgID), must be between 1 and $nbf")
    end

    # Solve: Z × I = V
    I = Z_matrix \ V_excitation

    # Get current at port location
    I_port = I[port.rwgID]

    if abs(I_port) < eps(FT)
        error("Zero current at port $(port.id). Cannot compute input impedance.")
    end

    # Zin = V_port / I_port
    Zin = port.V / I_port

    return Zin
end


"""
    computeS11(port::DeltaGapPort, Z_matrix, V_excitation; Z0=50.0)

Compute S11 for a DeltaGapPort.

Formula: S11 = (Zin - Z0) / (Zin + Z0)
"""
function computeS11(
    port::DeltaGapPort{FT, IT},
    Z_matrix::Matrix{Complex{FT}},
    V_excitation::Vector{Complex{FT}};
    Z0::FT = FT(50.0)
) where {FT<:Real, IT<:Integer}

    Zin = computeInputImpedance(port, Z_matrix, V_excitation)
    S11 = (Zin - Z0) / (Zin + Z0)
    return S11
end


# -----------------------------------------------------------------------------
# CurrentProbe (Current Source)
# -----------------------------------------------------------------------------

"""
    computeInputImpedance(port::CurrentProbe, Z_matrix, V_excitation)

Compute input impedance for a CurrentProbe.

Formula: Zin = V_port / I_probe
where V_port is computed from Z × I.
"""
function computeInputImpedance(
    port::CurrentProbe{FT, IT},
    Z_matrix::Matrix{Complex{FT}},
    V_excitation::Vector{Complex{FT}}
) where {FT<:Real, IT<:Integer}

    nbf = size(Z_matrix, 1)

    if port.rwgID <= 0 || port.rwgID > nbf
        error("Invalid port rwgID: $(port.rwgID)")
    end

    # Solve: Z × I_coeff = V_exc
    I_coeff = Z_matrix \ V_excitation

    # Compute resulting voltage at port
    V_port = zero(Complex{FT})
    for k in 1:nbf
        V_port += Z_matrix[port.rwgID, k] * I_coeff[k]
    end

    # Zin = V_port / I_probe
    Zin = V_port / port.I

    return Zin
end


"""
    computeS11(port::CurrentProbe, Z_matrix, V_excitation; Z0=50.0)

Compute S11 for a CurrentProbe.
"""
function computeS11(
    port::CurrentProbe{FT, IT},
    Z_matrix::Matrix{Complex{FT}},
    V_excitation::Vector{Complex{FT}};
    Z0::FT = FT(50.0)
) where {FT<:Real, IT<:Integer}

    Zin = computeInputImpedance(port, Z_matrix, V_excitation)
    S11 = (Zin - Z0) / (Zin + Z0)
    return S11
end



# -----------------------------------------------------------------------------
# Multi-Port S-Parameter Matrix
# -----------------------------------------------------------------------------

"""
    computeSParameters(ports::PortArray, Z_MoM, V_excitation; Z0=50.0)

Compute full S-parameter matrix for a multi-port system.

Algorithm:
1. Build Z_port[i,j] = V_i / I_j for each port combination
2. S = (Z_port - Z0*I) / (Z_port + Z0*I)
"""
function computeSParameters(
    ports::PortArray{FT, IT, PT},
    Z_MoM::Matrix{Complex{FT}},
    V_excitation::Vector{Complex{FT}};
    Z0::FT = FT(50.0),
    check_quality::Bool = true
) where {FT<:Real, IT<:Integer, PT<:PortType}

    nbf = size(Z_MoM, 1)
    num_ports = ports.numPorts

    # Single port case
    if num_ports == 1
        port = ports.ports[1]
        S11 = computeS11(port, Z_MoM, V_excitation; Z0 = Z0)
        return S11
    end

    # Multi-port: Build port impedance matrix
    Z_port = Matrix{Complex{FT}}(undef, num_ports, num_ports)

    for j in 1:num_ports
        port_j = ports.ports[j]
        pid_j = port_j.rwgID

        if pid_j <= 0 || pid_j > nbf
            error("Invalid port rwgID: $pid_j for port $(port_j.id)")
        end

        # Scaling factor for port j
        scale_j = if isa(port_j, DeltaGapPort)
            (port_j.triID_neg > 0) ? (port_j.edgel / 2) : port_j.edgel
        else
            FT(1.0)
        end

        # Create excitation for port j
        V_j = zeros(Complex{FT}, nbf)
        V_j[pid_j] = Complex{FT}(1.0) * scale_j

        # Solve: Z × I_j = V_j
        try
            I_j = Z_MoM \ V_j

            if abs(I_j[pid_j]) < eps(FT)
                @warn "Zero current at excited port $j"
                for i in 1:num_ports
                    Z_port[i, j] = Complex{FT}(Inf, 0)
                end
                continue
            end

            # Compute impedance for all ports
            for i in 1:num_ports
                port_i = ports.ports[i]
                pid_i = port_i.rwgID

                if pid_i <= 0 || pid_i > nbf
                    error("Invalid port rwgID: $pid_i for port $(port_i.id)")
                end

                scale_i = if isa(port_i, DeltaGapPort)
                    (port_i.triID_neg > 0) ? (port_i.edgel / 2) : port_i.edgel
                else
                    FT(1.0)
                end

                # Compute voltage response at port i
                V_exc_response = zero(Complex{FT})
                for k in 1:nbf
                    V_exc_response += Z_MoM[pid_i, k] * I_j[k]
                end
                V_port_i = V_exc_response / scale_i

                # Z_port[i,j] = V_i / I_j
                Z_port[i, j] = V_port_i / I_j[pid_j]
            end
        catch e
            @warn "Failed to solve for port $j excitation: $e"
            for i in 1:num_ports
                Z_port[i, j] = Complex{FT}(NaN, NaN)
            end
        end
    end

    # Convert to S-parameters: S = (Z - Z0*I) / (Z + Z0*I)
    Z0_matrix = Z0 * Matrix{Complex{FT}}(I, num_ports, num_ports)

    try
        S_matrix = (Z_port - Z0_matrix) / (Z_port + Z0_matrix)
        
        if check_quality && !any(isnan, S_matrix)
            check_sparameter_quality(S_matrix)
        end
        
        return S_matrix
    catch e
        @warn "Failed to compute S-parameter matrix: $e"
        return fill(Complex{FT}(NaN, NaN), num_ports, num_ports)
    end
end


"""
    computeSParameters(ports::PortArray, Z_MoM; Z0=50.0)

Compute S-parameters using only the impedance matrix.
"""
function computeSParameters(
    ports::PortArray{FT, IT, PT},
    Z_MoM::Matrix{Complex{FT}};
    Z0::FT = FT(50.0),
    check_quality::Bool = true
) where {FT<:Real, IT<:Integer, PT<:PortType}

    nbf = size(Z_MoM, 1)
    V_dummy = zeros(Complex{FT}, nbf)
    
    return computeSParameters(ports, Z_MoM, V_dummy; Z0 = Z0, check_quality = check_quality)
end


# -----------------------------------------------------------------------------
# Helper Functions
# -----------------------------------------------------------------------------

is_voltage_port(port::DeltaGapPort) = true
is_voltage_port(port::CurrentProbe) = false
is_voltage_port(port::T) where {T} = port isa DeltaGapPort


# -----------------------------------------------------------------------------
# Exports
# -----------------------------------------------------------------------------

export  computeInputImpedance,
        computeS11,
        computeSParameters,
        is_voltage_port
