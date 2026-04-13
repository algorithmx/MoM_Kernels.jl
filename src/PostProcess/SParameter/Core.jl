# =============================================================================
# Core S-Parameter Computation
# =============================================================================
#
# Main computeSParameters function that orchestrates N port excitations,
# solves, and measurements to build the S-parameter matrix.
#
# Designed to work with both dense matrices (direct solver) and 
# MLFMA operators (iterative solver).
#
# =============================================================================

using ProgressMeter

"""
    _port_numeric_types(port::PortType) -> (Type, Type)

Return the float precision type `FT` and integer type `IT` of a port by
inspecting its concrete type parameters.

All port types follow the convention that the first type parameter is
`FT<:Real` and the second is `IT<:Integer`. The function falls back to
`(Float64, Int)` when the parameters cannot be determined.
"""
function _port_numeric_types(port::PortType)
    T      = typeof(port)
    params = T.parameters
    FT = (length(params) >= 1 && params[1] isa Type && params[1] <: Real)    ? params[1] : Float64
    IT = (length(params) >= 2 && params[2] isa Type && params[2] <: Integer) ? params[2] : Int
    return FT, IT
end

"""
        ports::PortArray{FT,IT,PT},
        Z,
        nbf::Integer;
        solverT      = :direct,
        Z0           = 50.0,
        ieType       = :efie,
        trianglesInfo = [],
        cfie_alpha   = 0.5,
        check_quality = true,
        kwargs...
    )

Compute S-parameter matrix for multi-port system.

This is the main entry point for S-parameter calculation. It works with both
dense impedance matrices (for small problems) and MLFMA operators (for large
problems), providing full flexibility in solver choice.

# Algorithm

## Y-Matrix Approach
This function uses the short-circuit admittance (Y) matrix method for robust
multi-port S-parameter calculation:

1. For each port j (excitation port):
   - Build excitation vector V_j with port j active (using `ieType`)
   - Solve: Z × I_j = V_j (using specified solver)
   - Compute admittance column: Y[i,j] = I_i^port / V_j^src
     where I_i^port = extractPortCurrent(port_i, I_j)

2. Invert Y-matrix: Z_port = inv(Y_port)

3. Convert to S-parameters: S = (Z_port - Z0⋅I) / (Z_port + Z0⋅I)

## Why Y-Matrix?
The naive approach of computing V_i = (Z·I_j)[rwgID_i] fails for multi-port
analysis because Z·I_j = V_exc_j exactly (MoM identity), and V_exc_j[rwgID_i] = 0
for any passive port i≠j — so that approach always yields zero off-diagonal.

The Y-matrix method correctly captures coupling by measuring the actual
current at each port when another port is excited.

# Arguments
- `ports::PortArray`   : collection of ports
- `Z`                  : impedance matrix or operator
- `nbf::Integer`       : number of basis functions
- `solverT::Symbol`    : solver type (`:direct`, `:gmres`, `:bicgstab`, …)
- `Z0::Real`           : reference impedance (default: 50 Ω)
- `ieType::Symbol`     : integral equation type — `:efie` (default), `:mfie`, `:cfie`.
  Selects the excitation vector formula. Voltage ports (DeltaGap*) always
  use EFIE excitation regardless of this setting; `CurrentProbe` in CFIE mode
  applies the α-weighted combination.
- `trianglesInfo`      : mesh triangle data forwarded to excitation functions.
  May be omitted (defaults to `[]`) for port types that store all geometry on
  the struct itself.
- `cfie_alpha::Real`   : CFIE combination weight α (default: 0.5).
  Only effective when `ieType == :cfie` with `CurrentProbe` ports.
- `check_quality::Bool`: enable passivity / reciprocity check on result
- `kwargs...`          : additional solver parameters (rtol, restart, Pl, Pr, …)

# Returns
- `SParameterResult`: struct containing S-matrix, Z_port, and metadata

# Heterogeneous Port Arrays
The function supports heterogeneous `PortArray`s containing different port types
(e.g., mixing `DeltaGapPort`, `CurrentProbe`, `DeltaGapArrayPort`). All ports
must share the same float (`FT`) and integer (`IT`) type parameters.

# Examples

## EFIE, direct solver (small problems < 10K unknowns)
```julia
sparams = computeSParameters(ports, Zmat, nbf; solverT=:direct)
```

## CFIE, GMRES (resonant / closed cavity problems)
```julia
sparams = computeSParameters(ports, Zmat, nbf;
    ieType = :cfie, cfie_alpha = 0.5,
    trianglesInfo = trianglesInfo,
    solverT = :gmres, rtol = 1e-4)
```

## MLFMA (large problems > 100K unknowns)
```julia
sparams = computeSParameters(ports, Zopt, nbf;
    solverT = :gmres, Pl = Zprel, rtol = 1e-3)
```
"""
function computeSParameters(
    ports::PortArray{FT,IT,PT},
    Z,
    nbf::Integer;
    solverT::Symbol   = :direct,
    Z0::FT            = FT(50.0),
    ieType::Symbol    = :efie,
    trianglesInfo     = TriangleInfo{IT,FT}[],
    cfie_alpha::Real  = FT(0.5),
    check_quality::Bool = true,
    kwargs...
) where {FT<:Real, IT<:Integer, PT<:PortType}
    
    CT = Complex{FT}
    num_ports = ports.numPorts
    
    # Validate port count
    if num_ports == 0
        error("PortArray is empty")
    end
    
    # Validate all ports have consistent FT and IT types (for heterogeneous arrays)
    if PT === PortType  # Heterogeneous array
        for (i, port) in enumerate(ports.ports)
            pFT, pIT = _port_numeric_types(port)
            if pFT !== FT
                error("Port $i has different float type ($pFT) than PortArray ($FT). " *
                      "All ports must share the same FT and IT type parameters.")
            end
            if pIT !== IT
                error("Port $i has different integer type ($pIT) than PortArray ($IT). " *
                      "All ports must share the same FT and IT type parameters.")
            end
        end
    end
    
    # Short-circuit admittance matrix Y_port[i,j] = I_i^port / V_j^src.
    #
    # Each solve (port j excited, all others short-circuited V=0) directly
    # gives column j of Y via Y[i,j] = extractPortCurrent(port_i, I_j) / V_j^src.
    # The port impedance matrix is then Z_port = inv(Y_port), and S = z2s(Z_port).
    #
    # Why not read (Z·I_j)[rwgID_i] for the voltage at port i?
    # Because Z·I_j = V_exc_j exactly (MoM identity), and V_exc_j[rwgID_i] = 0
    # for any passive port i≠j — so that approach always yields zero off-diagonal.
    Y_port = Matrix{CT}(undef, num_ports, num_ports)
    
    # Progress tracking
    pmeter = Progress(num_ports; desc="S-parameters: solving $num_ports ports")
    
    # For each port: excite → solve → collect admittance column
    for j in 1:num_ports
        port_j = ports.ports[j]
        
        # Build excitation vector for port j
        V_j = buildExcitationVector(ports, j, nbf;
                  ieType=ieType, trianglesInfo=trianglesInfo,
                  cfie_alpha=cfie_alpha)
        
        # Solve: Z × I_j = V_j
        try
            I_j, converged = solve(Z, V_j; solverT=solverT, kwargs...)
            
            # converged is nothing for direct solver (always converged)
            # converged is a ConvergenceHistory for iterative solvers
            if converged !== nothing && !converged.isconverged
                @warn "Solver did not converge for port $j excitation"
            end
            
            # Admittance column j: Y[i,j] = extractPortCurrent(port_i, I_j) / V_j^src
            V_src_j = port_j.V
            for i in 1:num_ports
                Y_port[i,j] = extractPortCurrent(ports.ports[i], I_j) / V_src_j
            end
            
        catch e
            @warn "Failed to solve for port $j excitation: $e"
            for i in 1:num_ports
                Y_port[i,j] = CT(NaN, NaN)
            end
        end
        
        next!(pmeter)
    end
    
    # Invert admittance matrix → port impedance matrix Z_port = Y_port⁻¹
    Z_port = any(isnan, Y_port) ? fill(CT(NaN,NaN), num_ports, num_ports) : inv(Y_port)
    
    # Convert Z_port to S-parameters
    S = z2s(Z_port, Z0)
    
    # Quality check
    if check_quality && !any(isnan, S)
        check_sparameter_quality(S)
    end
    
    # Build result
    S_3d = reshape(S, size(S, 1), size(S, 2), 1)
    Z_port_3d = reshape(Z_port, size(Z_port, 1), size(Z_port, 2), 1)
    # Use FT to ensure type consistency (Params.frequency may be different precision)
    freqs = FT[Params.frequency]
    return SParameterResult(S_3d, Z_port_3d, freqs, Z0, [p.id for p in ports.ports])
end

"""
    computeS11(port::PortType, Z, nbf::Integer; Z0=50.0, kwargs...)

Compute S11 for a single port.

Convenience function for single-port analysis.  Accepts any concrete subtype
of `PortType`, including `DeltaGapPort`, `DeltaGapArrayPort` (with any
excitation distribution), `RectangularWaveguidePort`, and `CurrentProbe`.

# Arguments
- `port::PortType`: Single port
- `Z`: Impedance matrix or operator
- `nbf::Integer`: Number of basis functions
- `Z0::Real`: Reference impedance (default: 50 Ω).  Converted to the port's
  float precision automatically.
- `kwargs...`: Solver parameters (passed to `computeSParameters`)

# Returns
- `Complex{FT}`: S11 value

# Example
```julia
port = DeltaGapArrayPort(; ...)
bind_to_mesh!(port, pred, trianglesInfo, rwgsInfo)
S11 = computeS11(port, Zmat, nbf; solverT=:direct)
S11_dB = 20*log10(abs(S11))
```
"""
function computeS11(
    port::PortType,
    Z,
    nbf::Integer;
    Z0::Real          = 50.0,
    ieType::Symbol    = :efie,
    trianglesInfo     = [],
    cfie_alpha::Real  = 0.5,
    kwargs...
)
    T  = typeof(port)
    FT = (length(T.parameters) >= 1 && T.parameters[1] isa Type &&
          T.parameters[1] <: Real) ? T.parameters[1] : Float64
    ports = PortArray([port])
    sparams = computeSParameters(ports, Z, nbf;
                  Z0=FT(Z0), ieType=ieType,
                  trianglesInfo=trianglesInfo,
                  cfie_alpha=FT(cfie_alpha), kwargs...)
    return sparams.S[1,1,1]
end

"""
    computeInputImpedance(port::PortType, Z, nbf::Integer; Z0=50.0, kwargs...)

Compute input impedance for a single port.

Accepts any concrete subtype of `PortType`; `FT` is inferred from the port.

# Arguments
- `port::PortType`: Single port
- `Z`: Impedance matrix or operator
- `nbf::Integer`: Number of basis functions
- `Z0::Real`: Reference impedance (default: 50 Ω)
- `kwargs...`: Solver parameters

# Returns
- `Complex{FT}`: Input impedance Zin

# Example
```julia
port = DeltaGapArrayPort(; ...)
bind_to_mesh!(port, pred, trianglesInfo, rwgsInfo)
Zin = computeInputImpedance(port, Zmat, nbf)
```
"""
function computeInputImpedance(
    port::PortType,
    Z,
    nbf::Integer;
    Z0::Real          = 50.0,
    ieType::Symbol    = :efie,
    trianglesInfo     = [],
    cfie_alpha::Real  = 0.5,
    kwargs...
)
    T  = typeof(port)
    FT = (length(T.parameters) >= 1 && T.parameters[1] isa Type &&
          T.parameters[1] <: Real) ? T.parameters[1] : Float64
    ports = PortArray([port])
    sparams = computeSParameters(ports, Z, nbf;
                  Z0=FT(Z0), ieType=ieType,
                  trianglesInfo=trianglesInfo,
                  cfie_alpha=FT(cfie_alpha), kwargs...)
    return sparams.Z_port[1,1,1]
end

"""
    check_sparameter_quality(S::Matrix{CT}) where CT

Perform quality checks on S-parameter matrix.

Checks:
- Passivity: I - S'⋅S should be positive semi-definite
- Reciprocity (for reciprocal media): S should be symmetric

# Arguments
- `S::Matrix`: S-parameter matrix

# Returns
- `Bool`: true if quality checks pass

# Warnings
Issues warnings if passivity or reciprocity violations are detected.
"""
function check_sparameter_quality(S::Matrix{CT}) where {CT<:Complex}
    FT = real(CT)
    n = size(S, 1)
    
    # Check passivity: I - S'⋅S ≥ 0
    I_mat = Matrix{CT}(I, n, n)
    G = I_mat - S' * S  # Power gain matrix
    
    # Check eigenvalues
    eigs = eigvals(Hermitian(G))
    min_eig = minimum(real.(eigs))
    
    if min_eig < -sqrt(eps(FT))
        @warn "Passivity violation detected: minimum eigenvalue = $min_eig"
        return false
    end
    
    # Check reciprocity for 2-port (S12 should equal S21)
    if n == 2
        asymmetry = abs(S[1,2] - S[2,1])
        if asymmetry > sqrt(eps(FT))
            @warn "Reciprocity violation: |S12 - S21| = $asymmetry"
        end
    end
    
    return true
end

# =============================================================================
# Backward Compatibility Functions
# =============================================================================

"""
    computeSParameters(ports::PortArray, Z_MoM::Matrix{CT}, V_excitation::Vector{CT}; kwargs...) where CT

Backward-compatible version that accepts V_excitation (ignored).

This signature matches the old API for compatibility with existing code.
"""
function computeSParameters(
    ports::PortArray{FT,IT,PT},
    Z_MoM::Matrix{CT},
    V_excitation::Vector{CT};
    kwargs...
) where {CT<:Complex, FT<:Real, IT<:Integer, PT<:PortType}
    nbf = size(Z_MoM, 1)
    return computeSParameters(ports, Z_MoM, nbf; kwargs...)
end

"""
    computeSParameters(ports::PortArray, Z_MoM::Matrix{CT}; kwargs...) where CT

Backward-compatible version with only impedance matrix.
"""
function computeSParameters(
    ports::PortArray{FT,IT,PT},
    Z_MoM::Matrix{CT};
    kwargs...
) where {CT<:Complex, FT<:Real, IT<:Integer, PT<:PortType}
    nbf = size(Z_MoM, 1)
    return computeSParameters(ports, Z_MoM, nbf; kwargs...)
end

# =============================================================================
# Exports
# =============================================================================

export computeSParameters, computeS11, computeInputImpedance, check_sparameter_quality
