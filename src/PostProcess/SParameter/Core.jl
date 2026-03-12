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
    computeSParameters(
        ports::PortArray{FT,IT,PT},
        Z,
        nbf::Integer;
        solverT::Symbol = :direct,
        Z0::FT = FT(50.0),
        check_quality::Bool = true,
        kwargs...
    ) where {FT<:Real, IT<:Integer, PT<:PortType}

Compute S-parameter matrix for multi-port system.

This is the main entry point for S-parameter calculation. It works with both
dense impedance matrices (for small problems) and MLFMA operators (for large
problems), providing full flexibility in solver choice.

# Algorithm
1. For each port j (excitation port):
   a. Build excitation vector V_j with port j active
   b. Solve: Z × I_j = V_j (using specified solver)
   c. For each port i (measurement port):
      - Compute Z_port[i,j] = V_i / I_j
2. Convert Z_port to S-parameters: S = (Z_port - Z0⋅I) / (Z_port + Z0⋅I)

# Arguments
- `ports::PortArray`: Collection of ports
- `Z`: Impedance matrix or operator. Can be:
  - `Matrix{Complex{FT}}`: Dense matrix for direct solver
  - `MLFMAIterator`: MLFMA operator for iterative solver
  - Any type implementing the `solve()` interface
- `nbf::Integer`: Number of basis functions
- `solverT::Symbol`: Solver type (:direct, :gmres, :bicgstab, etc.)
- `Z0::Real`: Reference impedance (default: 50 Ω)
- `check_quality::Bool`: Enable quality checks on result
- `kwargs...`: Additional solver parameters (rtol, restart, Pl, Pr, etc.)

# Returns
- `SParameterResult`: Struct containing S-matrix, Z_port, and metadata

# Examples

## Direct solver (small problems < 10K unknowns)
```julia
Zmat = getImpedanceMatrix(geosInfo, nbf)
sparams = computeSParameters(ports, Zmat, nbf; solverT=:direct)
```

## GMRES with preconditioner (medium problems)
```julia
Zmat = getImpedanceMatrix(geosInfo, nbf)
sparams = computeSParameters(ports, Zmat, nbf; 
                            solverT=:gmres, 
                            rtol=1e-4, 
                            restart=50)
```

## MLFMA (large problems > 100K unknowns)
```julia
Zopt = MLFMAIterator(ZnearCSC, octree, geosInfo, bfsInfo)
Zprel = sparseApproximateInversePl(ZnearCSC, leafLevel)
sparams = computeSParameters(ports, Zopt, nbf;
                            solverT=:gmres,
                            Pl=Zprel,
                            rtol=1e-3)
```

## Accessing results
```julia
S = sparams.S[:,:,1]           # Full S-matrix at first frequency
S11 = sparams.S[1,1,1]         # S11 at first frequency
S21_dB = 20*log10(abs(S[2,1,1]))  # S21 in dB
```
"""
function computeSParameters(
    ports::PortArray{FT,IT,PT},
    Z,
    nbf::Integer;
    solverT::Symbol = :direct,
    Z0::FT = FT(50.0),
    check_quality::Bool = true,
    kwargs...
) where {FT<:Real, IT<:Integer, PT<:PortType}
    
    CT = Complex{FT}
    num_ports = ports.numPorts
    
    # Validate port count
    if num_ports == 0
        error("PortArray is empty")
    end
    
    # Initialize port impedance matrix
    Z_port = Matrix{CT}(undef, num_ports, num_ports)
    
    # Progress tracking
    pmeter = Progress(num_ports; desc="S-parameters: solving $num_ports ports")
    
    # For each port: excite → solve → measure all ports
    for j in 1:num_ports
        port_j = ports.ports[j]
        
        # Build excitation vector for port j
        V_j = buildExcitationVector(ports, j, nbf)
        
        # Solve: Z × I_j = V_j
        try
            I_j, converged = solve(Z, V_j; solverT=solverT, kwargs...)
            
            if !converged
                @warn "Solver did not converge for port $j excitation"
            end
            
            # Measure all ports with this excitation
            for i in 1:num_ports
                Z_port[i,j] = measurePortImpedance(
                    ports.ports[i], Z, I_j, ports.ports[j]
                )
            end
            
        catch e
            @warn "Failed to solve for port $j excitation: $e"
            for i in 1:num_ports
                Z_port[i,j] = CT(NaN, NaN)
            end
        end
        
        next!(pmeter)
    end
    
    # Convert Z_port to S-parameters
    S = z2s(Z_port, Z0)
    
    # Quality check
    if check_quality && !any(isnan, S)
        check_sparameter_quality(S)
    end
    
    # Build result
    return SParameterResult(
        S = reshape(S, size(S, 1), size(S, 2), 1),
        Z_port = reshape(Z_port, size(Z_port, 1), size(Z_port, 2), 1),
        frequencies = [Params.freq],
        Z0 = Z0,
        portIDs = [p.id for p in ports.ports]
    )
end

"""
    computeS11(port::PortType, Z, nbf::Integer; kwargs...)

Compute S11 for a single port.

Convenience function for single-port analysis.

# Arguments
- `port::PortType`: Single port
- `Z`: Impedance matrix or operator
- `nbf::Integer`: Number of basis functions
- `kwargs...`: Solver parameters (passed to computeSParameters)

# Returns
- `Complex{FT}`: S11 value

# Example
```julia
port = DeltaGapPort(...)
S11 = computeS11(port, Zmat, nbf; solverT=:direct)
S11_dB = 20*log10(abs(S11))
```
"""
function computeS11(
    port::DeltaGapPort{FT,IT},
    Z,
    nbf::Integer;
    Z0::FT = FT(50.0),
    kwargs...
) where {FT<:Real, IT<:Integer}
    ports = PortArray([port])
    sparams = computeSParameters(ports, Z, nbf; Z0=Z0, kwargs...)
    return sparams.S[1,1,1]
end

"""
    computeInputImpedance(port::PortType, Z, nbf::Integer; kwargs...)

Compute input impedance for a single port.

# Arguments
- `port::PortType`: Single port
- `Z`: Impedance matrix or operator
- `nbf::Integer`: Number of basis functions
- `kwargs...`: Solver parameters

# Returns
- `Complex{FT}`: Input impedance Zin

# Example
```julia
port = DeltaGapPort(...)
Zin = computeInputImpedance(port, Zmat, nbf)
```
"""
function computeInputImpedance(
    port::DeltaGapPort{FT,IT},
    Z,
    nbf::Integer;
    Z0::FT = FT(50.0),
    kwargs...
) where {FT<:Real, IT<:Integer}
    ports = PortArray([port])
    sparams = computeSParameters(ports, Z, nbf; Z0=Z0, kwargs...)
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
