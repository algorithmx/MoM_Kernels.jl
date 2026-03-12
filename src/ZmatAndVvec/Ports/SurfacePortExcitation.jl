# =============================================================================
# SurfacePortExcitation.jl - Port excitation vector implementations for surface cells
# =============================================================================
#
# This file implements excitation vector calculations for surface ports (RWG basis).
# Julia's multiple dispatch automatically selects these methods over the general
# ExcitingSources methods when the source is a PortType.
#
# Port Types Covered:
#   - DeltaGapPort: Single-edge voltage source (delta-gap model)
#   - CurrentProbe: Direct current injection
#   - DeltaGapArrayPort: Multi-edge array port with weighted distribution
#   - RectangularWaveguidePort: Rectangular waveguide port (delegates to DeltaGapArrayPort)
#
# IE Compatibility:
#   | Port Type             | EFIE | MFIE | CFIE |
#   |----------------------|------|------|------|
#   | DeltaGapPort         | ✓    | EFIE | ✓    |
#   | CurrentProbe         | ✓    | ✓    | ✓    |
#   | DeltaGapArrayPort    | ✓    | ERROR| ✓    |
#   | RectangularWaveguidePort | ✓ | ERROR| ✓    |
#
# Note: For MFIE with voltage ports (DeltaGapPort, DeltaGapArrayPort), we use EFIE
# excitation as industry standard, since voltage sources are naturally expressed
# in the EFIE formulation.
# =============================================================================

# =============================================================================
# EFIE - Surface Ports
# =============================================================================

# -----------------------------------------------------------------------------
# DeltaGapPort - Single edge voltage source
# -----------------------------------------------------------------------------

"""
    excitationVectorEFIE(port::DeltaGapPort, trianglesInfo, nbf)

Compute EFIE excitation vector for Delta-Gap port.

# Physical Formula
For Delta-Gap excitation, the excitation vector is:
```
V_m = V₀ × l_m / 2  (full RWG)
V_m = V₀ × l_m      (half RWG, boundary)
```

where V₀ is the port voltage and l_m is the edge length.

# Parameters
- `port::DeltaGapPort`: Delta-Gap port
- `trianglesInfo`: Triangle information array
- `nbf`: Total number of basis functions

# Returns
- Excitation vector (complex array of length nbf)
"""
function excitationVectorEFIE(
    port::DeltaGapPort{FT, IT},
    trianglesInfo::Vector{TriangleInfo{IT, FT}},
    nbf::Integer
) where {FT<:Real, IT<:Integer}

    V = zeros(Complex{FT}, nbf)
    excitationVectorEFIE!(V, port, trianglesInfo)
    return V
end


"""
    excitationVectorEFIE!(V, port::DeltaGapPort, trianglesInfo)

Add Delta-Gap port excitation to existing vector V.
"""
function excitationVectorEFIE!(
    V::Vector{Complex{FT}},
    port::DeltaGapPort{FT, IT},
    trianglesInfo::Vector{TriangleInfo{IT, FT}}
) where {FT<:Real, IT<:Integer}

    rwgID = port.rwgID
    (rwgID <= 0 || rwgID > length(V)) && return V

    # Full RWG: V × l/2, Half RWG: V × l
    factor = port.triID_neg > 0 ? (port.edgel / 2) : port.edgel
    V[rwgID] += port.V * factor

    return V
end


# -----------------------------------------------------------------------------
# CurrentProbe - Direct current injection
# -----------------------------------------------------------------------------

"""
    excitationVectorEFIE(probe::CurrentProbe, trianglesInfo, nbf)

Compute EFIE excitation vector for current probe.

For current probe, the excitation vector directly injects current:
```
V_m = I₀  (when m is the probe basis function)
V_m = 0   (otherwise)
```
"""
function excitationVectorEFIE(
    probe::CurrentProbe{FT, IT},
    trianglesInfo::Vector{TriangleInfo{IT, FT}},
    nbf::Integer
) where {FT<:Real, IT<:Integer}

    V = zeros(Complex{FT}, nbf)

    if probe.rwgID > 0 && probe.rwgID <= nbf
        V[probe.rwgID] = probe.I
    end

    return V
end


"""
    excitationVectorEFIE!(V, probe::CurrentProbe, trianglesInfo)

Add current probe excitation to existing vector V.
"""
function excitationVectorEFIE!(
    V::Vector{Complex{FT}},
    probe::CurrentProbe{FT, IT},
    trianglesInfo::Vector{TriangleInfo{IT, FT}}
) where {FT<:Real, IT<:Integer}

    if probe.rwgID > 0 && probe.rwgID <= length(V)
        V[probe.rwgID] += probe.I
    end

    return V
end


# -----------------------------------------------------------------------------
# DeltaGapArrayPort - Multi-edge with weights
# -----------------------------------------------------------------------------

"""
    excitationVectorEFIE(port::DeltaGapArrayPort, trianglesInfo, nbf)

Compute EFIE excitation vector for Delta-Gap array port.

Applies excitation to all RWG basis functions on the port perimeter identified
during `bind_to_mesh!`. Each boundary edge receives:

```
V[rwgID] += V_port × weight(edge) × edge_length / 2   (full RWG)
V[rwgID] += V_port × weight(edge) × edge_length       (half RWG, boundary)
```

Requires the port to be bound to mesh via `bind_to_mesh!` first.
"""
function excitationVectorEFIE(
    port::DeltaGapArrayPort{FT, IT, DT},
    trianglesInfo::Vector{TriangleInfo{IT, FT}},
    nbf::Integer
) where {FT<:Real, IT<:Integer, DT<:AbstractExcitationDistribution{FT}}

    port.isBound || error("Port must be bound to mesh before computing excitation")
    V = zeros(Complex{FT}, nbf)
    excitationVectorEFIE!(V, port, trianglesInfo)
    return V
end


"""
    excitationVectorEFIE!(V, port::DeltaGapArrayPort, trianglesInfo)

Add Delta-Gap array port excitation to existing vector V.
"""
function excitationVectorEFIE!(
    V::Vector{Complex{FT}},
    port::DeltaGapArrayPort{FT, IT, DT},
    trianglesInfo::Vector{TriangleInfo{IT, FT}}
) where {FT<:Real, IT<:Integer, DT<:AbstractExcitationDistribution{FT}}

    port.isBound || error("Port must be bound to mesh")

    if port.singleEdgeMode
        _excitation_single_edge!(V, port)
    else
        _excitation_array!(V, port)
    end

    return V
end


"""
    _excitation_array!(V, port)

Apply Delta-Gap excitation to all boundary edges (array excitation mode).

For each boundary edge identified during `bind_to_mesh!`, computes the excitation
contribution using the Delta-Gap formula: V × weight × length_factor.
"""
function _excitation_array!(
    V::Vector{Complex{FT}},
    port::DeltaGapArrayPort{FT, IT, DT}
) where {FT<:Real, IT<:Integer, DT<:AbstractExcitationDistribution{FT}}

    for i in eachindex(port.rwgIDs)
        rwgID = port.rwgIDs[i]
        (rwgID <= 0 || rwgID > length(V)) && continue

        factor = port.triID_neg[i] > 0 ? (port.edgeLengths[i] / 2) : port.edgeLengths[i]
        V[rwgID] += port.V * port.edgeWeights[i] * factor
    end
end


"""
    _excitation_single_edge!(V, port)

Apply Delta-Gap excitation to a single edge (like DeltaGapPort).

Uses the `primaryRwgID` field to select which edge to excite, falling back
to the first available edge if not found.
"""
function _excitation_single_edge!(
    V::Vector{Complex{FT}},
    port::DeltaGapArrayPort{FT, IT, DT}
) where {FT<:Real, IT<:Integer, DT<:AbstractExcitationDistribution{FT}}

    idx = findfirst(==(port.primaryRwgID), port.rwgIDs)
    idx === nothing && !isempty(port.rwgIDs) && (idx = 1)
    idx === nothing && return

    rwgID = port.rwgIDs[idx]
    (rwgID <= 0 || rwgID > length(V)) && return

    factor = port.triID_neg[idx] > 0 ? (port.edgeLengths[idx] / 2) : port.edgeLengths[idx]
    V[rwgID] += port.V * port.edgeWeights[idx] * factor
end


# -----------------------------------------------------------------------------
# RectangularWaveguidePort - Delegates to DeltaGapArrayPort
# -----------------------------------------------------------------------------

"""
    excitationVectorEFIE(port::RectangularWaveguidePort, trianglesInfo, nbf)

Compute EFIE excitation vector for rectangular waveguide port.

Delegates to DeltaGapArrayPort implementation.
"""
function excitationVectorEFIE(
    port::RectangularWaveguidePort{FT, IT, DT},
    trianglesInfo::Vector{TriangleInfo{IT, FT}},
    nbf::Integer
) where {FT<:Real, IT<:Integer, DT<:AbstractExcitationDistribution{FT}}

    V = zeros(Complex{FT}, nbf)
    excitationVectorEFIE!(V, port, trianglesInfo)
    return V
end


"""
    excitationVectorEFIE!(V, port::RectangularWaveguidePort, trianglesInfo)

Add rectangular waveguide port excitation to existing vector V.

Delegates to DeltaGapArrayPort implementation.
"""
function excitationVectorEFIE!(
    V::Vector{Complex{FT}},
    port::RectangularWaveguidePort{FT, IT, DT},
    trianglesInfo::Vector{TriangleInfo{IT, FT}}
) where {FT<:Real, IT<:Integer, DT<:AbstractExcitationDistribution{FT}}

    # Delegate to DeltaGapArrayPort implementation
    gap_port = _to_delta_gap_array_port(port)
    excitationVectorEFIE!(V, gap_port, trianglesInfo)
    return V
end


# =============================================================================
# MFIE - Surface Ports
# =============================================================================

# -----------------------------------------------------------------------------
# DeltaGapPort - Use EFIE (industry standard for voltage sources)
# -----------------------------------------------------------------------------

"""
    excitationVectorMFIE(port::DeltaGapPort, trianglesInfo, nbf; strategy)

Compute MFIE excitation vector for Delta-Gap port.

# Strategy Options
- `:efie_fallback` (default): Use EFIE excitation (industry standard)
- `:convert`: Convert voltage source to equivalent current source
- `:hybrid`: Use hybrid method with small EFIE contribution

# Note
DEPRECATED: Industry standard is to use EFIE excitation for voltage sources
in MFIE formulation. The `:efie_fallback` strategy is recommended.
"""
function excitationVectorMFIE(
    port::DeltaGapPort{FT, IT},
    trianglesInfo::Vector{TriangleInfo{IT, FT}},
    nbf::Integer;
    strategy::Symbol = :efie_fallback
) where {FT<:Real, IT<:Integer}

    if strategy == :efie_fallback
        # Industry standard: Use EFIE excitation for voltage sources
        @warn "MFIE with DeltaGapPort: Using EFIE excitation (industry standard)" maxlog=1
        return excitationVectorEFIE(port, trianglesInfo, nbf)

    elseif strategy == :convert
        # Strategy: Convert to equivalent current source
        # I_eq = V / Z_port, assuming 50 ohm port impedance
        Z_port = FT(50.0)
        I_eq = port.V / Z_port

        # Create equivalent current probe
        probe = CurrentProbe{FT, IT}(;
            id = port.id,
            I = I_eq,
            freq = port.freq,
            rwgID = port.rwgID,
            triID = port.triID_pos,
            edgel = port.edgel,
            center = port.center,
            isActive = port.isActive
        )

        return excitationVectorMFIE(probe, trianglesInfo, nbf)

    elseif strategy == :hybrid
        # Hybrid method: Add small EFIE contribution at port location
        V = zeros(Complex{FT}, nbf)

        if port.rwgID > 0 && port.rwgID <= nbf
            # Use small coefficient to avoid breaking MFIE numerical properties
            epsilon = FT(1e-3)
            V[port.rwgID] = port.V * port.edgel * epsilon
        end

        return V
    else
        error("Unknown MFIE excitation strategy: $strategy")
    end
end


# -----------------------------------------------------------------------------
# CurrentProbe - Natural for MFIE
# -----------------------------------------------------------------------------

"""
    excitationVectorMFIE(probe::CurrentProbe, trianglesInfo, nbf)

Compute MFIE excitation vector for current probe.

Current probe is natural for MFIE as it directly involves current.
```
V_m = I₀ (when m is the probe basis function)
V_m = 0  (otherwise)
```
"""
function excitationVectorMFIE(
    probe::CurrentProbe{FT, IT},
    trianglesInfo::Vector{TriangleInfo{IT, FT}},
    nbf::Integer
) where {FT<:Real, IT<:Integer}

    V = zeros(Complex{FT}, nbf)

    if probe.rwgID > 0 && probe.rwgID <= nbf
        V[probe.rwgID] = probe.I
    end

    return V
end


# -----------------------------------------------------------------------------
# DeltaGapArrayPort - Not supported
# -----------------------------------------------------------------------------

"""
    excitationVectorMFIE(port::DeltaGapArrayPort, trianglesInfo, nbf; strategy)

MFIE excitation for DeltaGapArrayPort is not implemented.

# Physical Consideration
MFIE formulation with multi-edge voltage ports requires conversion from
voltage excitation to equivalent magnetic current distribution.

Use EFIE formulation or implement a conversion strategy if needed.
"""
function excitationVectorMFIE(
    port::DeltaGapArrayPort{FT, IT, DT},
    trianglesInfo::Vector{TriangleInfo{IT, FT}},
    nbf::Integer;
    strategy::Symbol = :efie_fallback
) where {FT<:Real, IT<:Integer, DT<:AbstractExcitationDistribution{FT}}

    if strategy == :efie_fallback
        @warn "MFIE with DeltaGapArrayPort: Using EFIE excitation (industry standard)" maxlog=1
        return excitationVectorEFIE(port, trianglesInfo, nbf)
    else
        error("excitationVectorMFIE is not implemented for DeltaGapArrayPort with strategy=$strategy. " *
              "MFIE formulation with multi-edge voltage ports requires physical consideration. " *
              "Use EFIE formulation or strategy=:efie_fallback.")
    end
end


# -----------------------------------------------------------------------------
# RectangularWaveguidePort - Delegates to DeltaGapArrayPort
# -----------------------------------------------------------------------------

"""
    excitationVectorMFIE(port::RectangularWaveguidePort, trianglesInfo, nbf; strategy)

MFIE excitation for RectangularWaveguidePort.

Delegates to DeltaGapArrayPort implementation.
"""
function excitationVectorMFIE(
    port::RectangularWaveguidePort{FT, IT, DT},
    trianglesInfo::Vector{TriangleInfo{IT, FT}},
    nbf::Integer;
    strategy::Symbol = :efie_fallback
) where {FT<:Real, IT<:Integer, DT<:AbstractExcitationDistribution{FT}}

    # Delegate to DeltaGapArrayPort implementation
    gap_port = _to_delta_gap_array_port(port)
    return excitationVectorMFIE(gap_port, trianglesInfo, nbf; strategy = strategy)
end


# =============================================================================
# CFIE - Surface Ports
# =============================================================================

# -----------------------------------------------------------------------------
# DeltaGapPort - EFIE portion only for voltage sources
# -----------------------------------------------------------------------------

"""
    excitationVectorCFIE(port::DeltaGapPort, trianglesInfo, nbf; alpha, mfie_strategy)

Compute CFIE excitation vector for Delta-Gap port.

For voltage sources, we use EFIE excitation as the primary contribution.
The MFIE contribution can be optionally included using the specified strategy.

# Parameters
- `alpha::FT`: EFIE weight coefficient (default: 0.5)
- `mfie_strategy::Symbol`: MFIE excitation strategy (default: :efie_fallback)
"""
function excitationVectorCFIE(
    port::DeltaGapPort{FT, IT},
    trianglesInfo::Vector{TriangleInfo{IT, FT}},
    nbf::Integer;
    alpha::FT = FT(0.5),
    mfie_strategy::Symbol = :efie_fallback
) where {FT<:Real, IT<:Integer}

    # For voltage sources, use EFIE excitation
    V_efie = excitationVectorEFIE(port, trianglesInfo, nbf)

    if mfie_strategy == :efie_fallback
        # Use EFIE only for voltage sources in CFIE
        return V_efie
    else
        # Get MFIE excitation with specified strategy
        V_mfie = excitationVectorMFIE(port, trianglesInfo, nbf; strategy = mfie_strategy)

        # CFIE combination: V_CFIE = α × V_EFIE + (1-α) × η × V_MFIE
        one_minus_alpha = FT(1.0) - alpha
        return alpha .* V_efie .+ one_minus_alpha .* FT(η_0) .* V_mfie
    end
end


# -----------------------------------------------------------------------------
# CurrentProbe - Uses standard CFIE combination
# -----------------------------------------------------------------------------

"""
    excitationVectorCFIE(probe::CurrentProbe, trianglesInfo, nbf; alpha)

Compute CFIE excitation vector for current probe.

Uses standard CFIE combination:
```
V_CFIE = α × V_EFIE + (1-α) × η × V_MFIE
```
"""
function excitationVectorCFIE(
    probe::CurrentProbe{FT, IT},
    trianglesInfo::Vector{TriangleInfo{IT, FT}},
    nbf::Integer;
    alpha::FT = FT(0.5)
) where {FT<:Real, IT<:Integer}

    # Get EFIE and MFIE excitation vectors
    V_efie = excitationVectorEFIE(probe, trianglesInfo, nbf)
    V_mfie = excitationVectorMFIE(probe, trianglesInfo, nbf)

    # CFIE combination
    one_minus_alpha = FT(1.0) - alpha
    return alpha .* V_efie .+ one_minus_alpha .* FT(η_0) .* V_mfie
end


# -----------------------------------------------------------------------------
# DeltaGapArrayPort - EFIE only for voltage sources
# -----------------------------------------------------------------------------

"""
    excitationVectorCFIE(port::DeltaGapArrayPort, trianglesInfo, nbf; alpha)

Compute CFIE excitation vector for DeltaGapArrayPort.

For voltage sources, EFIE excitation is the natural choice.
"""
function excitationVectorCFIE(
    port::DeltaGapArrayPort{FT, IT, DT},
    trianglesInfo::Vector{TriangleInfo{IT, FT}},
    nbf::Integer;
    alpha::FT = FT(0.5),
    kwargs...
) where {FT<:Real, IT<:Integer, DT<:AbstractExcitationDistribution{FT}}

    # For voltage sources, EFIE is the natural choice
    return excitationVectorEFIE(port, trianglesInfo, nbf)
end


# -----------------------------------------------------------------------------
# RectangularWaveguidePort - Delegates to DeltaGapArrayPort
# -----------------------------------------------------------------------------

"""
    excitationVectorCFIE(port::RectangularWaveguidePort, trianglesInfo, nbf; alpha)

Compute CFIE excitation vector for RectangularWaveguidePort.

For voltage sources, EFIE excitation is the natural choice.
"""
function excitationVectorCFIE(
    port::RectangularWaveguidePort{FT, IT, DT},
    trianglesInfo::Vector{TriangleInfo{IT, FT}},
    nbf::Integer;
    alpha::FT = FT(0.5),
    mfie_strategy::Symbol = :efie_fallback
) where {FT<:Real, IT<:Integer, DT<:AbstractExcitationDistribution{FT}}

    # For voltage sources, EFIE is the natural choice
    return excitationVectorEFIE(port, trianglesInfo, nbf)
end
