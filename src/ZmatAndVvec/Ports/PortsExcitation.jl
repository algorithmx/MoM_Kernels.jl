# =============================================================================
# PortsExcitation.jl - Port excitation vector implementations
# =============================================================================
#
# This directory contains port-specific excitation vector calculations that
# use Julia's multiple dispatch to override the general ExcitingSources methods.
#
# Architecture:
#   - PortType <: ExcitingSource (ports ARE sources)
#   - Julia's dispatch automatically prefers PortType methods over ExcitingSources
#
# File Organization:
#   - SurfacePortExcitation.jl: Delta-Gap, Current Probe, Array Port (RWG basis)
#   - VolumePortExcitation.jl:  Wave Port (SWG/PWC basis)
#
# IE Compatibility Matrix:
#   | Port Type               | EFIE | MFIE | CFIE |
#   |-------------------------|------|------|------|
#   | DeltaGapPort            | ✓    | EFIE | ✓    |
#   | CurrentProbe            | ✓    | ✓    | ✓    |
#   | DeltaGapArrayPort       | ✓    | EFIE | ✓    |
#   | RectangularEdgePort| ✓    | EFIE | ✓    |
#   | WavePort (planned)      | ✓    | N/A  | N/A  |
#
# Key Design Principles:
#   1. Same interface: Uses excitationVectorEFIE/MFIE/CFIE - no new API
#   2. Automatic dispatch: Julia selects port methods via type constraints
#   3. No entry point changes: getExcitationVector works unchanged
#
# Extending with New Port Types:
#   1. Define struct in MoM_Basics.jl/src/Sources/ subtyping PortType
#   2. Add excitation methods in this directory
#   3. No changes needed to entry points or assembly logic
#
# =============================================================================

# Surface port excitation (RWG basis functions)
include("SurfacePortExcitation.jl")

# Volume port excitation (SWG/PWC basis functions)
include("VolumePortExcitation.jl")
