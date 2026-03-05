## SanityChecks.jl  (MoM_Kernels level)
## Condition-number check for the MoM impedance matrix.
##
## This is the one check that lives here rather than MoM_Basics, because
## it references MoM_Kernels state (SimulationParams.ieT) and is naturally
## called from getImpedanceMatrix.
##
## NOTE: This file now uses the tiered O(n²) condition checker from
## CheapConditionCheck.jl instead of naive O(n³) SVD. The algorithm
## automatically selects:
##   - n ≤ 1000:   Exact 2-norm via SVD
##   - 1000 < n ≤ 5000: Exact 1-norm
##   - n > 5000:   LAPACK RCOND estimate (Hager-Higham)

# ─────────────────────────────────────────────────────────────────────────────
# Z-matrix condition number
# ─────────────────────────────────────────────────────────────────────────────

"""
    check_matrix_condition(Z; warn_threshold=1e6, error_threshold=1e10, 
                           force::Bool=false, check_quality::Bool=true) -> Bool

Check the condition number of the impedance matrix Z and emit diagnostics if ill-conditioned.

**Algorithm Selection (Tiered):**
- `n ≤ 1000`:   Exact 2-norm via SVD (cond(Z, 2)) — O(n³)
- `1000 < n ≤ 5000`: Exact 1-norm (cond(Z, 1)) — O(n²)
- `n > 5000`:   LAPACK RCOND estimate — O(n²)

This tiered approach ensures that condition number checks are computationally feasible
for matrices of any size, unlike the naive O(n³) SVD approach.

**Arguments:**
- `Z`: Impedance matrix (typically Complex{Float32} or Complex{Float64})
- `warn_threshold`: Condition number threshold for warning (default: 1e6)
- `error_threshold`: Condition number threshold for error-level warning (default: 1e10)
- `force`: If true, use exact 2-norm even for large matrices (expensive!)
- `check_quality`: If false, skip the check entirely and return true

**Returns:**
- `true` if condition number is acceptable (κ < error_threshold or check skipped)
- `false` if severely ill-conditioned (κ ≥ error_threshold)

**Example:**
```julia
Z = getImpedanceMatrix(geosInfo, nbf)
ok = check_matrix_condition(Z; check_quality=true)
!ok && @warn "Results may be unreliable due to ill-conditioning"
```

**References:**
- Hager (1984): "Condition estimates", SIAM J. Sci. Stat. Comput.
- Higham (1988): "FORTRAN codes for estimating the one-norm of a real or complex matrix"
- LAPACK: xGECON routine (reciprocal condition estimator)
"""
function check_matrix_condition(Z::AbstractMatrix;
                                warn_threshold::Real  = 1e6,
                                error_threshold::Real = 1e10,
                                force::Bool           = false,
                                check_quality::Bool   = true)
    !check_quality && return true
    
    n = size(Z, 1)
    
    # Handle force flag: use exact 2-norm regardless of size
    if force && n > 5000
        @info "[MoM quality] Forcing exact 2-norm condition check (expensive for n=$n)"
        κ = cond(Z, 2)
        method = "2-norm (SVD exact, forced)"
        
        # Evaluate and report
        if κ < warn_threshold
            @info "[MoM quality] Z-matrix condition OK ($(method)): κ ≈ $(round(κ; sigdigits=4))"
            return true
        elseif κ < error_threshold
            @warn """[MoM quality] Z-matrix ill-conditioned ($(method))
  κ(Z) ≈ $(round(κ; sigdigits=4)) (threshold = $warn_threshold)
  hint: refine mesh, avoid resonance, or apply iterative regularization"""
            return true
        else
            @warn """[MoM quality] Z-matrix severely ill-conditioned ($(method)) — results likely unreliable
  κ(Z) ≈ $(round(κ; sigdigits=4)) (threshold = $error_threshold)
  hint: refine mesh, apply low-frequency regularization, or use MLFMA partitioning"""
            return false
        end
    end
    
    # Use the tiered cheap checker (O(n²) for large matrices)
    # This handles all three tiers automatically
    FT = eltype(Z) <: Complex ? real(eltype(Z)) : eltype(Z)
    
    # Convert to Complex{FT} if needed for type stability
    if eltype(Z) <: Complex
        return CheapConditionCheck.check_matrix_condition_tiered(
            Z;
            warn_threshold = warn_threshold,
            error_threshold = error_threshold,
            check_quality = true,
            verbose = true
        )
    else
        # For real matrices, convert to complex (MoM matrices are typically complex)
        Z_complex = Complex{FT}.(Z)
        return CheapConditionCheck.check_matrix_condition_tiered(
            Z_complex;
            warn_threshold = warn_threshold,
            error_threshold = error_threshold,
            check_quality = true,
            verbose = true
        )
    end
end

# Legacy compatibility: keep the old signature but delegate to new implementation
# The 'force' parameter behavior is preserved for backward compatibility
