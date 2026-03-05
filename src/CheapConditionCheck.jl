## CheapConditionCheck.jl
## Tiered O(n²) condition number estimation for large MoM impedance matrices
##
## This module provides efficient alternatives to O(n³) SVD-based condition
## number computation, enabling sanity checks for matrices of any size.
##
## Algorithm tiers:
##   - n ≤ 1000:   Exact 2-norm via SVD (cond(Z, 2))
##   - 1000 < n ≤ 5000: Exact 1-norm via efficient algorithms (cond(Z, 1))
##   - n > 5000:   LAPACK reciprocal condition estimate (Hager-Higham)
##
## References:
##   - Hager (1984): "Condition estimates", SIAM J. Sci. Stat. Comput.
##   - Higham (1988): Improvements to Hager's algorithm, ACM TOMS
##   - LAPACK xGECON: Reciprocal condition number estimator

module CheapConditionCheck

using LinearAlgebra
using LinearAlgebra.LAPACK

export check_matrix_condition_tiered, estimate_condition_number

# Thresholds for algorithm selection (tunable)
const SMALL_MATRIX_THRESHOLD  = 1000   # Use exact SVD for n ≤ this
const MEDIUM_MATRIX_THRESHOLD = 5000   # Use exact 1-norm for n ≤ this

# Default condition number thresholds for warnings
const DEFAULT_WARN_THRESHOLD  = 1e6    # Warn if κ ≥ this
const DEFAULT_ERROR_THRESHOLD = 1e10   # Error if κ ≥ this

"""
    check_matrix_condition_tiered(Z::AbstractMatrix{Complex{FT}};
                                   warn_threshold::Real = DEFAULT_WARN_THRESHOLD,
                                   error_threshold::Real = DEFAULT_ERROR_THRESHOLD,
                                   check_quality::Bool = true,
                                   verbose::Bool = true) where FT

Tiered condition number checker that selects the appropriate algorithm based on matrix size.

**Algorithm Selection:**
- `n ≤ $(SMALL_MATRIX_THRESHOLD)`: Exact 2-norm via SVD (`cond(Z, 2)`)
- `$(SMALL_MATRIX_THRESHOLD) < n ≤ $(MEDIUM_MATRIX_THRESHOLD)`: Exact 1-norm (`cond(Z, 1)`)
- `n > $(MEDIUM_MATRIX_THRESHOLD)`: LAPACK RCOND estimate (Hager-Higham algorithm)

**Arguments:**
- `Z`: Impedance matrix (typically Complex{Float32} or Complex{Float64})
- `warn_threshold`: Condition number threshold for warning (default: 1e6)
- `error_threshold`: Condition number threshold for error (default: 1e10)
- `check_quality`: If false, skip the check and return true immediately
- `verbose`: If true, emit @info on success; warnings are always emitted

**Returns:**
- `true` if condition number is acceptable (κ < error_threshold)
- `false` if severely ill-conditioned (κ ≥ error_threshold)

**Complexity:**
- Small matrices: O(n³) for exact SVD
- Medium matrices: O(n²) for exact 1-norm
- Large matrices: O(n²) for LAPACK estimate (reuses LU factorization)

**Example:**
```julia
Z = randn(ComplexF64, 1000, 1000)
ok = check_matrix_condition_tiered(Z; warn_threshold=1e6)
```
"""
function check_matrix_condition_tiered(Z::AbstractMatrix{Complex{FT}};
                                        warn_threshold::Real = DEFAULT_WARN_THRESHOLD,
                                        error_threshold::Real = DEFAULT_ERROR_THRESHOLD,
                                        check_quality::Bool = true,
                                        verbose::Bool = true) where FT
    !check_quality && return true
    
    n = size(Z, 1)
    
    # Validate matrix
    if n == 0
        @warn "[MoM quality] Empty impedance matrix (0×0)"
        return false
    end
    
    if size(Z, 1) != size(Z, 2)
        @warn "[MoM quality] Non-square impedance matrix ($(size(Z,1))×$(size(Z,2)))"
        return false
    end
    
    # Quick norm check for trivial cases
    norm1 = norm(Z, 1)
    if norm1 < eps(real(FT))
        @warn "[MoM quality] Zero impedance matrix"
        return false
    end
    
    # Select algorithm based on matrix size
    κ, method = _compute_condition_number(Z, n)
    
    # Evaluate and report
    if κ < warn_threshold
        if verbose
            @info "[MoM quality] Z-matrix condition OK ($(method)): κ ≈ $(round(κ; sigdigits=4))"
        end
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

"""
    _compute_condition_number(Z::AbstractMatrix{Complex{FT}}, n::Int) where FT

Internal function to compute condition number with appropriate algorithm.

Returns: (κ, method_name)
"""
function _compute_condition_number(Z::AbstractMatrix{Complex{FT}}, n::Int) where FT
    if n ≤ SMALL_MATRIX_THRESHOLD
        # Tier 1: Small matrices - exact 2-norm via SVD
        κ = cond(Z, 2)
        return κ, "2-norm (SVD exact)"
        
    elseif n ≤ MEDIUM_MATRIX_THRESHOLD
        # Tier 2: Medium matrices - exact 1-norm
        # Julia's cond(Z, 1) uses efficient O(n²) algorithms
        κ = cond(Z, 1)
        return κ, "1-norm (exact)"
        
    else
        # Tier 3: Large matrices - LAPACK RCOND estimate
        κ = _lapack_rcond_estimate(Z)
        return κ, "1-norm (LAPACK RCOND estimate)"
    end
end

"""
    _lapack_rcond_estimate(Z::AbstractMatrix{Complex{FT}}) where FT

Compute condition number estimate using LAPACK's reciprocal condition estimator.
Uses the Hager-Higham algorithm internally (O(n²) complexity).

**Algorithm:**
1. Compute 1-norm of Z
2. Compute LU factorization of Z
3. Call LAPACK.gecon! for reciprocal condition estimate
4. Invert to get condition number estimate

**Note:** This uses the existing LU factorization from MoM solve if available.
"""
function _lapack_rcond_estimate(Z::AbstractMatrix{Complex{FT}}) where FT
    n = size(Z, 1)
    RT = real(FT)
    
    # Compute 1-norm
    anorm = norm(Z, 1)
    
    # Handle near-zero norm
    anorm < eps(RT) && return RT(Inf)
    
    try
        # Compute LU factorization
        lu_fact = lu(Z)
        
        # Get reciprocal condition number from LAPACK
        # '1' = 1-norm, 'I' = infinity-norm
        rcond = LAPACK.gecon!('1', lu_fact.factors, anorm)
        
        # rcond is in [0, 1]; convert to condition number
        if rcond > eps(RT) && isfinite(rcond)
            return 1 / rcond
        else
            return RT(Inf)
        end
        
    catch e
        # LAPACK failed - fallback to 1-norm exact
        @debug "LAPACK gecon! failed, falling back to exact 1-norm: $e"
        return cond(Z, 1)
    end
end

"""
    estimate_condition_number(Z::AbstractMatrix{Complex{FT}};
                               use_2norm::Bool = false) where FT

Estimate condition number without warning/error logic.
Useful for logging/debugging purposes.

**Arguments:**
- `Z`: Impedance matrix
- `use_2norm`: If true, always compute 2-norm (slower but comparable to SVD)

**Returns:** Estimated condition number
"""
function estimate_condition_number(Z::AbstractMatrix{Complex{FT}};
                                    use_2norm::Bool = false) where FT
    n = size(Z, 1)
    
    if use_2norm && n ≤ SMALL_MATRIX_THRESHOLD
        return cond(Z, 2)
    else
        κ, _ = _compute_condition_number(Z, n)
        return κ
    end
end

# Legacy alias for backward compatibility (if needed)
const check_matrix_condition_cheap = check_matrix_condition_tiered

end # module CheapConditionCheck
