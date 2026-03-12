# =============================================================================
# S-Parameter Conversion Functions
# =============================================================================
#
# This file provides conversion between different parameter representations:
# - Z-parameters (port impedance)
# - S-parameters (scattering)
# - T-parameters (transmission, for cascading)
#
# =============================================================================

"""
    z2s(Z_port::Matrix{CT}, Z0::FT) where {CT,FT}

Convert port impedance matrix to S-parameters.

# Formula
S = (Z - Z₀⋅I) ⋅ (Z + Z₀⋅I)⁻¹

# Arguments
- `Z_port::Matrix{Complex{FT}}`: Port impedance matrix (N×N)
- `Z0::Real`: Reference impedance (default: 50 Ω)

# Returns
- `Matrix{Complex{FT}}`: S-parameter matrix (N×N)

# Example
```julia
Z_port = [50+10im  5-2im;
           5-2im  50+10im]
S = z2s(Z_port, 50.0)
```
"""
function z2s(Z_port::Matrix{CT}, Z0::FT) where {CT<:Complex, FT<:Real}
    n = size(Z_port, 1)
    @assert size(Z_port, 2) == n "Z_port must be square"
    
    # Z0 * Identity matrix
    Z0_mat = Z0 * Matrix{CT}(I, n, n)
    
    # Compute S = (Z - Z0*I) / (Z + Z0*I)
    numerator = Z_port - Z0_mat
    denominator = Z_port + Z0_mat
    
    # Solve for S: S ⋅ (Z + Z0*I) = (Z - Z0*I)
    # S = (Z - Z0*I) ⋅ (Z + Z0*I)⁻¹
    S = numerator / denominator
    
    return S
end

"""
    s2z(S::Matrix{CT}, Z0::FT) where {CT,FT}

Convert S-parameters to port impedance matrix.

# Formula
Z = Z₀ ⋅ (I + S) ⋅ (I - S)⁻¹

# Arguments
- `S::Matrix{Complex{FT}}`: S-parameter matrix (N×N)
- `Z0::Real`: Reference impedance (default: 50 Ω)

# Returns
- `Matrix{Complex{FT}}`: Port impedance matrix (N×N)

# Example
```julia
S = [0.5+0.1im  0.1-0.05im;
     0.1-0.05im  0.5+0.1im]
Z_port = s2z(S, 50.0)
```
"""
function s2z(S::Matrix{CT}, Z0::FT) where {CT<:Complex, FT<:Real}
    n = size(S, 1)
    @assert size(S, 2) == n "S must be square"
    
    # Identity matrix
    I_mat = Matrix{CT}(I, n, n)
    
    # Compute Z = Z0 * (I + S) / (I - S)
    numerator = I_mat + S
    denominator = I_mat - S
    
    # Z = Z0 * (I + S) ⋅ (I - S)⁻¹
    Z_port = Z0 * (numerator / denominator)
    
    return Z_port
end

"""
    s2t(S::Matrix{CT}) where CT

Convert 2-port S-parameters to T-parameters (transmission matrix).

T-parameters are useful for cascading networks:
T_total = T₁ ⋅ T₂ ⋅ T₃ ...

# Formula (2-port only)
```
T = 1/S₂₁ ⋅ [ -det(S)    S₁₁
              -S₂₂       1   ]
```

# Arguments
- `S::Matrix{Complex{FT}}`: 2×2 S-parameter matrix

# Returns
- `Matrix{Complex{FT}}`: 2×2 T-parameter matrix

# Note
T-parameters are only defined for 2-port networks. For N-port networks
with N > 2, use different cascading techniques.

# Example
```julia
S = [0.1+0.05im  0.9-0.1im;
     0.9-0.1im   0.1+0.05im]
T = s2t(S)
```
"""
function s2t(S::Matrix{CT}) where {CT<:Complex}
    @assert size(S) == (2, 2) "s2t() only defined for 2-port networks"
    
    S11, S12, S21, S22 = S[1,1], S[1,2], S[2,1], S[2,2]
    
    # Check for singular case
    if abs(S21) < eps(real(CT))
        error("S₂₁ ≈ 0, T-parameters undefined")
    end
    
    # T-matrix elements
    det_S = S11*S22 - S12*S21
    
    T = Matrix{CT}(undef, 2, 2)
    T[1,1] = -det_S / S21
    T[1,2] = S11 / S21
    T[2,1] = -S22 / S21
    T[2,2] = 1 / S21
    
    return T
end

"""
    t2s(T::Matrix{CT}) where CT

Convert 2-port T-parameters to S-parameters.

# Formula (2-port only)
```
S = 1/T₂₂ ⋅ [ T₁₂       det(T)
              1        -T₂₁    ]
```

# Arguments
- `T::Matrix{Complex{FT}}`: 2×2 T-parameter matrix

# Returns
- `Matrix{Complex{FT}}`: 2×2 S-parameter matrix
"""
function t2s(T::Matrix{CT}) where {CT<:Complex}
    @assert size(T) == (2, 2) "t2s() only defined for 2-port networks"
    
    T11, T12, T21, T22 = T[1,1], T[1,2], T[2,1], T[2,2]
    
    # Check for singular case
    if abs(T22) < eps(real(CT))
        error("T₂₂ ≈ 0, S-parameters undefined")
    end
    
    # S-matrix elements
    det_T = T11*T22 - T12*T21
    
    S = Matrix{CT}(undef, 2, 2)
    S[1,1] = T12 / T22
    S[1,2] = det_T / T22
    S[2,1] = 1 / T22
    S[2,2] = -T21 / T22
    
    return S
end

"""
    cascadeT(T1::Matrix{CT}, T2::Matrix{CT}, Trest::Matrix{CT}...) where CT

Cascade T-parameters: T_total = T₁ ⋅ T₂ ⋅ T₃ ...

# Arguments
- `T1, T2, ...`: T-parameter matrices (2×2 each)

# Returns
- `Matrix{Complex{FT}}`: Cascaded T-parameter matrix

# Example
```julia
T_total = cascadeT(T_network1, T_network2, T_network3)
S_total = t2s(T_total)
```
"""
function cascadeT(T1::Matrix{CT}, T2::Matrix{CT}, Trest::Matrix{CT}...) where {CT<:Complex}
    # Start with T1 ⋅ T2
    T_result = T1 * T2
    
    # Multiply remaining matrices
    for T in Trest
        T_result = T_result * T
    end
    
    return T_result
end

"""
    cascadeS(S1::Matrix{CT}, S2::Matrix{CT}, Z0::Real=50.0) where CT

Cascade two 2-port S-parameter networks.

# Arguments
- `S1::Matrix`: S-parameters of first network
- `S2::Matrix`: S-parameters of second network  
- `Z0::Real`: Reference impedance (must be same for both)

# Returns
- `Matrix{Complex{FT}}`: S-parameters of cascaded network

# Note
This converts to T-parameters for the cascade operation:
S_cascade = t2s(s2t(S1) ⋅ s2t(S2))
"""
function cascadeS(S1::Matrix{CT}, S2::Matrix{CT}, Z0::Real=50.0) where {CT<:Complex}
    # Convert to T, cascade, convert back to S
    T1 = s2t(S1)
    T2 = s2t(S2)
    T_cascade = T1 * T2
    return t2s(T_cascade)
end

# =============================================================================
# Exports
# =============================================================================

export z2s, s2z, s2t, t2s, cascadeT, cascadeS
