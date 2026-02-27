# Report: Z₀ (η₀) Factor Handling in CFIE Formulation

## Executive Summary

This report provides a comprehensive analysis of how the free-space impedance factor Z₀ (η₀) is used in the Combined Field Integral Equation (CFIE) implementation within this Method of Moments (MoM) codebase. The implementation follows the standard CFIE formulation where the magnetic field contribution is scaled by η₀ to ensure proper dimensional consistency and numerical accuracy.

---

## 1. Theoretical Background

### 1.1 CFIE Formulation

According to authoritative sources in computational electromagnetics, the Combined Field Integral Equation (CFIE) is defined as a linear combination of the Electric Field Integral Equation (EFIE) and the Magnetic Field Integral Equation (MFIE):

**Standard CFIE formulation:**
```
α·EFIE + (1-α)·MFIE
```

Where:
- **α** is the combination parameter (typically 0.5-0.6)
- **EFIE** is the Electric Field Integral Equation
- **MFIE** is the Magnetic Field Integral Equation

**Key theoretical insight from the literature:**
The normalized magnetic field **h** is defined as:
```
h = H/η  where η = √(μ/ε) is the wave impedance
```

This normalization means that when combining EFIE and MFIE, **dimensional consistency requires the MFIE portion to be scaled by η₀** when working with unnormalized magnetic fields.

### 1.2 Physical Constants

The codebase defines the following physical constants in `MoM_Basics.jl/src/ParametersSet.jl`:

| Constant | Symbol | Value | Description |
|----------|--------|-------|-------------|
| Speed of light | C₀ | 299792458.0 m/s | Vacuum speed of light |
| Vacuum permeability | μ₀ | 4π×10⁻⁷ H/m | Magnetic constant |
| Vacuum permittivity | ε₀ | 1/(μ₀·C₀²) | Electric constant |
| **Wave impedance of free space** | **η₀** | **√(μ₀/ε₀) ≈ 376.73 Ω** | **Key constant for CFIE** |
| Admittance | Y₀ | 1/η₀ | Free-space admittance |

---

## 2. Code Locations Where η₀ (Z₀) is Used

### 2.1 Constants Definition

**File:** `MoM_Basics.jl/src/ParametersSet.jl:45-51`

```julia
const η_0 = sqrt(μ_0/ε_0)           # Wave impedance of free space
const Y_0 = 1/η_0                    # Admittance of free space
const ηdiv16π = η_0/(16π)            # η₀/(16π) pre-computed constant
```

### 2.2 Pre-computed Constants

**File:** `MoM_Basics.jl/src/ParametersSet.jl:109-117`

```julia
JKη_0       = im*K_0*η_0             # j·k·η₀
Jη_0divK    = im*η_0/K_0             # j·η₀/k
JKηdiv16π   = 1im*K_0*η_0/16π        # j·k·η₀/(16π)
```

These pre-computed constants are used throughout the codebase for efficiency.

---

## 3. η₀ Usage in MFIE (Magnetic Field Integral Equation)

### 3.1 MFIE Impedance Matrix - Self Triangle

**File:** `MoM_Kernels.jl/src/ZmatAndVvec/MFIE/MFIERWGTri.jl:205-253`

For self-triangle (coincident field and source triangles):
```julia
# Line 214
ηdiv8S = η_0/(8*tri.area)

# Line 253
Ztemp *= lm*ln*ηdiv8S
```

**Analysis:** The impedance matrix element for self-interaction in MFIE is scaled by **η₀/(8S)**, where S is the triangle area.

### 3.2 MFIE Impedance Matrix - Different Triangles

**File:** `MoM_Kernels.jl/src/ZmatAndVvec/MFIE/MFIERWGTri.jl:8-100`

For non-coincident triangles:
```julia
# Lines 89-91
temp = lm*ln*ηdiv16π
Zmn *= temp
Znm *= temp
```

**Analysis:** The impedance matrix elements are scaled by **η₀/(16π)** multiplied by edge lengths.

### 3.3 MFIE Excitation Vector

**File:** `MoM_Kernels.jl/src/ZmatAndVvec/MFIE/MFIEExcitedVectors.jl:91`

```julia
Vtri[mi] *= η_0 * tri.edgel[mi]/2
```

**Analysis:** The excitation vector for MFIE is scaled by **η₀** times the edge length divided by 2.

---

## 4. η₀ Usage in CFIE (Combined Field Integral Equation)

### 4.1 CFIE Impedance Matrix - Self Triangle

**File:** `MoM_Kernels.jl/src/ZmatAndVvec/CFIE/CFIERWGTri.jl:235-346`

```julia
# Line 246
ηdiv8S = η_0/(8*tri.area)

# Lines 332-334
ZmnM *= (1-α)*lm*ln*ηdiv8S        # MFIE portion scaled by η₀
ZmnE *= α*lm*ln*JKηdiv16π         # EFIE portion scaled by j·k·η₀/(16π)
Zmn = ZmnE + ZmnM                 # Combined: α·EFIE + (1-α)·MFIE
```

**Key insight:**
- **EFIE portion:** Scaled by `α·jk·η₀/(16π)`
- **MFIE portion:** Scaled by `(1-α)·η₀/(8S)`
- The η₀ factor appears in **both** EFIE and MFIE portions

### 4.2 CFIE Impedance Matrix - Different Triangles

**File:** `MoM_Kernels.jl/src/ZmatAndVvec/CFIE/CFIERWGTri.jl:8-108`

```julia
# Lines 94-99
temp = (1-α)*lm*ln*ηdiv16π
ZmnM *= temp                       # MFIE portion
ZnmM *= temp
ZmnE *= α*lm*ln*JKηdiv16π         # EFIE portion
Zmn = ZmnE + ZmnM                  # Combined
```

### 4.3 CFIE Excitation Vector

**File:** `MoM_Kernels.jl/src/ZmatAndVvec/CFIE/CFIEExcitedVectors.jl:97-98`

```julia
# Line 97
Vtri[mi] = α*VE + (1-α)*η_0*VM

# Line 98
Vtri[mi] *= tri.edgel[mi]/2
```

**Where:**
- `VE` = Electric field contribution (ρ ⋅ E)
- `VM` = Magnetic field contribution (ρ ⋅ (n̂ × H))

**Critical observation:** The magnetic field contribution `VM` is multiplied by η₀ before the linear combination with EFIE.

---

## 5. η₀ Usage in EFIE (Electric Field Integral Equation)

### 5.1 EFIE Impedance Matrix

**File:** `MoM_Kernels.jl/src/ZmatAndVvec/EFIE/EFIERWGTri.jl:75`

```julia
Ztemp *= lm*ln*JKηdiv16π
```

**Analysis:** EFIE impedance elements are scaled by `jk·η₀/(16π)`, which includes η₀.

---

## 6. η₀ Usage in Field Calculations

### 6.1 Source Magnetic Field Calculation

**File:** `MoM_Basics.jl/src/Sources/Planewave.jl:76-78`

```julia
function sourceHfield(plw::PlaneWave, r)
    1/η_0*(plw.k̂ × sourceEfield(plw, r))
end
```

**Analysis:** The magnetic field is computed from the electric field using the plane wave relation:
```
H = (1/η₀)(k̂ × E)
```

This is the correct physical relationship for plane waves in free space.

### 6.2 Far-Field Calculation

**File:** `MoM_Kernels.jl/src/PostProcess/FarField.jl:37`

```julia
farEθϕ = (-Params.JK_0*η_0*div4π) .* Nθϕ
```

**Analysis:** The far-field calculation uses the factor `-jk·η₀/(4π)` which is the standard radiation formula.

---

## 7. Summary of η₀ Scaling Patterns

| Operation | EFIE Scaling | MFIE Scaling | CFIE Combined |
|-----------|--------------|--------------|---------------|
| Impedance (self) | - | η₀/(8S) | α·EFIE + (1-α)·MFIE |
| Impedance (different) | jk·η₀/(16π) | η₀/(16π) | α·EFIE + (1-α)·MFIE |
| Excitation | 1 | η₀ | α·VE + (1-α)·η₀·VM |

### Key Formula for CFIE:

**Impedance Matrix Element Zmn:**
```
Zmn = α · Zmn_EFIE + (1-α) · Zmn_MFIE
    = α · (jk·η₀/(16π)) · ⟨fₘ, T·fₙ⟩ + (1-α) · (η₀/(16π)) · ⟨fₘ, M·fₙ⟩
```

**Excitation Vector Element Vm:**
```
Vm = α · Vm_EFIE + (1-α) · Vm_MFIE
   = α · (lₘ/2) · ⟨ρₘ, E_inc⟩ + (1-α) · η₀ · (lₘ/2) · ⟨ρₘ, n̂ × H_inc⟩
```

---

## 8. Consistency Analysis

### 8.1 Theoretical Consistency

The implementation is **theoretically consistent** with the standard CFIE formulation:

1. **EFIE portion** uses `jk·η₀/(16π)` scaling
2. **MFIE portion** uses `η₀/(16π)` or `η₀/(8S)` scaling
3. **CFIE excitation** correctly scales the magnetic field term by η₀

### 8.2 Dimensional Analysis

- **EFIE terms** have units of Electric Field (V/m)
- **MFIE terms** (with η₀ scaling) also have units of Electric Field (V/m)
- The combination yields dimensionally consistent results

### 8.3 Code Consistency

The η₀ factor is applied consistently across:
- Direct solver implementations
- MLFMA (fast solver) implementations
- Near-field and far-field computations
- Both impedance matrix and excitation vector calculations

---

## 9. Conclusion

The codebase correctly implements the CFIE formulation with proper η₀ scaling:

1. **The MFIE contribution is consistently multiplied by η₀** in both the impedance matrix and excitation vector calculations
2. **The EFIE contribution contains η₀** within the `jk·η₀/(16π)` factor
3. **The CFIE linear combination** follows the standard formula: `α·EFIE + (1-α)·MFIE`
4. **Source field calculations** correctly use the plane wave relationship `H = (1/η₀)(k̂ × E)`

The implementation is mathematically sound and follows established practices in computational electromagnetics.

---

## References

### Academic Sources Consulted:
1. [Electromagnetic Integral Equations: Insights in Conditioning and Preconditioning](https://hal.science/hal-03472123/file/Electromagnetic_Integral_Equations_Insights_in_Conditioning_and_Preconditioning.pdf) - Adrian et al., IEEE Open Journal of Antennas and Propagation, 2021
2. [Altair Feko Documentation - CFIE](https://help.altair.com/feko/topics/feko/user_guide/solver_solution_methods/factor_cfie_feko_c.htm)
3. [Combined Field Integral Equation Based Theory](https://hub.hku.hk/bitstream/10722/247444/1/content.pdf)
4. [Fundamentals of Numerical and Asymptotic Methods](https://www.ursi.org/content/commissions/ComB/ComB_School/1st_EMTS_2013.pdf) - URSI

### Key Implementation Files:
- `MoM_Basics.jl/src/ParametersSet.jl` - Physical constants
- `MoM_Kernels.jl/src/ZmatAndVvec/CFIE/CFIERWGTri.jl` - CFIE impedance
- `MoM_Kernels.jl/src/ZmatAndVvec/CFIE/CFIEExcitedVectors.jl` - CFIE excitation
- `MoM_Kernels.jl/src/ZmatAndVvec/MFIE/MFIERWGTri.jl` - MFIE impedance
- `MoM_Kernels.jl/src/ZmatAndVvec/MFIE/MFIEExcitedVectors.jl` - MFIE excitation
- `MoM_Kernels.jl/src/ZmatAndVvec/EFIE/EFIERWGTri.jl` - EFIE impedance
- `MoM_Basics.jl/src/Sources/Planewave.jl` - Source field definitions
- `MoM_Kernels.jl/src/PostProcess/FarField.jl` - Far-field computation
