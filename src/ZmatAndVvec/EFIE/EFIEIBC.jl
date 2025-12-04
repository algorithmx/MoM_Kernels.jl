"""
Calculates the IBC correction term for a single triangle.
This is the Gram matrix of the RWG basis functions on the triangle, scaled by Zs.
"""
function IBCOnTri(tri::TriangleInfo{IT, FT}) where {IT, FT}
    CT = Complex{FT}
    Z_ibc = zeros(CT, 3, 3)
    
    if iszero(tri.Zs)
        return Z_ibc
    end

    rgt = getGQPTri(tri)
    # Factor: Zs * Area / (4 * Area^2) = Zs / (4 * Area)
    # Integral = Sum( (lm/2A * rho_m) . (ln/2A * rho_n) * Zs ) * Area
    #          = (lm * ln * Zs / (4 * A^2)) * Sum( rho_m . rho_n * weight ) * Area
    #          = (lm * ln * Zs / (4 * A)) * Sum( rho_m . rho_n * weight )
    factor = tri.Zs / (4 * tri.area)

    for g in 1:GQPNTri
        # Quadrature point vector
        rg = view(rgt, :, g)
        weight = TriGQInfo.weight[g]
        
        for ni in 1:3
            freeVn = view(tri.vertices, :, ni)
            ρn = rg - freeVn
            
            for mi in 1:3
                freeVm = view(tri.vertices, :, mi)
                ρm = rg - freeVm
                
                Z_ibc[mi, ni] += (ρm ⋅ ρn) * weight
            end
        end
    end

    # Apply constant factors
    for ni in 1:3
        ln = tri.edgel[ni]
        for mi in 1:3
            lm = tri.edgel[mi]
            Z_ibc[mi, ni] *= factor * lm * ln
        end
    end

    return Z_ibc
end

"""
Computes EFIE Impedance Matrix with IBC corrections.
Z_total = Z_PEC - Z_IBC
"""
function impedancemat4EFIE4IBC(trianglesInfo::Vector{TriangleInfo{IT, FT}}, nrwg::Integer, bfT::Type{BFT}) where {IT, FT, BFT<:RWG}
    # 1. Compute standard PEC matrix
    Zmat = impedancemat4EFIE4PEC(trianglesInfo, nrwg, bfT)
    
    lockZ = SpinLock()

    # 2. Apply IBC Correction
    Threads.@threads for tri in trianglesInfo
        # Check if Zs is zero to avoid unnecessary computation
        if !iszero(tri.Zs)
            Z_correction = IBCOnTri(tri)
            
            lock(lockZ)
            for ni in 1:3, mi in 1:3
                m = tri.inBfsID[mi]
                n = tri.inBfsID[ni]
                
                # If edge is boundary (0), skip
                (m == 0 || n == 0) && continue
                
                # Subtract correction: Z_total = Z_PEC - Z_IBC
                Zmat[m, n] -= Z_correction[mi, ni]
            end
            unlock(lockZ)
        end
    end
    
    return Zmat
end
