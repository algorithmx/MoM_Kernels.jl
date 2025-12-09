"""
    calSurfaceCurrents(geosInfo, bfsInfo, ICoeff)

Compute surface current density J at centroids using `geoElectricJCal`. Returns `FieldData`.
`bfsInfo` is provided for API compatibility but is not used (as `geosInfo` contains basis IDs).
"""
function calSurfaceCurrents(geosInfo, bfsInfo, ICoeff)
    # Dispatch based on type of geosInfo
    # Note: geoElectricJCal returns 3xN Matrix. FieldData needs Vector{SVec3D}.
    
    J_matrix = _compute_J_matrix(geosInfo, ICoeff)
    
    # Flatten geometry to get centers for FieldData
    geos_flat = _flatten_geos(geosInfo)
    npoints = length(geos_flat)
    
    FT = Precision.FT
    CT = Complex{FT}

    positions = Vector{SVec3D{FT}}(undef, npoints)
    J_vec     = Vector{SVec3D{CT}}(undef, npoints)
    
    for i in 1:npoints
        positions[i] = SVec3D{FT}(geos_flat[i].center)
        J_vec[i]     = SVec3D{CT}(J_matrix[1, i], J_matrix[2, i], J_matrix[3, i])
    end
    
    fd = FieldData{FT, CT}(npoints, positions)
    fd.fields[:J] = J_vec
    return fd
end

function _compute_J_matrix(geosInfo::AbstractVector{<:TriangleInfo}, ICoeff)
    return geoElectricJCal(ICoeff, geosInfo)
end

function _compute_J_matrix(geosInfo::AbstractVector{<:TetrahedraInfo}, ICoeff)
    return geoElectricJCal(ICoeff, geosInfo)
end

function _compute_J_matrix(geosInfo::AbstractVector{<:HexahedraInfo}, ICoeff)
    return geoElectricJCal(ICoeff, geosInfo)
end

function _compute_J_matrix(geosInfo::AbstractVector{<:AbstractVector}, ICoeff)
    # Mixed mesh: assumes [triangles, tetras/hexas]
    # geoElectricJCal can be called on each part
    # We need to concatenate results.
    # geosInfo is Vector of Vectors.
    Js = map(g -> geoElectricJCal(ICoeff, g), geosInfo)
    return reduce(hcat, Js)
end

function _compute_J_matrix(geosInfo::OffsetVector, ICoeff)
    # If using OffsetVector, it usually wraps a Vector. geoElectricJCal handles OffsetVector?
    # geoElectricJCal handles OffsetVector internally if passed directly.
    # But if geosInfo is a single OffsetVector, dispatch hits here.
    return geoElectricJCal(ICoeff, geosInfo)
end


function _flatten_geos(geosInfo::AbstractVector{<:VSCellType})
    return geosInfo
end

function _flatten_geos(geosInfo::AbstractVector{<:AbstractVector})
    return reduce(vcat, geosInfo)
end

"""
    saveSurfaceCurrents(filename::String, data::FieldData)

Save precomputed surface current field data to disk. This is a thin wrapper
around `saveFieldData` for API compatibility with existing code/tests.
"""
function saveSurfaceCurrents(filename::String, data::FieldData)
    saveFieldData(filename, data)
end

function saveSurfaceCurrents(filename::String, geosInfo, bfsInfo, ICoeff)
    data = calSurfaceCurrents(geosInfo, bfsInfo, ICoeff)
    saveFieldData(filename, data)
end
