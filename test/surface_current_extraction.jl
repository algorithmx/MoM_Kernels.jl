using MoM_Basics
using MoM_Kernels
using Test
using StaticArrays
using LinearAlgebra

@testset "Surface Current Extraction" begin
    # Setup dummy data
    FT = Float64
    CT = Complex{FT}
    setPrecision!(FT)
    
    # Mock geometry - Using TriangleInfo from MoM_Basics
    # We need to construct real TriangleInfo objects because calSurfaceCurrents now expects specific types
    
    # Manually constructing TriangleInfo is complex because of internal structure.
    # It is easier to test with a small real mesh or rely on the type dispatch if we can construct minimal objects.
    
    tri1 = TriangleInfo{Int, FT}(1)
    tri1.center = SVec3D{FT}(0,0,0)
    tri1.area = 1.0
    tri1.inBfsID = MVec3D{Int}(1, 2, 0)
    # Set necessary fields for geoElectricJCal
    # geoElectricJCal uses vertices, edgel for calculation.
    tri1.vertices = MMatrix{3,3,FT,9}([0.0 1.0 0.0; 0.0 0.0 1.0; 0.0 0.0 0.0])
    tri1.edgel = MVec3D{FT}(1.0, 1.0, 1.0)
    
    tri2 = TriangleInfo{Int, FT}(2)
    tri2.center = SVec3D{FT}(1,0,0)
    tri2.area = 1.0
    tri2.inBfsID = MVec3D{Int}(1, 0, 0)
    tri2.vertices = MMatrix{3,3,FT,9}([1.0 2.0 1.0; 0.0 0.0 1.0; 0.0 0.0 0.0])
    tri2.edgel = MVec3D{FT}(1.0, 1.0, 1.0)

    geos = [tri1, tri2]
    
    # Create mock coefficients
    ICoeff = zeros(Complex{FT}, 2)
    ICoeff[1] = 1.0 + 0.5im
    ICoeff[2] = 0.2 + 0.1im
    
    # Initialize Parameters
    inputBasicParameters(frequency=1e8)
    updateVSBFTParams!(;sbfT = :RWG)
    
    # Test Data Structure
    # SurfaceCurrentData is deprecated/removed in basics, using FieldData now.
    # @test SurfaceCurrentData isa UnionAll
    
    # Test Calculation
    # Note: bfsInfo can be empty as it is ignored in new implementation
    data = calSurfaceCurrents(geos, [], ICoeff)
    
    @test data isa FieldData
    @test data.npoints == 2
    @test data.positions[1] ≈ SVec3D{FT}(0,0,0)
    @test data.positions[2] ≈ SVec3D{FT}(1,0,0)
    
    # Check that currents are computed (non-zero for non-zero coefficients)
    @test norm(data.fields[:J][1]) > 0 || norm(data.fields[:J][2]) > 0
    
    # Test Saving CSV
    filename_csv = "test_currents.csv"
    saveSurfaceCurrents(filename_csv, data)
    @test isfile(filename_csv)
    
    # Verify CSV content
    lines = readlines(filename_csv)
    @test length(lines) == 3 # Header + 2 rows
    @test startswith(lines[1], "rx,ry,rz")
    
    rm(filename_csv)
    
    # Test convenience wrapper (geosInfo, bfsInfo, ICoeff) -> save
    filename_conv = "test_currents_conv.csv"
    saveSurfaceCurrents(filename_conv, geos, [], ICoeff)
    @test isfile(filename_conv)
    rm(filename_conv)
end
