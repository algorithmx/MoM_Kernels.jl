using MoM_Kernels
using MoM_Basics
using Test
using StaticArrays
using LinearAlgebra

# Import internal (non-exported) kernel routines needed for these tests
import MoM_Kernels: IBCOnTri, impedancemat4EFIE4PEC, impedancemat4EFIE4IBC

@testset "IBC Kernel Implementation" begin
    # Helper to create a dummy triangle with valid geometry
    function create_test_triangle()
        tri = TriangleInfo{Int64, Float64}(1)
        # Define a simple right-angled triangle in XY plane
        # Area = 0.5 * 1.0 * 1.0 = 0.5
        tri.vertices[:, 1] = [0.0, 0.0, 0.0]
        tri.vertices[:, 2] = [1.0, 0.0, 0.0]
        tri.vertices[:, 3] = [0.0, 1.0, 0.0]
        
        # Manually set geometry parameters usually done by setTriParam!
        tri.area = 0.5
        tri.center = [1/3, 1/3, 0.0]
        tri.edgel = [1.0, sqrt(2.0), 1.0] # lengths of edges opposite to vertices 1, 2, 3
        # Edge 1 (opp to v1): v2->v3. len=sqrt(2)
        # Edge 2 (opp to v2): v3->v1. len=1
        # Edge 3 (opp to v3): v1->v2. len=1
        # The struct stores edgel[i] as length of edge opposite to vertex i? 
        # Let's check TriangleInfo definition in MoM_Basics: 
        # "rational arrangement... three basis functions free ends are the three vertices"
        # Usually edge i is opposite to vertex i.
        # Let's assume standard MoM setup.
        
        # Set dummy Basis IDs so it's not treated as boundary
        tri.inBfsID = [1, 2, 3]
        
        return tri
    end

    @testset "IBCOnTri Calculation" begin
        tri = create_test_triangle()
        
        # Case 1: Zs = 0
        tri.Zs = 0.0 + 0.0im
        Z_ibc = IBCOnTri(tri)
        @test all(iszero, Z_ibc)

        # Case 2: Non-zero Zs
        Zs_val = 100.0 + 50.0im
        tri.Zs = Zs_val
        Z_ibc = IBCOnTri(tri)
        
        # Property checks
        # 1. Symmetry: Matrix should be symmetric (Gram matrix of real functions)
        @test Z_ibc[1, 2] ≈ Z_ibc[2, 1]
        @test Z_ibc[1, 3] ≈ Z_ibc[3, 1]
        @test Z_ibc[2, 3] ≈ Z_ibc[3, 2]
        
        # 2. Scaling: If we double Zs, matrix should double
        tri.Zs = 2.0 * Zs_val
        Z_ibc_2 = IBCOnTri(tri)
        @test Z_ibc_2 ≈ 2.0 .* Z_ibc
    end

    @testset "impedancemat4EFIE4IBC Integration" begin
        # Setup a minimal system with 1 triangle (physically impossible for closed surface but ok for unit test of matrix assembly)
        tri = create_test_triangle()
        tris = [tri]
        nrwg = 3 # max ID in tri.inBfsID
        
        # 1. Compute PEC matrix (Standard)
        Z_pec = impedancemat4EFIE4PEC(tris, nrwg, RWG)
        
        # 2. Compute IBC matrix with Zs=0 (Should match PEC)
        tri.Zs = 0.0 + 0.0im
        Z_ibc_zero = impedancemat4EFIE4IBC(tris, nrwg, RWG)
        @test Z_ibc_zero ≈ Z_pec
        
        # 3. Compute IBC matrix with Zs != 0
        # Should equal PEC - Correction
        Zs_val = 377.0 + 0.0im
        tri.Zs = Zs_val
        
        # Manually compute correction
        Z_correction = IBCOnTri(tri)
        
        Z_ibc_full = impedancemat4EFIE4IBC(tris, nrwg, RWG)
        
        # Verify: Z_total[m,n] = Z_pec[m,n] - Z_correction[mi, ni]
        # Since we have only 1 triangle, indices match directly 1-to-1 for this simple case
        for i in 1:3, j in 1:3
            if i == j 
                 # Diagonal term in global matrix might involve multiple triangles in real mesh, 
                 # but here just 1 triangle.
            end
            # Check the subtraction logic
            @test Z_ibc_full[i, j] ≈ Z_pec[i, j] - Z_correction[i, j]
        end
    end

end
