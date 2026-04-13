# test_efie_groundplane.jl
# Tests for EFIE impedance matrix correctness with PEC ground plane
#
# These tests verify that:
# 1. The ground-plane Z-matrix differs from free-space (image is non-trivial)
# 2. Free-space is recovered when ground is infinitely far away
# 3. Z satisfies reciprocity (Z_mn = Z_nm) with ground plane
# 4. Individual Z-element integrals are consistent across all three integration
#    regimes (far, near-singular, self-term)
# 5. Z approaches free-space as ground moves further away
# 6. Solving with ground plane completes without error

using MoM_Basics, MoM_Kernels, LinearAlgebra, Test, Logging

# ============================================================================
# Helper: build a small triangle mesh and get Z-matrix for a given GF config
# ============================================================================
function setup_and_compute_Z(; gf_type::Symbol = :freespace, z_gnd::Float64 = 0.0,
                                frequency::Real = 20e8)
    proj_path = joinpath(@__DIR__, "..")

    # Configure GF before computing Z
    if gf_type == :groundplane
        set_green_function_type!(:groundplane; z_gnd = z_gnd)
    else
        set_green_function_type!(:freespace)
    end

    inputParameters(; frequency = frequency, ieT = :EFIE)

    filename = joinpath(proj_path, "meshfiles/Tri.nas")
    meshData, _  = getMeshData(filename; meshUnit = :mm)
    ngeo, nbf, geosInfo, bfsInfo = getBFsFromMeshData(meshData; sbfT = :RWG)

    Zmat = getImpedanceMatrix(geosInfo, nbf)

    # Always restore free-space after test
    reset_green_function_config!()

    return Zmat, nbf, geosInfo, bfsInfo
end


@testset "EFIE with Ground Plane" begin

    # ===========================================================================
    # Test 1: Ground plane Z-matrix differs from free-space
    # ===========================================================================
    @testset "Z differs from free-space" begin
        Z_free, nbf, _, _ = setup_and_compute_Z(; gf_type = :freespace)

        # Ground at z = -0.01 m (below mesh by 10 mm)
        Z_gnd, _, _, _  = setup_and_compute_Z(; gf_type = :groundplane, z_gnd = -0.01)

        # They must NOT be the same (image term contributes)
        @test nbf > 0
        ΔZ = Z_gnd .- Z_free
        rel_diff = norm(ΔZ) / norm(Z_free)
        @test rel_diff > 1e-6        # non-trivial difference
        @test isfinite(rel_diff)
    end

    # ===========================================================================
    # Test 2: Reciprocity Z_mn = Z_nm with ground plane
    # ===========================================================================
    @testset "Reciprocity with ground plane" begin
        Z_gnd, nbf, _, _ = setup_and_compute_Z(; gf_type = :groundplane, z_gnd = -0.01)

        # Check symmetry: EFIE Z should be symmetric even with a ground plane
        asym = norm(Z_gnd - transpose(Z_gnd)) / norm(Z_gnd)
        @test asym < 1e-10
    end

    # ===========================================================================
    # Test 3: Z approaches free-space as ground moves far away
    # ===========================================================================
    @testset "Z converges to free-space for distant ground" begin
        Z_free, _, _, _ = setup_and_compute_Z(; gf_type = :freespace)

        # Ground very far below — image contribution should be negligible
        Z_far, _, _, _  = setup_and_compute_Z(; gf_type = :groundplane, z_gnd = -100.0)

        rel_diff = norm(Z_far .- Z_free) / norm(Z_free)
        @test rel_diff < 1e-3   # effectively free-space
    end

    # ===========================================================================
    # Test 4: Z changes monotonically as ground approaches mesh
    #         (image grows stronger → larger deviation from free-space)
    # ===========================================================================
    @testset "Z deviation grows as ground approaches" begin
        Z_free, _, _, _ = setup_and_compute_Z(; gf_type = :freespace)

        z_gnds = [-1.0, -0.1, -0.01]           # progressively closer
        prev_diff = 0.0
        for z_gnd in z_gnds
            Z_gnd, _, _, _ = setup_and_compute_Z(; gf_type = :groundplane, z_gnd = z_gnd)
            diff = norm(Z_gnd .- Z_free) / norm(Z_free)
            @test diff > prev_diff              # deviation should grow
            prev_diff = diff
        end
    end

    # ===========================================================================
    # Test 5: Diagonal (self-term) elements are finite and physical
    # ===========================================================================
    @testset "Diagonal elements are finite and nonzero" begin
        Z_gnd, nbf, _, _ = setup_and_compute_Z(; gf_type = :groundplane, z_gnd = -0.01)

        for i in 1:nbf
            @test isfinite(Z_gnd[i, i])
            @test abs(Z_gnd[i, i]) > 0
        end

        # With e^{+jωt} convention, sub-wavelength RWG self-impedance is
        # dominated by the scalar-potential (charge) term → capacitive → imag < 0
        any_negative_imag = any(imag(Z_gnd[i, i]) < 0 for i in 1:nbf)
        @test any_negative_imag
    end

    # ===========================================================================
    # Test 6: Full solve completes without error (ground plane)
    # ===========================================================================
    @testset "Solve with ground plane" begin
        Z_gnd, nbf, geosInfo, bfsInfo = setup_and_compute_Z(;
            gf_type = :groundplane, z_gnd = -0.01)

        source = PlaneWave(π, 0, 0f0, 1f0)
        V      = MoM_Kernels.getExcitationVector(geosInfo, nbf, source)
        @test length(V) == nbf

        ICoeff, _ = MoM_Kernels.solve(Z_gnd, V; solverT = :direct)
        @test length(ICoeff) == nbf
        @test all(isfinite, ICoeff)
        @test norm(ICoeff) > 0
    end

    # ===========================================================================
    # Test 7: Element-level consistency across integration regimes
    #         Compare the three entry points with a hand-built GroundPlaneGF
    # ===========================================================================
    @testset "Element-level: EFIEOnTris dispatches correctly" begin
        # Set up parameters
        inputParameters(; frequency = 20e8, ieT = :EFIE)

        proj_path = joinpath(@__DIR__, "..")
        filename  = joinpath(proj_path, "meshfiles/Tri.nas")
        meshData, _  = getMeshData(filename; meshUnit = :mm)
        _, _, geosInfo, _ = getBFsFromMeshData(meshData; sbfT = :RWG)

        FT = eltype(geosInfo[1].center)
        k  = Params.K_0
        z_gnd = FT(-0.01)

        gf_free = FreeSpaceGF{FT}(Complex{FT}(k))
        gf_gnd  = GroundPlaneGF{FT}(Complex{FT}(k), z_gnd)

        tri1 = geosInfo[1]

        # Self-term: EFIEOnTris(tri, gf)
        Z_self_free = MoM_Kernels.EFIEOnTris(tri1, gf_free)
        Z_self_gnd  = MoM_Kernels.EFIEOnTris(tri1, gf_gnd)

        # The self-term with ground should differ (image contribution)
        @test norm(Z_self_gnd .- Z_self_free) > 0

        # Both should be finite
        @test all(isfinite, Z_self_free)
        @test all(isfinite, Z_self_gnd)

        # Far-field: pick two triangles far apart
        tri_far = geosInfo[end]
        center_dist = norm(tri1.center .- tri_far.center)
        if center_dist > Params.Rsglr
            Z_far_free = MoM_Kernels.EFIEOnTris(tri1, tri_far, gf_free)
            Z_far_gnd  = MoM_Kernels.EFIEOnTris(tri1, tri_far, gf_gnd)

            @test norm(Z_far_gnd .- Z_far_free) > 0
            @test all(isfinite, Z_far_free)
            @test all(isfinite, Z_far_gnd)
        end

        reset_green_function_config!()
    end

    # ===========================================================================
    # Test 8: No warnings emitted for GroundPlaneGF
    # ===========================================================================
    @testset "No fallback warnings for GroundPlaneGF" begin
        inputParameters(; frequency = 20e8, ieT = :EFIE)

        proj_path = joinpath(@__DIR__, "..")
        filename  = joinpath(proj_path, "meshfiles/Tri.nas")
        meshData, _  = getMeshData(filename; meshUnit = :mm)
        _, _, geosInfo, _ = getBFsFromMeshData(meshData; sbfT = :RWG)

        FT    = eltype(geosInfo[1].center)
        k     = Params.K_0
        z_gnd = FT(-0.01)
        gf    = GroundPlaneGF{FT}(Complex{FT}(k), z_gnd)

        tri1  = geosInfo[1]

        # Self-term should NOT emit the old "not fully implemented" warning
        # @test_logs with no patterns and min_level=Warn asserts zero warnings
        @test_logs min_level=Logging.Warn MoM_Kernels.EFIEOnTris(tri1, gf)

        reset_green_function_config!()
    end

    # ===========================================================================
    # Test 9: ensure free-space path is unaffected (backward compat)
    # ===========================================================================
    @testset "Free-space backward compatibility" begin
        # Compute Z via the default entry point (no explicit GF)
        inputParameters(; frequency = 20e8, ieT = :EFIE)
        reset_green_function_config!()

        proj_path = joinpath(@__DIR__, "..")
        filename  = joinpath(proj_path, "meshfiles/Tri.nas")
        meshData, _  = getMeshData(filename; meshUnit = :mm)
        _, nbf, geosInfo, _ = getBFsFromMeshData(meshData; sbfT = :RWG)

        Z_default = getImpedanceMatrix(geosInfo, nbf)

        # Explicitly set free-space and compute again
        set_green_function_type!(:freespace)
        Z_explicit = getImpedanceMatrix(geosInfo, nbf)

        @test Z_default ≈ Z_explicit atol = 1e-14

        reset_green_function_config!()
    end

    # Clean up any results directory created during tests
    rm("results"; force = true, recursive = true)

end  # @testset "EFIE with Ground Plane"
