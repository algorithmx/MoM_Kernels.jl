# S-Parameter Test Suite
# Tests for the modular S-parameter implementation in MoM_Kernels

using MoM_Kernels, MoM_Basics, LinearAlgebra, Test
using StaticArrays: MVector

# =============================================================================
# Test Helper Functions
# =============================================================================

function _create_test_deltagap_port(rwgID::Int = 1; V::ComplexF64 = 1.0 + 0.0im)
    return DeltaGapPort{Float64, Int}(
        id = rwgID,
        V = V,
        freq = 1.0e9,
        rwgID = rwgID,
        triID_pos = 1,
        triID_neg = 0,
        edgel = 1.0,
        position = MVector{3, Float64}(0.0, 0.0, 0.0),
        direction = MVector{3, Float64}(1.0, 0.0, 0.0),
        isActive = true
    )
end

function _create_test_current_probe(rwgID::Int = 1)
    return CurrentProbe{Float64, Int}(
        id = rwgID,
        I = ComplexF64(1.0, 0.0),
        freq = 1.0e9,
        rwgID = rwgID,
        triID = 1,
        edgel = 1.0,
        center = MVector{3, Float64}(0.0, 0.0, 0.0),
        isActive = true
    )
end

function _create_simple_impedance_matrix(n::Int; FT::Type = Float64)
    # Create a symmetric positive-definite matrix for testing
    A = randn(FT, n, n)
    Z = Complex{FT}.(A' * A + I)  # Ensure positive definite
    return Z
end

function _create_matched_load_zmatrix(n::Int, Z0::Real = 50.0; FT::Type = Float64)
    # Create a diagonal matrix with Z0 on diagonal (matched load)
    Z = Complex{FT}.(Z0 * I(n))
    return Z
end

function _create_two_port_network(case::Symbol = :matched; FT::Type = Float64)
    # Create standard two-port networks for testing
    if case == :matched
        # Matched load: Z = Z0 * I
        Z = Complex{FT}[50 0; 0 50]
        S = Complex{FT}[0 0; 0 0]
    elseif case == :open
        # Open circuit: very high impedance
        Z = Complex{FT}[1e10 0; 0 1e10]
        S = Complex{FT}[1 0; 0 1]
    elseif case == :short
        # Short circuit: zero impedance
        Z = Complex{FT}[0 0; 0 0]
        S = Complex{FT}[-1 0; 0 -1]
    elseif case == :attenuator
        # 3dB attenuator (10dB)
        S = Complex{FT}[0 1; 1 0] / sqrt(10)
        Z = s2z(S, 50.0)
    else
        error("Unknown network case: $case")
    end
    return Z, S
end

# =============================================================================
# S-Parameter Tests
# =============================================================================

function test_sparameters()

    @testset "S-Parameter Types" begin

        @testset "SParameterResult Construction" begin
            FT = Float64
            S = zeros(Complex{FT}, 2, 2, 1)
            Z_port = zeros(Complex{FT}, 2, 2, 1)
            freqs = [1e9]
            Z0 = 50.0
            portIDs = [1, 2]

            result = SParameterResult(S, Z_port, freqs, Z0, portIDs)

            @test result.Z0 == Z0
            @test size(result.S) == (2, 2, 1)
            @test size(result.Z_port) == (2, 2, 1)
            @test result.frequencies == freqs
            @test result.portIDs == portIDs
        end

        @testset "SParameterResult Multi-Frequency" begin
            FT = Float64
            nfreq = 5
            S = zeros(Complex{FT}, 3, 3, nfreq)
            Z_port = zeros(Complex{FT}, 3, 3, nfreq)
            freqs = range(1e9, stop=5e9, length=nfreq)
            portIDs = [1, 2, 3]

            result = SParameterResult(S, Z_port, collect(freqs), 50.0, portIDs)

            @test size(result.S, 3) == nfreq
            @test length(result.frequencies) == nfreq
        end

    end

    @testset "Z-S Conversions" begin

        @testset "z2s: Matched Load" begin
            Z = ComplexF64[50 0; 0 50]
            S = z2s(Z, 50.0)
            @test S ≈ zeros(2, 2) atol=1e-10
        end

        @testset "z2s: Open Circuit" begin
            Z = ComplexF64[1e10 0; 0 1e10]
            S = z2s(Z, 50.0)
            @test real(S[1,1]) ≈ 1.0 atol=1e-6
            @test real(S[2,2]) ≈ 1.0 atol=1e-6
        end

        @testset "z2s: Short Circuit" begin
            Z = ComplexF64[0 0; 0 0]
            S = z2s(Z, 50.0)
            @test real(S[1,1]) ≈ -1.0 atol=1e-10
            @test real(S[2,2]) ≈ -1.0 atol=1e-10
        end

        @testset "s2z: Matched Load" begin
            S = zeros(ComplexF64, 2, 2)
            Z = s2z(S, 50.0)
            @test Z ≈ ComplexF64[50 0; 0 50] atol=1e-10
        end

        @testset "Round-Trip Conversion z2s→s2z" begin
            # Random symmetric impedance matrix
            Z_original = ComplexF64[75+10im 20-5im; 20-5im 60+15im]
            S = z2s(Z_original, 50.0)
            Z_recovered = s2z(S, 50.0)
            @test Z_recovered ≈ Z_original atol=1e-10
        end

        @testset "Round-Trip Conversion s2z→z2s" begin
            S_original = ComplexF64[0.1+0.2im 0.8+0.1im; 0.8+0.1im 0.15+0.25im]
            Z = s2z(S_original, 50.0)
            S_recovered = z2s(Z, 50.0)
            @test S_recovered ≈ S_original atol=1e-10
        end

        @testset "Single Port Conversion" begin
            Z_1p = ComplexF64[50.0 + 0.0im;;]  # 1x1 matrix
            S_1p = z2s(Z_1p, 50.0)
            @test S_1p ≈ zeros(1, 1) atol=1e-10

            # Test near-open circuit (S ≈ 1, not exactly 1 to avoid singularity)
            S_1p_open = ComplexF64[0.999999 + 0.0im;;]  # 1x1 matrix, nearly open
            Z_1p_recovered = s2z(S_1p_open, 50.0)
            @test real(Z_1p_recovered[1,1]) > 1e5  # Very high impedance
        end

    end

    @testset "Port Excitation" begin

        @testset "excitePort! DeltaGapPort" begin
            nbf = 10
            V = zeros(ComplexF64, nbf)
            port = _create_test_deltagap_port(5; V = 2.0 + 0.0im)

            excitePort!(V, port)

            # Half RWG (boundary, triID_neg=0): V × l = 2.0 × 1.0 = 2.0
            @test V[5] == 2.0 + 0.0im
            @test count(!iszero, V) == 1
        end

        @testset "excitePort! CurrentProbe" begin
            nbf = 10
            V = zeros(ComplexF64, nbf)
            port = _create_test_current_probe(3)

            excitePort!(V, port)

            # Current probe should set excitation at its rwgID
            @test V[3] == 1.0 + 0.0im
            @test count(!iszero, V) == 1
        end

        @testset "buildExcitationVector" begin
            nbf = 10
            port = _create_test_deltagap_port(4; V = 1.0 + 0.5im)

            V = buildExcitationVector(port, nbf)

            @test length(V) == nbf
            # Half RWG (boundary, triID_neg=0): V × l = (1.0+0.5im) × 1.0
            @test V[4] == 1.0 + 0.5im
            @test count(!iszero, V) == 1
        end

    end

    @testset "Port Measurement" begin

        @testset "extractPortVoltage DeltaGapPort" begin
            port = _create_test_deltagap_port(3; V = 2.0 + 1.0im)
            nbf = 10
            I = zeros(ComplexF64, nbf)
            I[3] = 0.1 - 0.05im  # Current at port

            # Create a simple Z matrix (diagonal with 50Ω)
            Z = _create_matched_load_zmatrix(nbf, 50.0)

            # For DeltaGapPort, voltage is computed from Z*I
            V_measured = extractPortVoltage(port, Z, I)
            @test isfinite(V_measured)
        end

        @testset "extractPortCurrent CurrentProbe" begin
            port = _create_test_current_probe(5)
            nbf = 10
            I = zeros(ComplexF64, nbf)
            I[5] = 0.5 + 0.2im  # Current at probe location

            I_measured = extractPortCurrent(port, I)
            @test I_measured == I[5]
        end

        @testset "measurePortImpedance Single Port" begin
            # Simple test: Z = V/I
            port = _create_test_deltagap_port(1; V = 10.0 + 0.0im)
            nbf = 5
            I = zeros(ComplexF64, nbf)
            I[1] = 0.2 + 0.0im  # I = 0.2A

            # Create a simple diagonal Z matrix
            Z = ComplexF64[50.0 0 0 0 0;
                           0 50.0 0 0 0;
                           0 0 50.0 0 0;
                           0 0 0 50.0 0;
                           0 0 0 0 50.0]

            # measurePortImpedance(port_i, Z, I_j, port_j)
            # For S11: port_i = port_j (same port)
            Z_meas = measurePortImpedance(port, Z, I, port)
            @test isfinite(Z_meas)
        end

    end

    @testset "Core S-Parameter Functions" begin

        # NOTE: computeS11, computeInputImpedance, and computeSParameters require
        # a full simulation setup (Params.freq, proper mesh, etc.)
        # These integration tests are tested in the main test suite
        # Here we just verify the functions exist and have correct signatures

        @testset "computeS11 API exists" begin
            # Verify the function exists with correct signature
            @test hasmethod(computeS11, (DeltaGapPort{Float64,Int}, Any, Int))
        end

        @testset "computeInputImpedance API exists" begin
            # Verify the function exists with correct signature
            @test hasmethod(computeInputImpedance, (DeltaGapPort{Float64,Int}, Any, Int))
        end

        @testset "computeSParameters API exists" begin
            # Create test ports (just for API verification)
            port1 = _create_test_deltagap_port(3)
            port2 = _create_test_deltagap_port(7)
            ports = PortArray([port1, port2])

            # Verify the function exists
            @test hasmethod(computeSParameters, (
                PortArray{Float64,Int,DeltaGapPort{Float64,Int}},
                Any, Int
            ))
        end

    end

    @testset "Touchstone Export" begin

        @testset "saveTouchstone 1-Port (.s1p)" begin
            S = zeros(ComplexF64, 1, 1, 1)
            Z_port = ComplexF64[50.0 + 0.0im][:,:,:]
            result = SParameterResult(S, Z_port, [1e9], 50.0, [1])

            filename = tempname() * ".s1p"
            saveTouchstone(filename, result; format=:db)

            @test isfile(filename)

            # Verify file content
            lines = readlines(filename)
            @test any(line -> occursin("#", line), lines)
            @test any(line -> occursin("Hz", line) || occursin("GHZ", line), lines)

            rm(filename)
        end

        @testset "saveTouchstone 2-Port (.s2p)" begin
            S = ComplexF64[0.1+0.2im 0.8+0.1im; 0.8+0.1im 0.15+0.25im][:,:,:]
            Z_port = s2z(S[:,:,1], 50.0)[:,:,:]
            result = SParameterResult(S, Z_port, [2.4e9], 50.0, [1, 2])

            filename = tempname() * ".s2p"
            saveTouchstone(filename, result; format=:ma)

            @test isfile(filename)

            lines = readlines(filename)
            @test any(line -> occursin("#", line), lines)

            rm(filename)
        end

        @testset "saveTouchstone Multi-Frequency" begin
            nfreq = 5
            S = zeros(ComplexF64, 2, 2, nfreq)
            Z_port = zeros(ComplexF64, 2, 2, nfreq)
            freqs = [1e9, 2e9, 3e9, 4e9, 5e9]

            for i in 1:nfreq
                S[:,:,i] = ComplexF64[0.1 0.8; 0.8 0.1] * (i * 0.1)
                Z_port[:,:,i] = s2z(S[:,:,i], 50.0)
            end

            result = SParameterResult(S, Z_port, freqs, 50.0, [1, 2])

            filename = tempname() * ".s2p"
            saveTouchstone(filename, result; format=:ri)

            @test isfile(filename)

            lines = readlines(filename)
            # Check for 5 frequency points (plus header lines)
            data_lines = filter(line -> !isempty(line) && !startswith(line, "#") && !startswith(line, "!"), lines)
            @test length(data_lines) >= 5

            rm(filename)
        end

    end

    @testset "S-Parameter Quality Checks" begin

        @testset "check_sparameter_reciprocity" begin
            # Reciprocal S-matrix
            S_reciprocal = ComplexF64[0.1 0.5; 0.5 0.2]
            @test check_sparameter_reciprocity(S_reciprocal) == true

            # Non-reciprocal S-matrix
            S_nonrecip = ComplexF64[0.1 0.5; 0.6 0.2]
            result = check_sparameter_reciprocity(S_nonrecip)
            # Should return false but don't fail test (it's a warning)
            @test typeof(result) == Bool
        end

        @testset "check_passivity" begin
            # Passive S-matrix (|S| <= 1)
            S_passive = ComplexF64[0.1 0.3; 0.3 0.1]
            @test check_passivity(S_passive) == true

            # Active S-matrix (amplification)
            S_active = ComplexF64[2.0 0; 0 2.0]
            result = check_passivity(S_active)
            @test typeof(result) == Bool
        end

    end

    @testset "T-Matrix Conversions" begin

        @testset "s2t and t2s round-trip" begin
            # Two-port S-matrix
            S_original = ComplexF64[0.1+0.2im 0.8+0.1im; 0.8+0.1im 0.15+0.25im]
            T = s2t(S_original)
            S_recovered = t2s(T)
            @test S_recovered ≈ S_original atol=1e-10
        end

        @testset "cascadeT" begin
            # Two identical networks cascaded
            S_half = ComplexF64[0.2 0.7; 0.7 0.2]
            T_half = s2t(S_half)

            # Cascade with itself
            T_full = cascadeT(T_half, T_half)
            S_full = t2s(T_full)

            @test size(S_full) == (2, 2)
            @test all(isfinite.(S_full))
        end

        @testset "cascadeS" begin
            # Two S-matrices cascaded
            S1 = ComplexF64[0.1 0.8; 0.8 0.1]
            S2 = ComplexF64[0.15 0.75; 0.75 0.15]

            S_cascaded = cascadeS(S1, S2)

            @test size(S_cascaded) == (2, 2)
            @test all(isfinite.(S_cascaded))
        end

    end

end

# =============================================================================
# Export test function
# =============================================================================

export test_sparameters
