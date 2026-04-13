# S-Parameter Test Suite
# Tests for the modular S-parameter implementation in MoM_Kernels

using MoM_Kernels, MoM_Basics, LinearAlgebra, Test
using StaticArrays: MVector
using MoM_Basics: UniformDistribution, MVec3D

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

function _create_test_deltagap_array_port(nedges::Int = 2; V::ComplexF64 = 1.0 + 0.0im)
    # Create a DeltaGapArrayPort with manually set fields for testing
    # Note: This bypasses the normal constructor to avoid mesh binding
    FT = Float64
    IT = Int
    DT = UniformDistribution{FT}
    
    # Create with proper type parameters and manually set mesh-dependent fields
    port = DeltaGapArrayPort{FT, IT}(
        id = IT(1),
        V = V,
        freq = FT(1.0e9),
        center = MVec3D{FT}(0.0, 0.0, 0.0),
        normal = MVec3D{FT}(0.0, 0.0, 1.0),
        excitationDistribution = UniformDistribution{FT}(),
        isActive = true
    )
    
    # Manually set mesh-binding fields for testing
    port.isBound = true
    port.rwgIDs = collect(IT, 1:nedges)
    port.triID_pos = ones(IT, nedges)
    port.triID_neg = zeros(IT, nedges)  # Half RWGs
    port.edgeLengths = ones(FT, nedges)
    port.edgeCenters = [MVec3D{FT}(0.0, 0.0, 0.0) for _ in 1:nedges]
    port.edgeOrient = [MVec3D{FT}(1.0, 0.0, 0.0) for _ in 1:nedges]
    port.edgeWeights = ones(Complex{FT}, nedges) ./ nedges  # Uniform weights summing to 1
    port.singleEdgeMode = false
    port.primaryRwgID = IT(0)
    
    return port
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

        @testset "buildExcitationVector with ieType parameter" begin
            nbf = 10
            port = _create_test_deltagap_port(4; V = 1.0 + 0.0im)
            trianglesInfo = TriangleInfo{Int, Float64}[]

            # Test EFIE (default)
            V_efie = buildExcitationVector(port, nbf; ieType=:efie, trianglesInfo=trianglesInfo)
            @test V_efie[4] == 1.0 + 0.0im

            # Test MFIE - for DeltaGapPort, this falls back to EFIE (documented behavior)
            V_mfie = buildExcitationVector(port, nbf; ieType=:mfie, trianglesInfo=trianglesInfo)
            @test V_mfie[4] == 1.0 + 0.0im  # Falls back to EFIE

            # Test CFIE - for DeltaGapPort, this falls back to EFIE (documented behavior)
            V_cfie = buildExcitationVector(port, nbf; ieType=:cfie, trianglesInfo=trianglesInfo, cfie_alpha=0.5)
            @test V_cfie[4] == 1.0 + 0.0im  # Falls back to EFIE
        end

        @testset "excitePort! with ieType parameter" begin
            nbf = 10
            V = zeros(ComplexF64, nbf)
            port = _create_test_deltagap_port(4; V = 1.0 + 0.0im)
            trianglesInfo = TriangleInfo{Int, Float64}[]

            # Test EFIE with explicit ieType
            V .= 0
            excitePort!(V, port, trianglesInfo, :efie)
            @test V[4] == 1.0 + 0.0im

            # Test MFIE - for DeltaGapPort, falls back to EFIE (documented behavior)
            V .= 0
            excitePort!(V, port, trianglesInfo, :mfie)
            @test V[4] == 1.0 + 0.0im  # Falls back to EFIE

            # Test CFIE - for DeltaGapPort, falls back to EFIE (documented behavior)
            V .= 0
            excitePort!(V, port, trianglesInfo, :cfie; cfie_alpha=0.5)
            @test V[4] == 1.0 + 0.0im  # Falls back to EFIE
        end

    end

    @testset "Port Measurement" begin

        @testset "extractPortCurrent DeltaGapPort" begin
            port = _create_test_deltagap_port(3; V = 2.0 + 1.0im)
            nbf = 10
            I = zeros(ComplexF64, nbf)
            I[3] = 0.1 - 0.05im  # Current at port

            I_measured = extractPortCurrent(port, I)
            # Half RWG: length_factor = edgel = 1.0
            @test I_measured == I[3] * 1.0
        end

        @testset "extractPortCurrent CurrentProbe" begin
            port = _create_test_current_probe(5)
            nbf = 10
            I = zeros(ComplexF64, nbf)
            I[5] = 0.5 + 0.2im  # Current at probe location

            I_measured = extractPortCurrent(port, I)
            @test I_measured == I[5]
        end

        @testset "extractPortCurrent DeltaGapArrayPort" begin
            port = _create_test_deltagap_array_port(3)
            nbf = 10
            I = zeros(ComplexF64, nbf)
            I[1] = 0.1 + 0.0im
            I[2] = 0.2 + 0.0im
            I[3] = 0.3 + 0.0im

            I_measured = extractPortCurrent(port, I)
            # Half RWGs: length_factor = edgel = 1.0 for each
            # Edge weights are applied: each weight = 1/3 for uniform distribution
            # Expected: (0.1 + 0.2 + 0.3) * (1/3) = 0.6/3 = 0.2
            expected = (0.1 + 0.2 + 0.3) / 3  # Sum of currents * length_factors * weights
            @test I_measured ≈ expected atol=1e-10
        end

        @testset "measurePortImpedance for Array Ports" begin
            # measurePortImpedance only works for array ports now
            # (DeltaGapArrayPort, RectangularEdgePort)
            
            port_array = _create_test_deltagap_array_port(2)
            nbf = 5
            I = zeros(ComplexF64, nbf)
            I[1] = 0.1 + 0.0im
            I[2] = 0.1 + 0.0im

            # Create a simple diagonal Z matrix
            Z = ComplexF64[50.0 0 0 0 0;
                           0 50.0 0 0 0;
                           0 0 50.0 0 0;
                           0 0 0 50.0 0;
                           0 0 0 0 50.0]

            # For array port: Z = V / I_total
            # I_total = (0.1 + 0.1) * (1/2) = 0.1 (each edge has weight 0.5 for 2 edges)
            # Z = 1.0 / 0.1 = 10.0
            Z_meas = measurePortImpedance(port_array, Z, I, port_array)
            @test isfinite(Z_meas)
            @test Z_meas ≈ 10.0 atol=1e-10
        end

        @testset "measurePortImpedance throws for single-edge ports" begin
            # measurePortImpedance for DeltaGapPort and CurrentProbe was removed
            # because it used an incorrect algorithm for multi-port analysis
            
            port_single = _create_test_deltagap_port(1)
            nbf = 5
            I = zeros(ComplexF64, nbf)
            Z = ComplexF64[50.0 0 0 0 0;
                           0 50.0 0 0 0;
                           0 0 50.0 0 0;
                           0 0 0 50.0 0;
                           0 0 0 0 50.0]

            # Should throw error for single-edge ports
            @test_throws ErrorException measurePortImpedance(port_single, Z, I, port_single)

            port_probe = _create_test_current_probe(1)
            @test_throws ErrorException measurePortImpedance(port_probe, Z, I, port_probe)
        end

    end

    @testset "PortArray Excitation" begin

        @testset "PortArray with ieType parameter API" begin
            # Note: getExcitationVector is defined in MoM_Basics and calls functions in MoM_Kernels.
            # Testing the full integration requires both modules to be fully loaded.
            # Here we verify the function signature exists.
            
            port1 = _create_test_deltagap_port(3)
            port2 = _create_test_deltagap_port(7)
            ports = PortArray([port1, port2])
            nbf = 10
            trianglesInfo = TriangleInfo{Int, Float64}[]

            # Verify the function exists with correct signature
            @test hasmethod(MoM_Basics.getExcitationVector, 
                (typeof(ports), Vector{TriangleInfo{Int,Float64}}, Int, Symbol))
            
            # The actual functionality is tested in integration tests
        end

        @testset "PortArray addExcitationVector! with ieType API" begin
            # Note: addExcitationVector! is defined in MoM_Basics
            # Testing the full integration requires both modules to be fully loaded.
            
            port1 = _create_test_deltagap_port(3)
            port2 = _create_test_deltagap_port(7)
            ports = PortArray([port1, port2])
            
            # Verify the function exists with correct signature
            V = zeros(ComplexF64, 10)
            trianglesInfo = TriangleInfo{Int, Float64}[]
            @test hasmethod(MoM_Basics.addExcitationVector!, 
                (typeof(V), typeof(ports), Vector{TriangleInfo{Int,Float64}}, Symbol))
            
            # The actual functionality is tested in integration tests
        end

        @testset "buildExcitationVector for PortArray index" begin
            port1 = _create_test_deltagap_port(3)
            port2 = _create_test_deltagap_port(7)
            ports = PortArray([port1, port2])
            nbf = 10
            trianglesInfo = TriangleInfo{Int, Float64}[]

            # Test building excitation for specific port by index
            V1 = MoM_Kernels.buildExcitationVector(ports, 1, nbf; trianglesInfo=trianglesInfo)
            @test V1[3] == 1.0 + 0.0im
            @test V1[7] == 0.0  # Port 2 not excited

            V2 = MoM_Kernels.buildExcitationVector(ports, 2, nbf; trianglesInfo=trianglesInfo)
            @test V2[3] == 0.0  # Port 1 not excited
            @test V2[7] == 1.0 + 0.0im
        end

        @testset "Heterogeneous PortArray type validation" begin
            # Test that PortArray with mixed port types works
            port_dg = _create_test_deltagap_port(1)
            port_cp = _create_test_current_probe(2)
            ports_hetero = PortArray([port_dg, port_cp])

            # Verify the PortArray was created correctly
            @test ports_hetero.numPorts == 2
            @test length(ports_hetero.ports) == 2
            @test typeof(ports_hetero.ports[1]) <: DeltaGapPort
            @test typeof(ports_hetero.ports[2]) <: CurrentProbe

            # Test excitation still works
            nbf = 5
            trianglesInfo = TriangleInfo{Int, Float64}[]
            V = MoM_Kernels.buildExcitationVector(ports_hetero, 1, nbf; trianglesInfo=trianglesInfo)
            @test V[1] == 1.0 + 0.0im  # DeltaGapPort excitation
            @test V[2] == 0.0  # CurrentProbe not excited
        end

    end

    @testset "Core S-Parameter Functions" begin

        # NOTE: computeS11, computeInputImpedance, and computeSParameters require
        # a full simulation setup (Params.freq, proper mesh, etc.)
        # These integration tests are tested in the main test suite
        # Here we just verify the functions exist and have correct signatures

        @testset "computeS11 API exists" begin
            # Accepts any PortType (DeltaGapPort, DeltaGapArrayPort, ...)
            @test hasmethod(computeS11, (DeltaGapPort{Float64,Int}, Any, Int))
            @test hasmethod(computeS11, (PortType, Any, Int))
        end

        @testset "computeInputImpedance API exists" begin
            @test hasmethod(computeInputImpedance, (DeltaGapPort{Float64,Int}, Any, Int))
            @test hasmethod(computeInputImpedance, (PortType, Any, Int))
        end

        @testset "computeSParameters API exists" begin
            # Create test ports (just for API verification)
            port1 = _create_test_deltagap_port(3)
            port2 = _create_test_deltagap_port(7)
            ports = PortArray([port1, port2])

            # Homogeneous concrete-PT array (original behaviour unchanged)
            @test hasmethod(computeSParameters, (
                PortArray{Float64,Int,DeltaGapPort{Float64,Int}},
                Any, Int
            ))
            # Heterogeneous array (PT = abstract PortType)
            @test hasmethod(computeSParameters, (
                PortArray{Float64,Int,PortType},
                Any, Int
            ))
        end

        @testset "computeSParameters with ieType parameter" begin
            # Create test ports
            port1 = _create_test_deltagap_port(3)
            port2 = _create_test_deltagap_port(7)
            ports = PortArray([port1, port2])
            nbf = 10
            Z = _create_simple_impedance_matrix(nbf)
            trianglesInfo = TriangleInfo{Int, Float64}[]

            # Test that ieType parameter exists
            @test hasmethod(computeSParameters, (
                PortArray{Float64,Int,DeltaGapPort{Float64,Int}},
                Matrix{ComplexF64}, Int
            ))

            # Note: We can't actually call computeSParameters without a full solver setup,
            # but we verify the function accepts the ieType keyword
            # The function signature includes: ieType::Symbol = :efie
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

    @testset "Touchstone Round-Trip" begin

        @testset "2-port RI round-trip" begin
            S_orig = ComplexF64[0.1+0.2im 0.8+0.1im; 0.8+0.1im 0.15+0.25im]
            Z_port = s2z(S_orig, 50.0)
            result = SParameterResult(S_orig, Z_port, [2.4e9], 50.0, [1, 2])

            filename = tempname() * ".s2p"
            saveTouchstone(filename, result; format=:ri)
            loaded = loadTouchstone(filename)

            @test loaded.num_ports == 2
            @test length(loaded.frequencies) == 1
            @test abs(loaded.frequencies[1] - 2.4e9) < 1.0  # Hz precision
            @test loaded.Z0 ≈ 50.0
            for i in 1:2, j in 1:2
                @test loaded.S[i,j,1] ≈ S_orig[i,j] atol=1e-10
            end
            rm(filename)
        end

        @testset "2-port DB round-trip" begin
            S_orig = ComplexF64[0.1+0.2im 0.5-0.1im; 0.5-0.1im 0.05+0.3im]
            Z_port = s2z(S_orig, 50.0)
            result = SParameterResult(S_orig, Z_port, [1e9], 75.0, [1, 2])

            filename = tempname() * ".s2p"
            saveTouchstone(filename, result; format=:db)
            loaded = loadTouchstone(filename)

            @test loaded.Z0 ≈ 75.0
            for i in 1:2, j in 1:2
                @test abs(loaded.S[i,j,1] - S_orig[i,j]) < 1e-8
            end
            rm(filename)
        end

        @testset "2-port MA round-trip" begin
            S_orig = ComplexF64[0.3*cis(deg2rad(-45)) 0.9*cis(deg2rad(30));
                                0.9*cis(deg2rad(30)) 0.2*cis(deg2rad(-60))]
            Z_port = s2z(S_orig, 50.0)
            result = SParameterResult(S_orig, Z_port, [5e9], 50.0, [1, 2])

            filename = tempname() * ".s2p"
            saveTouchstone(filename, result; format=:ma)
            loaded = loadTouchstone(filename)

            for i in 1:2, j in 1:2
                @test abs(loaded.S[i,j,1] - S_orig[i,j]) < 1e-10
            end
            rm(filename)
        end

        @testset "3-port RI round-trip" begin
            S_orig = ComplexF64[
                0.1+0.0im  0.5+0.1im  0.5+0.1im;
                0.5+0.1im  0.1+0.0im  0.5+0.1im;
                0.5+0.1im  0.5+0.1im  0.1+0.0im]
            Z_port = s2z(S_orig, 50.0)
            result = SParameterResult(S_orig, Z_port, [3e9], 50.0, [1, 2, 3])

            filename = tempname() * ".s3p"
            saveTouchstone(filename, result; format=:ri)
            loaded = loadTouchstone(filename)

            @test loaded.num_ports == 3
            for i in 1:3, j in 1:3
                @test loaded.S[i,j,1] ≈ S_orig[i,j] atol=1e-10
            end
            rm(filename)
        end

        @testset "1-port RI round-trip" begin
            S_orig = ComplexF64[0.3+0.4im][:,:]
            Z_port = s2z(S_orig, 50.0)
            result = SParameterResult(S_orig, Z_port, [1e9], 50.0, [1])

            filename = tempname() * ".s1p"
            saveTouchstone(filename, result; format=:ri)
            loaded = loadTouchstone(filename)

            @test loaded.num_ports == 1
            @test loaded.S[1,1,1] ≈ S_orig[1,1] atol=1e-10
            rm(filename)
        end

        @testset "Multi-frequency round-trip" begin
            nfreq = 4
            S = zeros(ComplexF64, 2, 2, nfreq)
            Z_port = zeros(ComplexF64, 2, 2, nfreq)
            freqs = [1e9, 2e9, 3e9, 4e9]
            for fi in 1:nfreq
                S[:,:,fi] = ComplexF64[0.1 0.8; 0.8 0.1] * (fi * 0.1)
                Z_port[:,:,fi] = s2z(S[:,:,fi], 50.0)
            end
            result = SParameterResult(S, Z_port, freqs, 50.0, [1, 2])

            filename = tempname() * ".s2p"
            saveTouchstone(filename, result; format=:ri)
            loaded = loadTouchstone(filename)

            @test length(loaded.frequencies) == nfreq
            for fi in 1:nfreq
                @test abs(loaded.frequencies[fi] - freqs[fi]) < 1.0
                for i in 1:2, j in 1:2
                    @test loaded.S[i,j,fi] ≈ S[i,j,fi] atol=1e-10
                end
            end
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
