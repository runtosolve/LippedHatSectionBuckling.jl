include("../src/LippedHatSectionBuckling.jl")
using Test

# ── shared inputs ────────────────────────────────────────────────────────────
E = 29500.0
ν = 0.30 

t = 0.0188
D = 0.276
L = 0.551

B = 2.01
H = 1.54
r = t


material = LippedHatSectionBuckling.Material(E, ν)
dimensions = LippedHatSectionBuckling.Dimensions(t, H, B, L, D, r)

# ── cross-section geometry ───────────────────────────────────────────────────
@testset "get_section_coordinates" begin
    coords = LippedHatSectionBuckling.get_section_coordinates(dimensions)

    @test length(coords.X) == length(coords.Y)
    @test length(coords.X) > 0
    @test all(isfinite, coords.X)
    @test all(isfinite, coords.Y)
    # origin shifted so minimum is at t/2
    @test minimum(coords.X) ≈ t / 2  atol=1e-6
    @test minimum(coords.Y) ≈ t / 2  atol=1e-6
end

# ── distortional wavelength ──────────────────────────────────────────────────
@testset "calculate_Lcrd" begin
    Lcrd_P = LippedHatSectionBuckling.calculate_Lcrd(dimensions, material, "P")
    Lcrd_M = LippedHatSectionBuckling.calculate_Lcrd(dimensions, material, "M")

    @test isfinite(Lcrd_P)
    @test isfinite(Lcrd_M)
    @test Lcrd_P > 0
    @test Lcrd_M > 0
end

# ── axial buckling ───────────────────────────────────────────────────────────
@testset "calculate_Pcrℓ" begin
    sec = LippedHatSectionBuckling.calculate_Pcrℓ(dimensions, material)

    @test sec.label == "Pcrℓ"
    @test sec.results.Rcr > 0
    @test sec.results.Lcr > 0
    @test isfinite(sec.results.Rcr)
end

@testset "calculate_Pcrd" begin
    sec = LippedHatSectionBuckling.calculate_Pcrd(dimensions, material)

    @test sec.label == "Pcrd"
    @test sec.results.Rcr > 0
    @test sec.results.Lcr > 0
    @test isfinite(sec.results.Rcr)
    # distortional load factor should be >= local
    sec_l = LippedHatSectionBuckling.calculate_Pcrℓ(dimensions, material)
    @test sec.results.Rcr >= sec_l.results.Rcr
end

# ── strong-axis bending buckling ─────────────────────────────────────────────
@testset "calculate_Mcrd_xx" begin
    sec = LippedHatSectionBuckling.calculate_Mcrd_xx(dimensions, material)

    @test sec.label == "Mcrd_xx"
    @test sec.results.Rcr > 0
    @test isfinite(sec.results.Rcr)
end

@testset "calculate_Mcrℓ_xx" begin
    sec = LippedHatSectionBuckling.calculate_Mcrℓ_xx(dimensions, material)

    @test sec.label == "Mcrℓ_xx"
    @test sec.results.Rcr > 0
    @test isfinite(sec.results.Rcr)
end

# ── weak-axis bending buckling ───────────────────────────────────────────────
@testset "calculate_Mcrℓ_yy_neg" begin
    sec = LippedHatSectionBuckling.calculate_Mcrℓ_yy_neg(dimensions, material)

    @test sec.label == "Mcrℓ_yy_neg"
    @test sec.results.Rcr > 0
    @test isfinite(sec.results.Rcr)
end
