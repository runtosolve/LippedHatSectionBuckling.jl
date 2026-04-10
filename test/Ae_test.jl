include("../src/LippedHatSectionBuckling.jl")

# ---- section definition ---------------------------------------------------
E = 29500.0
ν = 0.30

t = 0.0219
D = 0.276
L = 0.551
B = 2.01
H = 1.54
h = H + t*2
r = 0.0790

material   = LippedHatSectionBuckling.Material(E, ν)
dimensions = LippedHatSectionBuckling.Dimensions(t, h, B, L, D, r)

# ---- compressive stress ---------------------------------------------------
# f = Fn from Chapter E; use Fy here as an upper bound (fully yielded case)
Fy = 70.0   # ksi — change to actual Fn as needed
f  = Fy

# ---- flat widths ----------------------------------------------------------
w = LippedHatSectionBuckling.flat_widths(dimensions)
println("Flat widths (in):")
println("  D1 (return lip, left)  = ", round(w.D1; digits=4))
println("  L1 (lip, left)         = ", round(w.L1; digits=4))
println("  B1 (flange, left)      = ", round(w.B1; digits=4))
println("  H  (web)               = ", round(w.H;  digits=4))
println("  B2 (flange, right)     = ", round(w.B2; digits=4))
println("  L2 (lip, right)        = ", round(w.L2; digits=4))
println("  D2 (return lip, right) = ", round(w.D2; digits=4))
println()

# ---- effective area -------------------------------------------------------
ae = LippedHatSectionBuckling.calculate_Ae(dimensions, material, f, Fy)

println("Effective widths at f = ", f, " ksi:")
println("  Element                  beff (in)   Fcrl (ksi)   lambda")
for (el, b, fcrl, lam) in (
    ("D1 (unstiffened)         ", ae.beff.D1, ae.frcl.D1, ae.lambda.D1),
    ("L1 (Partially stiffened) ", ae.beff.L1, ae.frcl.L1, ae.lambda.L1),
    ("B1 (stiffened)           ", ae.beff.B1, ae.frcl.B1, ae.lambda.B1),
    ("H  (stiffened)           ", ae.beff.H,  ae.frcl.H,  ae.lambda.H),
    ("B2 (stiffened)           ", ae.beff.B2, ae.frcl.B2, ae.lambda.B2),
    ("L2 (Partially stiffened) ", ae.beff.L2, ae.frcl.L2, ae.lambda.L2),
    ("D2 (unstiffened)         ", ae.beff.D2, ae.frcl.D2, ae.lambda.D2),
)
    println("  ", el, "   ", rpad(round(b;    digits=4), 8),
                             "   ", rpad(round(fcrl; digits=2), 10),
                             "   ", round(lam; digits=4))
end
println()
# println("Section 1.3 lip stiffener (L1 = L2 by symmetry):")
# println("  S   = ", round(ae.sec13.L1.S;  digits=4))
# println("  Ia  = ", round(ae.sec13.L1.Ia; digits=8), " in⁴")
# println("  Is  = ", round(ae.sec13.L1.Is; digits=8), " in⁴")
# println("  R1  = ", round(ae.sec13.L1.R1; digits=4))
# println("  k   = ", round(ae.sec13.L1.k;  digits=4))
# println()
println("Total fillet area (fully effective) = ", round(ae.A_corners; digits=4), " in²")
println("Effective area, Ae                       = ", round(ae.Ae;        digits=4), " in²")

