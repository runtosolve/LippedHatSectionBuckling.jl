include("../src/LippedHatSectionBuckling.jl")
using cFSM
using CUFSM

# ---- section definition ---------------------------------------------------
E = 29500.0
ν = 0.30

t = 0.0188
D = 0.276
L = 0.551
B = 2.01
H = 1.54
h = H + t*2          # centerline web depth (accounts for thickness at flanges)
r = 0.0790

material   = LippedHatSectionBuckling.Material(E, ν)
dimensions = LippedHatSectionBuckling.Dimensions(t, h, B, L, D, r)

# ---- local buckling half-wavelength estimate (dominant plate dimension) ---
Lcrl_estimate = B > H ? B : H
println("Local Lcrl estimate = ", round(Lcrl_estimate; digits=3), " in")

# ---- logarithmically-spaced lengths around the expected Lcrl --------------
L_min   = min(B, H) / 2
L_max   = 2.0 * Lcrl_estimate
lengths = exp.(range(log(L_min), log(L_max); length=50))


# ==========================================================================
# constrained FSM (cFSM), local-only modes
# ==========================================================================
cfsm_section = LippedHatSectionBuckling.get_cFSM_section(
    dimensions, material, LippedHatSectionBuckling.Load(0.0, -1.0, 0.0, 0.0, 0.0);
    n_per_segment=4
)

result_cfsm = cFSM.cfsm_analysis(
    cfsm_section.node, cfsm_section.elem, cfsm_section.prop,
    lengths, 1, ["L"]; BC="S-S"
)
Rcr_per_length = [minimum(result_cfsm.curve[l][:, 2]) for l in eachindex(result_cfsm.curve)]
Lcrl_cfsm = result_cfsm.lengths[argmin(Rcr_per_length)]
Rcrl_cfsm = minimum(Rcr_per_length)

# ==========================================================================
# CUFSM (unconstrained FSM), full signature curve
# ==========================================================================
coordinates = LippedHatSectionBuckling.get_section_coordinates(dimensions)
model_cufsm = CUFSM.Tools.open_section_analysis(
    coordinates.X, coordinates.Y, t, lengths, E, ν,
    0.0, -1.0, 0.0, 0.0, 0.0, [], [], 1
)
Rcrl_curve_cufsm = CUFSM.Tools.get_load_factor(model_cufsm, 1)
Lcrl_cufsm = lengths[argmin(Rcrl_curve_cufsm)]
Rcrl_cufsm = minimum(Rcrl_curve_cufsm)

# ==========================================================================
# Comparison
# ==========================================================================
println("\n--- Mcrl_xx comparison (Mxx = 1 reference load) ---")
println("          CUFSM      cFSM")
println("Lcrl  = ", rpad(round(Lcrl_cufsm; digits=3), 10), round(Lcrl_cfsm; digits=3), " in")
println("Rcrl  = ", rpad(round(Rcrl_cufsm; digits=4), 10), round(Rcrl_cfsm; digits=4), "  (load factor)")
println("Diff  = ", round(100*(Rcrl_cfsm - Rcrl_cufsm)/Rcrl_cufsm; digits=2), " %  (cFSM relative to CUFSM)")
