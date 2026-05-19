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

# ---- AISIS100 estimate for Lcrd (sets the length scan range) --------------
Lcrd_estimate = LippedHatSectionBuckling.calculate_Lcrd(dimensions, material, "M")
println("AISIS100 Lcrd estimate = ", round(Lcrd_estimate; digits=3), " in")

# ---- logarithmically-spaced lengths around the expected Lcrd --------------
L_min   = min(B, H) / 4
L_max   = 4.0 * Lcrd_estimate
lengths = exp.(range(log(L_min), log(L_max); length=60))

# ==========================================================================
# constrained FSM (cFSM), distortional-only modes
# ==========================================================================
cfsm_section = LippedHatSectionBuckling.get_cFSM_section(
    dimensions, material, LippedHatSectionBuckling.Load(0.0, -1.0, 0.0, 0.0, 0.0);
    n_per_segment=4
)

result_cfsm = cFSM.cfsm_analysis(
    cfsm_section.node, cfsm_section.elem, cfsm_section.prop,
    lengths, 1, ["D"]; BC="S-S"
)
Rcr_per_length = [minimum(result_cfsm.curve[l][:, 2]) for l in eachindex(result_cfsm.curve)]
Lcrd_cfsm = result_cfsm.lengths[argmin(Rcr_per_length)]
Rcrd_cfsm = minimum(Rcr_per_length)