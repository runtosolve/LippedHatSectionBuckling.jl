include("../src/LippedHatSectionBuckling.jl")
using cFSM

# ---- section definition ---------------------------------------------------
E = 29500.0
ν = 0.30

t = 0.0188
D = 0.276
L = 0.551
B = 2.01
H = 1.54
h = H + t*2
r = 0.0790

material   = LippedHatSectionBuckling.Material(E, ν)
dimensions = LippedHatSectionBuckling.Dimensions(t, h, B, L, D, r)

# ---- build straight-sided cFSM model (no arc nodes) ----------------------
# P=1 axial reference load; n_per_segment=4 gives 3 collinear sub-nodes/segment
cfsm_section = LippedHatSectionBuckling.get_cFSM_section(
    dimensions, material, LippedHatSectionBuckling.Load(1.0, 0.0, 0.0, 0.0, 0.0);
    n_per_segment=4
)

# ---- AISIS100 estimate for Lcrd (used only to set the length scan range) --
Lcrd_estimate = LippedHatSectionBuckling.calculate_Lcrd(dimensions, material, "P")
println("AISIS100 Lcrd estimate = ", round(Lcrd_estimate; digits=3), " in")
hi

# ---- logarithmically-spaced lengths around the expected Lcrd --------------
L_min = min(B, H) / 2
L_max = 2.0 * Lcrd_estimate
lengths = exp.(range(log(L_min), log(L_max); length=60))

# ---- constrained FSM distortional-only signature curve --------------------
result = cFSM.cfsm_analysis(cfsm_section.node, cfsm_section.elem, cfsm_section.prop, lengths, modes = ["D"], BC = "S-S")

println("cFSM  Lcrd = ", round(result.Lcr; digits=3), " in")
println("cFSM  Rcrd = ", round(result.Rcr; digits=4), "  (load factor, P=1 reference)")


