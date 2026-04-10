"""
    LippedHatSectionBuckling

Elastic buckling analysis for cold-formed steel lipped hat sections using the finite strip
method (CUFSM) and AISIS100 analytical formulas.

The hat section geometry is defined by seven segments: `[D, L, B, H, B, L, D]`
with angles `[-π/2, 0, π/2, 0, -π/2, 0, π/2]`, where:
- `H` — web (top horizontal plate)
- `B` — flange (vertical side)
- `L` — lip (horizontal, extending outward from flange bottom)
- `D` — return lip (vertical edge stiffener at tip)


# Main workflow
```julia
material   = LippedHatSectionBuckling.Material(E, ν)
dimensions = LippedHatSectionBuckling.Dimensions(t, H, B, L, D, r)
section    = LippedHatSectionBuckling.calculate_Pcrℓ(dimensions, material)
Rcr        = section.results.Rcr
```
"""
module LippedHatSectionBuckling


using CrossSectionGeometry, CUFSM, AISIS100, cFSM


struct SectionInput

    E
    ν
    t
    H
    B
    L
    D
    r

end

"""
    Material(E, ν)

Elastic material properties.

# Fields
- `E` — Young's modulus (ksi)
- `ν` — Poisson's ratio
"""
struct Material

    E
    ν

end

"""
    Dimensions(t, H, B, L, D, r)

Cross-section dimensions for a lipped hat section.

# Fields
- `t` — thickness (in)
- `H` — web width (in)
- `B` — flange height (in)
- `L` — lip width (in)
- `D` — return lip width (in)
- `r` — inside bend radius (in)
"""
struct Dimensions

    t
    H  #web    
    B  #flange
    L  #lip
    D  #return lip
    r  #inside radius

end


"""
    Load(P, Mxx, Mzz, M11, M22)

Applied reference loads for a buckling analysis. Set the relevant component to 1.0
and all others to 0.0 to isolate a particular buckling mode.

# Fields
- `P`   — axial load
- `Mxx` — strong-axis bending moment
- `Mzz` — weak-axis bending moment
- `M11` — principal-axis moment 1
- `M22` — principal-axis moment 2
"""
struct Load

    P
    Mxx
    Mzz
    M11
    M22

end

"""
    Results

Output from a finite-strip buckling analysis.

# Fields
- `model` — raw CUFSM model object
- `Lcr`   — critical half-wavelength corresponding to minimum load factor (in)
- `Rcr`   — minimum buckling load factor (dimensionless)
"""
struct Results

    model
    Lcr
    Rcr

end


"""
    Section

Container bundling all inputs and outputs for a single buckling calculation.

# Fields
- `label`      — string identifier for the buckling mode (e.g. `"Pcrℓ"`, `"Mcrd_xx"`)
- `material`   — `Material` object
- `dimensions` — `Dimensions` object
- `load`       — `Load` object used in the analysis
- `results`    — `Results` object with `Lcr` and `Rcr`
"""
struct Section

    label
    material
    dimensions
    load
    results

end




"""
    get_section_coordinates(dimensions) -> NamedTuple{(:X, :Y)}

Build the centerline node coordinates for a lipped hat section from a `Dimensions` object.
The geometry follows the segment sequence `[D, L, B, H, B, L, D]`.
Coordinates are shifted so the lower-left corner starts at `(t/2, t/2)`.

# Returns
A named tuple `(X=..., Y=...)` of centerline node coordinate vectors.
"""
function get_section_coordinates(dimensions)

    (;
    t, 

    H,
    B,
    L,
    D, 
    r  #inside radius 
    ) = dimensions 

    section_dimensions = [D, L, B, H, B, L, D]
    r = [t+r, t+r, r, r, t+r, t+r]
    n = [3, 3, 3, 3, 3, 3, 3]
    n_r = [4, 4, 4, 4, 4, 4];
    θ = [-π/2, 0, π/2, 0, -π/2, 0, π/2]

    centerline = "to left"
    offset = (0.0, 0.0)
    coordinates = CrossSectionGeometry.create_thin_walled_cross_section_geometry(section_dimensions, θ, n, r, n_r, t, centerline=centerline, offset=offset)

    X = [coordinates.centerline_node_XY[i][1] for i in eachindex(coordinates.centerline_node_XY)]
    Y = [coordinates.centerline_node_XY[i][2] for i in eachindex(coordinates.centerline_node_XY)]

    ΔX = -minimum(X) + t[1] / 2
    ΔY = -minimum(Y) + t[1] / 2

    X .+= ΔX
    Y .+= ΔY

    coordinates = (X=X, Y=Y)

    return coordinates

end




"""
    get_cFSM_section(dimensions, material, load; n_per_segment=4) -> NamedTuple

Build a straight-sided (no radius) CUFSM node/elem/prop model suitable for
constrained FSM (cFSM) analysis.

The section geometry is approximated by straight elements only — bend radii are
ignored. This is required because cFSM (GBT-based) mode classification assumes
straight plate elements meeting at sharp corners. Intermediate nodes within each
straight segment are collinear, so `_meta_elems` correctly merges them into
sub-nodes, giving exactly 8 main nodes and 7 meta-elements for the lipped hat.

# Arguments
- `dimensions`      — `Dimensions` object
- `material`        — `Material` object
- `load`            — `Load` object (e.g. `Load(1.0,0.,0.,0.,0.)` for axial)
- `n_per_segment`   — nodes per straight segment including both endpoints (≥2)

# Returns
Named tuple `(node, elem, prop)` ready for `cFSM.cfsm_analysis`.
"""
function get_cFSM_section(dimensions, material, load; n_per_segment::Int=4)

    (; t, H, B, L, D) = dimensions
    (; E, ν) = material
    (; P, Mxx, Mzz, M11, M22) = load

    # --- 8 corner nodes of the lipped hat (centerline, sharp corners) ---
    # Segment order [D, L, B, H, B, L, D]
    # Return lips point UP (free end at top), lips at bottom, web at top — matches get_section_coordinates
    x0 = t/2;  y0 = t/2
    corners_X = [x0,        x0,        x0+L,      x0+L,      x0+L+H,    x0+L+H,    x0+2L+H,   x0+2L+H]
    corners_Y = [y0+D,      y0,        y0,        y0+B,      y0+B,      y0,        y0,        y0+D    ]

    # --- interpolate n_per_segment nodes along each of the 7 segments ---
    X = Float64[]
    Y = Float64[]
    for seg in 1:7
        xi = corners_X[seg];   xj = corners_X[seg+1]
        yi = corners_Y[seg];   yj = corners_Y[seg+1]
        ts = range(0.0, 1.0; length=n_per_segment)
        # omit last node of each segment (added as first of next) except final
        lim = seg == 7 ? n_per_segment : n_per_segment - 1
        for k in 1:lim
            push!(X, xi + ts[k]*(xj-xi))
            push!(Y, yi + ts[k]*(yj-yi))
        end
    end

    nnodes = length(X)
    nelems = nnodes - 1

    # --- CUFSM node matrix  [node# x z dofx dofz dofy dofrot stress] ---
    node = zeros(nnodes, 8)
    node[:, 1] = 1:nnodes
    node[:, 2] = X
    node[:, 3] = Y
    node[:, 4:7] .= 1.0    # all DOFs free (simply-supported ends handled by BC)

    # --- CUFSM elem matrix  [elem# nodei nodej t matnum] ---
    elem = zeros(nelems, 5)
    elem[:, 1] = 1:nelems
    elem[:, 2] = 1:nelems
    elem[:, 3] = 2:nnodes
    elem[:, 4] .= t
    elem[:, 5] .= 100.0

    # --- material ---
    G    = E / (2*(1+ν))
    prop = [100.0  E  E  ν  ν  G]

    # --- section properties from the same straight-sided mesh ---
    # Must be consistent with the node/elem model used in the cFSM analysis;
    # using curved-corner properties would produce a different Ixx / yc and
    # therefore a different reference stress distribution than CUFSM software.
    straight_coord = hcat(X, Y)
    straight_ends  = hcat(1.0:(nnodes-1), 2.0:nnodes, fill(t, nnodes-1))
    sp = CUFSM.cutwp_prop2(straight_coord, straight_ends)

    # --- apply reference stress distribution via stresgen ---
    node  = CUFSM.stresgen(node, P, Mxx, Mzz, M11, M22,
                           sp.A, sp.xc, sp.yc,
                           sp.Ixx, sp.Iyy, sp.Ixy,
                           sp.θ, sp.I1, sp.I2, 0)

    return (node=node, elem=elem, prop=prop)

end


"""
    get_cFSM_section_curved(dimensions, material, load; n_per_arc=4) -> NamedTuple

Build a CUFSM node/elem/prop model using the full curved-corner geometry (arc nodes
included) for use with constrained FSM (cFSM) analysis.

Unlike `get_cFSM_section`, bend radii are modelled explicitly via arc nodes generated
by `get_section_coordinates`. Arc nodes inside each corner are treated as sub-nodes
during mode classification, which requires a relaxed collinearity tolerance passed back
as `sub_node_tol`. Pass this value to `cfsm_analysis` via the `sub_node_tol` keyword:

```julia
sec = get_cFSM_section_curved(dimensions, material, load)
result = cFSM.cfsm_analysis(sec.node, sec.elem, sec.prop, lengths;
                             modes=["L"], sub_node_tol=sec.sub_node_tol)
```

# Arguments
- `dimensions`   — `Dimensions` object
- `material`     — `Material` object
- `load`         — `Load` object (e.g. `Load(1.0,0.,0.,0.,0.)` for axial)
- `n_per_arc`    — number of arc *elements* per 90° corner (must match `n_r` used in
                   `get_section_coordinates`; default 4)

# Returns
Named tuple `(node, elem, prop, sub_node_tol)`.
- `sub_node_tol` — angle threshold (rad) to pass to `cfsm_analysis`
"""
function get_cFSM_section_curved(dimensions, material, load; n_per_arc::Int=4)

    (; E, ν) = material
    (; P, Mxx, Mzz, M11, M22) = load
    t = dimensions.t

    # Full curved-corner centerline coordinates
    coordinates = get_section_coordinates(dimensions)
    X = coordinates.X
    Y = coordinates.Y
    nnodes = length(X)
    nelems = nnodes - 1

    # --- CUFSM node matrix  [node# x z dofx dofz dofy dofrot stress] ---
    node = zeros(nnodes, 8)
    node[:, 1] = 1:nnodes
    node[:, 2] = X
    node[:, 3] = Y
    node[:, 4:7] .= 1.0

    # --- CUFSM elem matrix  [elem# nodei nodej t matnum] ---
    elem = zeros(nelems, 5)
    elem[:, 1] = 1:nelems
    elem[:, 2] = 1:nelems
    elem[:, 3] = 2:nnodes
    elem[:, 4] .= t
    elem[:, 5] .= 100.0

    # --- material ---
    G    = E / (2*(1+ν))
    prop = [100.0  E  E  ν  ν  G]

    # --- section properties from the curved geometry ---
    curved_coord = hcat(X, Y)
    curved_ends  = hcat(1.0:(nnodes-1), 2.0:nnodes, fill(t, nnodes-1))
    sp = CUFSM.cutwp_prop2(curved_coord, curved_ends)

    # --- apply reference stress distribution ---
    node = CUFSM.stresgen(node, P, Mxx, Mzz, M11, M22,
                          sp.A, sp.xc, sp.yc,
                          sp.Ixx, sp.Iyy, sp.Ixy,
                          sp.θ, sp.I1, sp.I2, 0)

    # --- sub_node_tol: arc interior nodes subtend π/(2*n_per_arc) per element.
    #     Use 1.5× that as the threshold — safely above arc angle, below 90° corners.
    sub_node_tol = 1.5 * π / (2 * n_per_arc)

    return (node=node, elem=elem, prop=prop, sub_node_tol=sub_node_tol)

end


"""
    calculate_buckling_properties(coordinates, dimensions, loads, material, lengths) -> Results

Run a CUFSM finite-strip open-section analysis over `lengths` and return the minimum
buckling load factor `Rcr` and its corresponding half-wavelength `Lcr`.

# Arguments
- `coordinates` — `(X, Y)` centerline node coordinates from `get_section_coordinates`
- `dimensions`  — `Dimensions` object
- `loads`       — `Load` object (reference load state)
- `material`    — `Material` object
- `lengths`     — vector or range of half-wavelengths to evaluate (in)
"""
function calculate_buckling_properties(coordinates, dimensions, loads, material, lengths)

    (;E,
    ν
    ) = material 

    (;P,
    Mxx,
    Mzz,
    M11,
    M22,
    ) = loads

    t = dimensions.t

    constraints = []
    springs = []
    
    neigs = 1

    model = CUFSM.Tools.open_section_analysis(coordinates.X, coordinates.Y, t, lengths, E, ν, P, Mxx, Mzz, M11, M22, constraints, springs, neigs)

    eig = 1
    Rcr_curve = CUFSM.Tools.get_load_factor(model, eig)
    Rcr = minimum(Rcr_curve) 
    
    index = argmin(Rcr_curve)
    Lcr = lengths[index]

    results = Results(model, Lcr, Rcr)

    return results 

end

"""
    calculate_Lcrd(dimensions, material, load_type) -> Float64

Compute the AISIS100 distortional buckling half-wavelength for the hat section's
bottom flange + lip assembly (`H` as web depth, `B` as flange flat width, `L` as lip).

# Arguments
- `dimensions`  — `Dimensions` object
- `material`    — `Material` object
- `load_type`   — `"P"` for axial load or `"M"` for bending moment
"""
function calculate_Lcrd(dimensions, material, load_type)

    (;t,
    H,
    B,
    L,
    D,
    r ) = dimensions

    (;E,
    ν) = material 

    CorZ = 0
    θ_top = 90.0
    #Calculate top flange + lip section properties.
    CorZ = 0

    b = B - t      # flange
    d = L - t    # lip
    θ = 90.0
    ho = H         # web depth
    μ = ν
    E = 29500.0
    G = E / (2 * (1 + μ))
    kϕ = 0.0
    Af,Jf,Ixf,Iyf,Ixyf,Cwf,xof,xhf,yhf,yof = AISIS100.v16S3.table_2_3_3__1(CorZ,t,b,d,θ)

    #Calculate the purlin distortional buckling half-wavelength.
    Lm = 999999999.0

    if load_type == "P"
        Lcrd = AISIS100.v16S3.appendix2_2_3_3_1__7(ho, μ, t, Ixf, xof, xhf, Cwf, Ixyf, Iyf)
    elseif load_type == "M"
        Lcrd = AISIS100.v16S3.appendix2_2_3_3_2__4(ho, μ, t, Ixf, xof, xhf, Cwf, Ixyf, Iyf)
    end

    return Lcrd

end



"""
    calculate_Pcrℓ(dimensions, material) -> Section
    calculate_Pcrℓ(dimensions, coordinates, material) -> Section
    calculate_Pcrℓ(dimensions, coordinates, material, lengths) -> Section

Compute the local buckling load factor under uniform axial compression (`P = 1`).
The search length range is automatically set based on `B` (flange) and `H` (web).
Returns a `Section` with `label = "Pcrℓ"`.
"""
function calculate_Pcrℓ(dimensions, material)

    label = "Pcrℓ"

    load = Load(1.0, 0.0, 0.0, 0.0, 0.0)

    B = dimensions.B
    H = dimensions.H

    if B > H
        lengths = range(0.5 * B, 1.5 * B, 9)
    else
        lengths = range(0.5 * H, 1.5 * H, 9)
    end

    coordinates = get_section_coordinates(dimensions)

    results = calculate_buckling_properties(coordinates, dimensions, load, material, lengths)

    section = Section(label, material, dimensions, load, results)

    return section 

end



function calculate_Pcrℓ(dimensions, coordinates, material)

    label = "Pcrℓ"

    load = Load(1.0, 0.0, 0.0, 0.0, 0.0)

    B = dimensions.B
    H = dimensions.H

    if B > H
        lengths = range(0.5 * B, 1.5 * B, 9)
    else
        lengths = range(0.5 * H, 1.5 * H, 9)
    end

    # coordinates = get_section_coordinates(dimensions)

    results = calculate_buckling_properties(coordinates, dimensions, load, material, lengths)

    section = Section(label, material, dimensions, load, results)

    return section 

end


function calculate_Pcrℓ(dimensions, coordinates, material, lengths)

    label = "Pcrℓ"

    load = Load(1.0, 0.0, 0.0, 0.0, 0.0)

    # B = dimensions.B 
    # H = dimensions.H

    # if B > H
    #     lengths = range(0.5 * B, 1.5 * B, 9)
    # else
    #     lengths = range(0.5 * H, 1.5 * H, 9)
    # end

    # coordinates = get_section_coordinates(dimensions)

    results = calculate_buckling_properties(coordinates, dimensions, load, material, lengths)

    section = Section(label, material, dimensions, load, results)

    return section 

end



"""
    calculate_Pcrd(dimensions, material) -> Section
    calculate_Pcrd(dimensions, coordinates, material, lengths) -> Section

Compute the distortional buckling load factor under uniform axial compression (`P = 1`).
The search length is determined via `calculate_Lcrd` with `load_type = "P"`.
Returns a `Section` with `label = "Pcrd"`.
"""
function calculate_Pcrd(dimensions, material)

    label = "Pcrd"

    load = Load(1.0, 0.0, 0.0, 0.0, 0.0)

    load_type = "P"

    Lcrd = calculate_Lcrd(dimensions, material, load_type)
    
    lengths = range(1 * Lcrd, 2 * Lcrd, 9)

    coordinates = get_section_coordinates(dimensions)

    results = calculate_buckling_properties(coordinates, dimensions, load, material, lengths)

    section = Section(label, material, dimensions, load, results)

    return section 

end



function calculate_Pcrd(dimensions, coordinates, material, lengths)

    label = "Pcrd"

    load = Load(1.0, 0.0, 0.0, 0.0, 0.0)

    # load_type = "P"
    # lengths = [calculate_Lcrd(dimensions, material, load_type)]

    # coordinates = get_section_coordinates(dimensions)

    load_type = "P"
    
    Lcrd = calculate_Lcrd(dimensions, material, load_type)
    
    lengths = range(1 * Lcrd, 2 * Lcrd, 9)

    results = calculate_buckling_properties(coordinates, dimensions, load, material, lengths)

    section = Section(label, material, dimensions, load, results)

    return section 

end




"""
    calculate_Mcrd_xx(dimensions, material) -> Section
    calculate_Mcrd_xx(dimensions, coordinates, material) -> Section
    calculate_Mcrd_xx(dimensions, coordinates, material, lengths) -> Section

Compute the distortional buckling moment factor for strong-axis bending (`Mxx = 1`).
The search length is determined via `calculate_Lcrd` with `load_type = "M"`.
Returns a `Section` with `label = "Mcrd_xx"`.
"""
function calculate_Mcrd_xx(dimensions, material)

    label = "Mcrd_xx"

    load = Load(0.0, 1.0, 0.0, 0.0, 0.0)

    load_type = "M"
    lengths = [calculate_Lcrd(dimensions, material, load_type)]

    coordinates = get_section_coordinates(dimensions)

    results = calculate_buckling_properties(coordinates, dimensions, load, material, lengths)

    section = Section(label, material, dimensions, load, results)

    return section 

end


function calculate_Mcrd_xx(dimensions, coordinates, material)

    label = "Mcrd_xx"

    load = Load(0.0, 1.0, 0.0, 0.0, 0.0)

    load_type = "M"
    lengths = [calculate_Lcrd(dimensions, material, load_type)]

    # coordinates = get_section_coordinates(dimensions)

    results = calculate_buckling_properties(coordinates, dimensions, load, material, lengths)

    section = Section(label, material, dimensions, load, results)

    return section 

end


function calculate_Mcrd_xx(dimensions, coordinates, material, lengths)

    label = "Mcrd_xx"

    load = Load(0.0, 1.0, 0.0, 0.0, 0.0)

    load_type = "M"
    # lengths = [calculate_Lcrd(dimensions, material, load_type)]

    # coordinates = get_section_coordinates(dimensions)

    results = calculate_buckling_properties(coordinates, dimensions, load, material, lengths)

    section = Section(label, material, dimensions, load, results)

    return section 

end



"""
    calculate_Mcrℓ_xx(dimensions, material) -> Section
    calculate_Mcrℓ_xx(dimensions, coordinates, material) -> Section
    calculate_Mcrℓ_xx(dimensions, coordinates, material, lengths) -> Section

Compute the local buckling moment factor for strong-axis bending (`Mxx = 1`).
The search length range is based on `B` (flange) and `H` (web).
Returns a `Section` with `label = "Mcrℓ_xx"`.
"""
function calculate_Mcrℓ_xx(dimensions, material)

    label = "Mcrℓ_xx"

    load = Load(0.0, 1.0, 0.0, 0.0, 0.0)

    B = dimensions.B
    H = dimensions.H

    if B > H / 2
        lengths = range(0.5 * B, 1.5 * B, 9)
    else
        lengths = range(0.25 * H, 1.25 * H, 9)
    end

    coordinates = get_section_coordinates(dimensions)

    results = calculate_buckling_properties(coordinates, dimensions, load, material, lengths)

    section = Section(label, material, dimensions, load, results)

    return section 

end



function calculate_Mcrℓ_xx(dimensions, coordinates, material)

    label = "Mcrℓ_xx"

    load = Load(0.0, 1.0, 0.0, 0.0, 0.0)

    B = dimensions.B
    H = dimensions.H

    if B > H / 2
        lengths = range(0.5 * B, 1.5 * B, 9)
    else
        lengths = range(0.25 * H, 1.25 * H, 9)
    end

    # coordinates = get_section_coordinates(dimensions)

    results = calculate_buckling_properties(coordinates, dimensions, load, material, lengths)

    section = Section(label, material, dimensions, load, results)

    return section 

end




function calculate_Mcrℓ_xx(dimensions, coordinates, material, lengths)

    label = "Mcrℓ_xx"

    load = Load(0.0, 1.0, 0.0, 0.0, 0.0)

    # B = dimensions.B 
    # H = dimensions.H

    # if B > H / 2
    #     lengths = range(0.5 * B, 1.5 * B, 9)
    # else
    #     lengths = range(0.25 * H, 1.25 * H, 9)
    # end

    # coordinates = get_section_coordinates(dimensions)

    results = calculate_buckling_properties(coordinates, dimensions, load, material, lengths)

    section = Section(label, material, dimensions, load, results)

    return section 

end


#this is not the traditional strong axis distortional buckling, still important though in some cases for weak axis bending
"""
    calculate_Mcrd_yy_pos(dimensions, material) -> Section
    calculate_Mcrd_yy_pos(dimensions, coordinates, material) -> Section

Compute the distortional buckling moment factor for positive weak-axis bending (`Mzz = -1`).
Not traditional strong-axis distortional buckling; relevant for weak-axis bending cases.
Returns a `Section` with `label = "Mcrd_yy_pos"`.
"""
function calculate_Mcrd_yy_pos(dimensions, material)

    label = "Mcrd_yy_pos"

    load = Load(0.0, 0.0, -1.0, 0.0, 0.0)

    load_type = "M"
    Lcrd = calculate_Lcrd(dimensions, material, load_type)
    
    lengths = range(1 * Lcrd, 2 * Lcrd, 9)

    coordinates = get_section_coordinates(dimensions)

    results = calculate_buckling_properties(coordinates, dimensions, load, material, lengths)

    section = Section(label, material, dimensions, load, results)

    return section 

end


function calculate_Mcrd_yy_pos(dimensions, coordinates, material)

    label = "Mcrd_yy_pos"

    load = Load(0.0, 0.0, -1.0, 0.0, 0.0)

    load_type = "M"
    Lcrd = calculate_Lcrd(dimensions, material, load_type)
    
    lengths = range(1 * Lcrd, 2 * Lcrd, 9)

    # coordinates = get_section_coordinates(dimensions)

    results = calculate_buckling_properties(coordinates, dimensions, load, material, lengths)

    section = Section(label, material, dimensions, load, results)

    return section 

end


"""
    calculate_Mcrℓ_yy_pos(dimensions, material) -> Section
    calculate_Mcrℓ_yy_pos(dimensions, coordinates, material) -> Section

Compute the local buckling moment factor for positive weak-axis bending (`Mzz = -1`).
Returns a `Section` with `label = "Mcrℓ_yy_pos"`.
"""
function calculate_Mcrℓ_yy_pos(dimensions, material)

    label = "Mcrℓ_yy_pos"

    load = Load(0.0, 0.0, -1.0, 0.0, 0.0)

    B = dimensions.B
    H = dimensions.H

    # lengths = range(0.5*maximum([B, H]), 2.0*maximum([B, H]), 7)

    if B > 2 * H
        lengths = range(0.5 * B, 1.5 * B, 9)
    else
        lengths = range(0.25 * H, 1.25 * H, 9)
    end

    # B = dimensions.H
    # H = dimensions.H

    # if B > 2 * H
    #     lengths = range(0.5 * B, 1.5 * B, 9)
    # else
    #     lengths = range(0.25 * H, 1.25 * H, 9)
    # end

    # coordinates = get_section_coordinates(dimensions)

    results = calculate_buckling_properties(coordinates, dimensions, load, material, lengths)

    section = Section(label, material, dimensions, load, results)

    return section 

end


function calculate_Mcrℓ_yy_pos(dimensions, coordinates, material)

    label = "Mcrℓ_yy_pos"

    load = Load(0.0, 0.0, -1.0, 0.0, 0.0)

    B = dimensions.B
    # H = dimensions.H
    L = dimensions.L

    # lengths = range(0.5*maximum([B, D]), 2.0*maximum([B, D]), 7)

    lengths = range(L, 1.2 * B, 9)

    # if B > 2 * D 
    #     lengths = range(0.5 * B, 1.5 * B, 9)
    # else
    #     lengths = range(0.25 * D, 1.25 * D, 9)
    # end


    # B = dimensions.B 
    # D = dimensions.D

    # if B > 2 * D 
    #     lengths = range(0.5 * B, 1.5 * B, 9)
    # else
    #     lengths = range(0.25 * D, 1.25 * D, 9)
    # end

    # coordinates = get_section_coordinates(dimensions)

    results = calculate_buckling_properties(coordinates, dimensions, load, material, lengths)

    section = Section(label, material, dimensions, load, results)

    return section 

end






"""
    calculate_Mcrℓ_yy_neg(dimensions, material) -> Section
    calculate_Mcrℓ_yy_neg(dimensions, coordinates, material) -> Section
    calculate_Mcrℓ_yy_neg(dimensions, coordinates, material, lengths) -> Section

Compute the local buckling moment factor for negative weak-axis bending (`Mzz = +1`).
Returns a `Section` with `label = "Mcrℓ_yy_neg"`.
"""
function calculate_Mcrℓ_yy_neg(dimensions, material)

    label = "Mcrℓ_yy_neg"

    load = Load(0.0, 0.0, 1.0, 0.0, 0.0)

    B = dimensions.B
    H = dimensions.H


    lengths = range(0.5*maximum([B, H]), 2.0*maximum([B, H]), 7)



    # if B > 2 * D
    #     lengths = range(0.5 * B, 1.5 * B, 9)
    # else
    #     lengths = range(0.25 * D, 1.25 * D, 9)
    # end

    coordinates = get_section_coordinates(dimensions)

    results = calculate_buckling_properties(coordinates, dimensions, load, material, lengths)

    section = Section(label, material, dimensions, load, results)

    return section 

end




function calculate_Mcrℓ_yy_neg(dimensions, coordinates, material)

    label = "Mcrℓ_yy_neg"

    load = Load(0.0, 0.0, 1.0, 0.0, 0.0)

    B = dimensions.B
    H = dimensions.H

    # if B > 2 * H
    #     lengths = range(0.5 * B, 1.5 * B, 9)
    # else
    #     lengths = range(0.25 * H, 1.25 * H, 9)
    # end

    lengths = range(0.5*maximum([B, H]), 2.0*maximum([B, H]), 7)



    # coordinates = get_section_coordinates(dimensions)

    results = calculate_buckling_properties(coordinates, dimensions, load, material, lengths)

    section = Section(label, material, dimensions, load, results)

    return section

end


function calculate_Mcrℓ_yy_neg(dimensions, coordinates, material, lengths)

    label = "Mcrℓ_yy_neg"

    load = Load(0.0, 0.0, 1.0, 0.0, 0.0)

    B = dimensions.B
    H = dimensions.H

    # if B > 2 * H
    #     lengths = range(0.5 * B, 1.5 * B, 9)
    # else
    #     lengths = range(0.25 * H, 1.25 * H, 9)
    # end

    # lengths = range(0.5*maximum([B, H]), 2.0*maximum([B, H]), 7)



    # coordinates = get_section_coordinates(dimensions)

    results = calculate_buckling_properties(coordinates, dimensions, load, material, lengths)

    section = Section(label, material, dimensions, load, results)

    return section 

end





"""
    flat_widths(dimensions) -> NamedTuple

Return the flat widths (excluding corner-radius arcs) for each of the 7 plate
elements of the lipped hat section in segment order `[D, L, B, H, B, L, D]`.

Tangent-point setback per 90° corner (from the virtual corner intersection):
- D–L and L–B junctions (open/outer bends): `r + t/2`  (centerline arc radius)
- B–H junctions (closed/inner bends):        `r + t`    (outside arc radius)

`H` in `dimensions` is the overall **outside** hat height (`H_inside + 2t`), so
the `r + t` setback at each B–H junction automatically yields `H_inside − 2r`.

# Returns
Named tuple with fields `D1, L1, B1, H, B2, L2, D2` (left-to-right).
"""
function flat_widths(dimensions)

    (; t, H, B, L, D, r) = dimensions

    
    r_o = r + t     


    w_D = D - r_o

    w_L = L - r_o - r_o

    w_B = B - r_o - r_o

    w_H = H - r_o - r_o

    return (D1=w_D, L1=w_L, B1=w_B, H=w_H, B2=w_B, L2=w_L, D2=w_D)

end


"""
    effective_width_stiffened(w, t, f, E, ν) -> NamedTuple

Effective width of a **stiffened** plate element per AISI S100-16 Appendix 1,
Equations 1.1-1 through 1.1-4. Uses plate-buckling coefficient `k = 4`.

A stiffened element is one supported by connected plate elements along both
longitudinal edges (e.g. web H, flanges B, lips L of the lipped hat).

# Arguments
- `w` — flat width of the element (in)
- `t` — thickness (in)
- `f` — compressive stress in the element (ksi); for compression members, `f = Fn`
- `E` — Young's modulus (ksi)
- `ν` — Poisson's ratio

# Returns
Named tuple `(b, ρ, λ, F_crl)`.
- `b`     — effective width (in)
- `ρ`     — local reduction factor (dimensionless)
- `λ`     — slenderness factor (dimensionless)
- `F_crl` — local plate-buckling stress (ksi)
"""
function effective_width_stiffened(w, t, f, E, ν)

    k = 4.0
    F_crl = k * π^2 * E / (12 * (1 - ν^2)) * (t / w)^2

    λ = sqrt(f / F_crl)

    ρ = λ <= 0.673 ? 1.0 : (1 - 0.22 / λ) / λ
    ρ = min(ρ, 1.0)

    b = ρ * w

    return (b=b, ρ=ρ, λ=λ, F_crl=F_crl)

end


"""
    effective_width_lip_stiffened(w_L, w_D_flat, D_nom, t, f, E, ν) -> NamedTuple

Effective width of a uniformly compressed element with a **simple lip edge stiffener**
per AISI S100 Section 1.3.

Applied to each lip (L) element of the lipped hat, where the return lip (D) acts as
the edge stiffener.

# Arguments
- `w_L`      — flat width of the lip element (the element being analyzed, in)
- `w_D_flat` — flat width of the return lip (stiffener), used for `Is` (in)
- `D_nom`    — nominal (actual) return lip dimension, used for `D/w` in Table 1.3-1 (in)
- `t`        — thickness (in)
- `f`        — compressive stress in the element (ksi)
- `E`        — Young's modulus (ksi)
- `ν`        — Poisson's ratio

# Procedure (θ = 90° assumed — return lip perpendicular to lip)
1. `S = 1.28 √(E/f)`  (Eq. 1.3-7)
2. If `w_L/t ≤ 0.328 S`: element is fully effective, stiffener provides full restraint.
3. Else: compute `Ia` (Eq. 1.3-8), `Is = w_D_flat³ t / 12` (Eq. 1.3-10, sin²90° = 1),
   `R₁ = min(Is/Ia, 1)` (Eq. 1.3-9), `n` (Eq. 1.3-11), `k` (Table 1.3-1).
4. Apply Section 1.1.1 with the resulting `k`.

# Returns
Named tuple `(b, ρ, λ, F_crl, k, R1, Is, Ia)`.
"""
function effective_width_partially_stiffened(w_L, w_D_flat, D_nom, t, f, E, ν)

    S   = 1.28 * sqrt(E / f)           # Eq. 1.3-7
    w_t = w_L / t

    if w_t <= 0.328 * S
        # No stiffener reduction needed — treat as if R₁ = 1
        Ia = 0.0
        Is = w_D_flat^3 * t / 12       # θ = 90°, sin²θ = 1
        R1 = 1.0
        n  = max(0.582 - w_t / (4 * S), 1/3)   # Eq. 1.3-11
        Dw = D_nom / w_L
        if Dw <= 0.25
            k = min(3.57 * R1^n + 0.43, 4.0)
        elseif Dw <= 0.8
            k = min((4.82 - 5 * Dw) * R1^n + 0.43, 4.0)
        else
            k = 4.0
        end
    else
        # Adequate moment of inertia of stiffener (Eq. 1.3-8)
        Ia = min(399 * t^4 * (w_t / S - 0.328)^3,
                 t^4 * (115 * w_t / S + 5))

        # Is: unreduced MOI of stiffener about its own centroid (Eq. 1.3-10)
        # Corner between stiffener and element is NOT part of stiffener (AISI note)
        Is = w_D_flat^3 * t / 12       # θ = 90°, sin²θ = 1

        R1 = min(Is / Ia, 1.0)         # Eq. 1.3-9

        n  = max(0.582 - w_t / (4 * S), 1/3)   # Eq. 1.3-11

        Dw = D_nom / w_L               # Table 1.3-1: D is the nominal stiffener dimension
        if Dw <= 0.25
            k = min(3.57 * R1^n + 0.43, 4.0)
        elseif Dw <= 0.8
            k = min((4.82 - 5 * Dw) * R1^n + 0.43, 4.0)
        else
            k = 4.0
        end
    end

    # Section 1.1.1 with the resulting k
    F_crl = k * π^2 * E / (12 * (1 - ν^2)) * (t / w_L)^2
    λ     = sqrt(f / F_crl)
    ρ     = λ <= 0.673 ? 1.0 : (1 - 0.22 / λ) / λ
    ρ     = min(ρ, 1.0)
    b     = ρ * w_L

    return (b=b, ρ=ρ, λ=λ, F_crl=F_crl, k=k, R1=R1, Is=Is, Ia=Ia)

end


"""
    effective_width_stiffened_web_supported(w, t, fd, E, ν, Fy) -> NamedTuple

Improved effective width for a **stiffened element supported by a web on each
longitudinal edge** per AISI S100-16 Appendix 1, Equations 1.1-6 through 1.1-8.

This improved estimate applies to elements B1, H, and B2 of the lipped hat section.

# Arguments
- `w`  — flat width of the element (in)
- `t`  — thickness (in)
- `fd` — compressive stress in the element (ksi); `fd` is substituted for `f` in λ
- `E`  — Young's modulus (ksi)
- `ν`  — Poisson's ratio
- `Fy` — yield stress (ksi)

# Returns
Named tuple `(b, ρ, λ, λ_c, F_crl)`.
- `b`     — effective width (in)
- `ρ`     — local reduction factor (dimensionless)
- `λ`     — slenderness factor (dimensionless)
- `λ_c`  — transition slenderness (Eq. 1.1-8, dimensionless)
- `F_crl` — local plate-buckling stress (ksi)
"""
function effective_width_stiffened_web_supported(w, t, fd, E, ν, Fy)

    k = 4.0
    F_crl = k * π^2 * E / (12 * (1 - ν^2)) * (t / w)^2

    λ = sqrt(fd / F_crl)

    λ_c = 0.256 + 0.328 * (w / t) * sqrt(Fy / E)  # Eq. 1.1-8

    if λ <= 0.673
        ρ = 1.0
    elseif λ < λ_c
        ρ = (1.358 - 0.461 / λ) / λ               # Eq. 1.1-6
    else
        ρ = (0.41 + 0.59 * sqrt(Fy / fd) - 0.22 / λ) / λ  # Eq. 1.1-7
    end

    ρ = min(ρ, 1.0)

    b = ρ * w

    return (b=b, ρ=ρ, λ=λ, λ_c=λ_c, F_crl=F_crl)

end


"""
    effective_width_unstiffened(w, t, f, E, ν) -> NamedTuple

Effective width of an **unstiffened** plate element per AISI S100-16 Appendix 1,
Equations 1.1-1 through 1.1-4. Uses plate-buckling coefficient `k = 0.43`.

An unstiffened element has one free longitudinal edge (e.g. the return lips D
of the lipped hat section).

# Arguments
- `w` — flat width of the element (in)
- `t` — thickness (in)
- `f` — compressive stress in the element (ksi); for compression members, `f = Fn`
- `E` — Young's modulus (ksi)
- `ν` — Poisson's ratio

# Returns
Named tuple `(b, ρ, λ, F_crl)`.
- `b`     — effective width (in)
- `ρ`     — local reduction factor (dimensionless)
- `λ`     — slenderness factor (dimensionless)
- `F_crl` — local plate-buckling stress (ksi)
"""
function effective_width_unstiffened(w, t, f, E, ν)

    k = 0.43
    F_crl = k * π^2 * E / (12 * (1 - ν^2)) * (t / w)^2

    λ = sqrt(f / F_crl)

    ρ = λ <= 0.673 ? 1.0 : (1 - 0.22 / λ) / λ
    ρ = min(ρ, 1.0)

    b = ρ * w

    return (b=b, ρ=ρ, λ=λ, F_crl=F_crl)

end


"""
    calculate_Ae(dimensions, material, f, Fy) -> NamedTuple

Compute the **effective area** `Ae` for the lipped hat section under uniform
compressive stress `f`, per AISI S100-24.

Element classification for the 7-element lipped hat `[D, L, B, H, B, L, D]`:
- D (return lips): **unstiffened** — free tip, `k = 0.43` (Section 1.2.2)
- L (lips): **partially stiffened** by the return lip D — Section 1.3 (simple lip edge
  stiffener). `k` from Table 1.3-1 based on `D/w` ratio and stiffener adequacy `R₁`.
- B1, H, B2 (flanges and web): **stiffened**, `k = 4` (Section 1.1.2, uniform stress)

Curved corners are **fully effective** per AISI S100 (they do not buckle locally).
The section has 6 × 90° corners (3 per side); each contributes area `(π/2)(r + t/2)t`.

# Arguments
- `dimensions` — `Dimensions` object
- `material`   — `Material` object
- `f`          — compressive stress in the elements (ksi), e.g. `Fn` from Chapter E
- `Fy`         — yield stress (ksi); retained for API compatibility (not used internally)

# Returns
Named tuple with fields:
- `Ae`        — effective area including fully-effective corners (in²)
- `A_corners` — total corner area (fully effective, in²)
- `widths`    — flat widths named tuple `(D1, L1, B1, H, B2, L2, D2)`
- `beff`      — effective widths named tuple
- `frcl`      — critical local buckling stress per element
- `lambda`    — slenderness per element
- `sec13`     — Section 1.3 diagnostics for lips: `(S, Ia, Is, R1, k)` for L1 and L2
"""
function calculate_Ae(dimensions, material, f, _Fy=nothing)

    (; t, r, D) = dimensions
    (; E, ν) = material

    w = flat_widths(dimensions)

    r_D1 = effective_width_unstiffened(w.D1, t, f, E, ν)
    r_L1 = effective_width_partially_stiffened(w.L1, w.D1, D, t, f, E, ν)
    r_B1 = effective_width_stiffened(w.B1, t, f, E, ν)
    r_H  = effective_width_stiffened(w.H,  t, f, E, ν)
    r_B2 = effective_width_stiffened(w.B2, t, f, E, ν)
    r_L2 = effective_width_partially_stiffened(w.L2, w.D2, D, t, f, E, ν)
    r_D2 = effective_width_unstiffened(w.D2, t, f, E, ν)

    b_total = r_D1.b + r_L1.b + r_B1.b + r_H.b + r_B2.b + r_L2.b + r_D2.b

    # 6 fully-effective 90° corners (3 per side): area = (π/2)(r + t/2)t each
    A_corners = 6 * (π/2) * (r + t/2) * t

    Ae = b_total * t + A_corners

    beff   = (D1=r_D1.b,     L1=r_L1.b,     B1=r_B1.b,     H=r_H.b,      B2=r_B2.b,     L2=r_L2.b,     D2=r_D2.b)
    frcl   = (D1=r_D1.F_crl, L1=r_L1.F_crl, B1=r_B1.F_crl, H=r_H.F_crl,  B2=r_B2.F_crl, L2=r_L2.F_crl, D2=r_D2.F_crl)
    lambda = (D1=r_D1.λ,     L1=r_L1.λ,     B1=r_B1.λ,     H=r_H.λ,      B2=r_B2.λ,     L2=r_L2.λ,     D2=r_D2.λ)

    S_lip = 1.28 * sqrt(E / f)
    sec13 = (
        L1 = (S=S_lip, Ia=r_L1.Ia, Is=r_L1.Is, R1=r_L1.R1, k=r_L1.k),
        L2 = (S=S_lip, Ia=r_L2.Ia, Is=r_L2.Is, R1=r_L2.R1, k=r_L2.k),
    )

    return (Ae=Ae, A_corners=A_corners, widths=w, beff=beff, frcl=frcl, lambda=lambda, sec13=sec13)

end



# ---------------------------------------------------------------------------
# cFSM convenience wrappers
# ---------------------------------------------------------------------------

"""
    calculate_Lcrl(dimensions, material) -> Float64

Estimate the local buckling half-wavelength for the lipped hat section.
Returns `max(B, H)` as a characteristic plate-width estimate.

# Arguments
- `dimensions`  — `Dimensions` object
- `material`    — `Material` object (reserved for future use)
"""
function calculate_Lcrl(dimensions, _material)
    (; B, H) = dimensions
    return max(B, H)
end


"""
    calculate_Pcrℓ_cFSM(dimensions, material; n_per_segment=4) -> Section
    calculate_Pcrℓ_cFSM(dimensions, material, lengths; n_per_segment=4) -> Section

Compute the local buckling load factor under uniform axial compression (`P = 1`)
using constrained FSM (cFSM), isolating the local (L) mode.
The search length range mirrors `calculate_Pcrℓ`.
Returns a `Section` with `label = "Pcrℓ_cFSM"`.
"""
function calculate_Pcrℓ_cFSM(dimensions, material; n_per_segment::Int=4)

    label = "Pcrℓ_cFSM"
    load  = Load(1.0, 0.0, 0.0, 0.0, 0.0)

    B = dimensions.B
    H = dimensions.H

    if B > H
        lengths = range(0.5 * B, 1.5 * B, 9)
    else
        lengths = range(0.5 * H, 1.5 * H, 9)
    end

    sec     = get_cFSM_section(dimensions, material, load; n_per_segment)
    result  = cFSM.cfsm_analysis(sec.node, sec.elem, sec.prop, lengths; modes=["L"])
    results = Results(result, result.Lcr, result.Rcr)

    return Section(label, material, dimensions, load, results)

end

function calculate_Pcrℓ_cFSM(dimensions, material, lengths; n_per_segment::Int=4)

    label = "Pcrℓ_cFSM"
    load  = Load(1.0, 0.0, 0.0, 0.0, 0.0)

    sec     = get_cFSM_section(dimensions, material, load; n_per_segment)
    result  = cFSM.cfsm_analysis(sec.node, sec.elem, sec.prop, lengths; modes=["L"])
    results = Results(result, result.Lcr, result.Rcr)

    return Section(label, material, dimensions, load, results)

end


"""
    calculate_Pcrd_cFSM(dimensions, material; n_per_segment=4) -> Section
    calculate_Pcrd_cFSM(dimensions, material, lengths; n_per_segment=4) -> Section

Compute the distortional buckling load factor under uniform axial compression (`P = 1`)
using constrained FSM (cFSM), isolating the distortional (D) mode.
The search length range mirrors `calculate_Pcrd`.
Returns a `Section` with `label = "Pcrd_cFSM"`.
"""
function calculate_Pcrd_cFSM(dimensions, material; n_per_segment::Int=4)

    label = "Pcrd_cFSM"
    load  = Load(1.0, 0.0, 0.0, 0.0, 0.0)

    Lcrd    = calculate_Lcrd(dimensions, material, "P")
    lengths = range(1 * Lcrd, 2 * Lcrd, 9)

    sec     = get_cFSM_section(dimensions, material, load; n_per_segment)
    result  = cFSM.cfsm_analysis(sec.node, sec.elem, sec.prop, lengths; modes=["D"])
    results = Results(result, result.Lcr, result.Rcr)

    return Section(label, material, dimensions, load, results)

end

function calculate_Pcrd_cFSM(dimensions, material, lengths; n_per_segment::Int=4)

    label = "Pcrd_cFSM"
    load  = Load(1.0, 0.0, 0.0, 0.0, 0.0)

    sec     = get_cFSM_section(dimensions, material, load; n_per_segment)
    result  = cFSM.cfsm_analysis(sec.node, sec.elem, sec.prop, lengths; modes=["D"])
    results = Results(result, result.Lcr, result.Rcr)

    return Section(label, material, dimensions, load, results)

end


"""
    calculate_Mcrℓ_xx_cFSM(dimensions, material; n_per_segment=4) -> Section
    calculate_Mcrℓ_xx_cFSM(dimensions, material, lengths; n_per_segment=4) -> Section

Compute the local buckling moment factor for strong-axis bending (`Mxx = 1`)
using constrained FSM (cFSM), isolating the local (L) mode.
The search length range mirrors `calculate_Mcrℓ_xx`.
Returns a `Section` with `label = "Mcrℓ_xx_cFSM"`.
"""
function calculate_Mcrℓ_xx_cFSM(dimensions, material; n_per_segment::Int=4)

    label = "Mcrℓ_xx_cFSM"
    load  = Load(0.0, 1.0, 0.0, 0.0, 0.0)

    B = dimensions.B
    H = dimensions.H

    if B > H / 2
        lengths = range(0.5 * B, 1.5 * B, 9)
    else
        lengths = range(0.25 * H, 1.25 * H, 9)
    end

    sec     = get_cFSM_section(dimensions, material, load; n_per_segment)
    result  = cFSM.cfsm_analysis(sec.node, sec.elem, sec.prop, lengths; modes=["L"])
    results = Results(result, result.Lcr, result.Rcr)

    return Section(label, material, dimensions, load, results)

end

function calculate_Mcrℓ_xx_cFSM(dimensions, material, lengths; n_per_segment::Int=4)

    label = "Mcrℓ_xx_cFSM"
    load  = Load(0.0, 1.0, 0.0, 0.0, 0.0)

    sec     = get_cFSM_section(dimensions, material, load; n_per_segment)
    result  = cFSM.cfsm_analysis(sec.node, sec.elem, sec.prop, lengths; modes=["L"])
    results = Results(result, result.Lcr, result.Rcr)

    return Section(label, material, dimensions, load, results)

end


"""
    calculate_Mcrd_xx_cFSM(dimensions, material; n_per_segment=4) -> Section
    calculate_Mcrd_xx_cFSM(dimensions, material, lengths; n_per_segment=4) -> Section

Compute the distortional buckling moment factor for strong-axis bending (`Mxx = 1`)
using constrained FSM (cFSM), isolating the distortional (D) mode.
The search length mirrors `calculate_Mcrd_xx` (single AISIS100 Lcrd estimate).
Returns a `Section` with `label = "Mcrd_xx_cFSM"`.
"""
function calculate_Mcrd_xx_cFSM(dimensions, material; n_per_segment::Int=4)

    label = "Mcrd_xx_cFSM"
    load  = Load(0.0, 1.0, 0.0, 0.0, 0.0)

    Lcrd    = calculate_Lcrd(dimensions, material, "M")
    lengths = [Lcrd]

    sec     = get_cFSM_section(dimensions, material, load; n_per_segment)
    result  = cFSM.cfsm_analysis(sec.node, sec.elem, sec.prop, lengths; modes=["D"])
    results = Results(result, result.Lcr, result.Rcr)

    return Section(label, material, dimensions, load, results)

end

function calculate_Mcrd_xx_cFSM(dimensions, material, lengths; n_per_segment::Int=4)

    label = "Mcrd_xx_cFSM"
    load  = Load(0.0, 1.0, 0.0, 0.0, 0.0)

    sec     = get_cFSM_section(dimensions, material, load; n_per_segment)
    result  = cFSM.cfsm_analysis(sec.node, sec.elem, sec.prop, lengths; modes=["D"])
    results = Results(result, result.Lcr, result.Rcr)

    return Section(label, material, dimensions, load, results)

end


"""
    calculate_Mcrℓ_xx_neg_cFSM(dimensions, material; n_per_segment=4) -> Section
    calculate_Mcrℓ_xx_neg_cFSM(dimensions, material, lengths; n_per_segment=4) -> Section

Compute the local buckling moment factor for negative strong-axis bending (`Mxx = -1`)
using constrained FSM (cFSM), isolating the local (L) mode.
Negative Mxx compresses the lip/flange assembly (bottom of the hat section).
The search length range mirrors `calculate_Mcrℓ_xx`.
Returns a `Section` with `label = "Mcrℓ_xx_neg_cFSM"`.
"""
function calculate_Mcrℓ_xx_neg_cFSM(dimensions, material; n_per_segment::Int=4)

    label = "Mcrℓ_xx_neg_cFSM"
    load  = Load(0.0, -1.0, 0.0, 0.0, 0.0)

    B = dimensions.B
    H = dimensions.H

    if B > H / 2
        lengths = range(0.5 * B, 1.5 * B, 9)
    else
        lengths = range(0.25 * H, 1.25 * H, 9)
    end

    sec     = get_cFSM_section(dimensions, material, load; n_per_segment)
    result  = cFSM.cfsm_analysis(sec.node, sec.elem, sec.prop, lengths; modes=["L"])
    results = Results(result, result.Lcr, result.Rcr)

    return Section(label, material, dimensions, load, results)

end

function calculate_Mcrℓ_xx_neg_cFSM(dimensions, material, lengths; n_per_segment::Int=4)

    label = "Mcrℓ_xx_neg_cFSM"
    load  = Load(0.0, -1.0, 0.0, 0.0, 0.0)

    sec     = get_cFSM_section(dimensions, material, load; n_per_segment)
    result  = cFSM.cfsm_analysis(sec.node, sec.elem, sec.prop, lengths; modes=["L"])
    results = Results(result, result.Lcr, result.Rcr)

    return Section(label, material, dimensions, load, results)

end


"""
    calculate_Mcrd_xx_neg_cFSM(dimensions, material; n_per_segment=4) -> Section
    calculate_Mcrd_xx_neg_cFSM(dimensions, material, lengths; n_per_segment=4) -> Section

Compute the distortional buckling moment factor for negative strong-axis bending (`Mxx = -1`)
using constrained FSM (cFSM), isolating the distortional (D) mode.
Negative Mxx compresses the lip/flange assembly (bottom of the hat section).
The search length mirrors `calculate_Mcrd_xx` (single AISIS100 Lcrd estimate).
Returns a `Section` with `label = "Mcrd_xx_neg_cFSM"`.
"""
function calculate_Mcrd_xx_neg_cFSM(dimensions, material; n_per_segment::Int=4)

    label = "Mcrd_xx_neg_cFSM"
    load  = Load(0.0, -1.0, 0.0, 0.0, 0.0)

    Lcrd    = calculate_Lcrd(dimensions, material, "M")
    lengths = [Lcrd]

    sec     = get_cFSM_section(dimensions, material, load; n_per_segment)
    result  = cFSM.cfsm_analysis(sec.node, sec.elem, sec.prop, lengths; modes=["D"])
    results = Results(result, result.Lcr, result.Rcr)

    return Section(label, material, dimensions, load, results)

end

function calculate_Mcrd_xx_neg_cFSM(dimensions, material, lengths; n_per_segment::Int=4)

    label = "Mcrd_xx_neg_cFSM"
    load  = Load(0.0, -1.0, 0.0, 0.0, 0.0)

    sec     = get_cFSM_section(dimensions, material, load; n_per_segment)
    result  = cFSM.cfsm_analysis(sec.node, sec.elem, sec.prop, lengths; modes=["D"])
    results = Results(result, result.Lcr, result.Rcr)

    return Section(label, material, dimensions, load, results)

end


"""
    calculate_Mcrℓ_yy_pos_cFSM(dimensions, material; n_per_segment=4) -> Section
    calculate_Mcrℓ_yy_pos_cFSM(dimensions, material, lengths; n_per_segment=4) -> Section

Compute the local buckling moment factor for positive weak-axis bending (`Mzz = -1`)
using constrained FSM (cFSM), isolating the local (L) mode.
The search length range mirrors `calculate_Mcrℓ_yy_pos`.
Returns a `Section` with `label = "Mcrℓ_yy_pos_cFSM"`.
"""
function calculate_Mcrℓ_yy_pos_cFSM(dimensions, material; n_per_segment::Int=4)

    label = "Mcrℓ_yy_pos_cFSM"
    load  = Load(0.0, 0.0, -1.0, 0.0, 0.0)

    B = dimensions.B
    L = dimensions.L
    lengths = range(L, 1.2 * B, 9)

    sec     = get_cFSM_section(dimensions, material, load; n_per_segment)
    result  = cFSM.cfsm_analysis(sec.node, sec.elem, sec.prop, lengths; modes=["L"])
    results = Results(result, result.Lcr, result.Rcr)

    return Section(label, material, dimensions, load, results)

end

function calculate_Mcrℓ_yy_pos_cFSM(dimensions, material, lengths; n_per_segment::Int=4)

    label = "Mcrℓ_yy_pos_cFSM"
    load  = Load(0.0, 0.0, -1.0, 0.0, 0.0)

    sec     = get_cFSM_section(dimensions, material, load; n_per_segment)
    result  = cFSM.cfsm_analysis(sec.node, sec.elem, sec.prop, lengths; modes=["L"])
    results = Results(result, result.Lcr, result.Rcr)

    return Section(label, material, dimensions, load, results)

end


"""
    calculate_Mcrd_yy_pos_cFSM(dimensions, material; n_per_segment=4) -> Section
    calculate_Mcrd_yy_pos_cFSM(dimensions, material, lengths; n_per_segment=4) -> Section

Compute the distortional buckling moment factor for positive weak-axis bending (`Mzz = -1`)
using constrained FSM (cFSM), isolating the distortional (D) mode.
The search length range mirrors `calculate_Mcrd_yy_pos`.
Returns a `Section` with `label = "Mcrd_yy_pos_cFSM"`.
"""
function calculate_Mcrd_yy_pos_cFSM(dimensions, material; n_per_segment::Int=4)

    label = "Mcrd_yy_pos_cFSM"
    load  = Load(0.0, 0.0, -1.0, 0.0, 0.0)

    Lcrd    = calculate_Lcrd(dimensions, material, "M")
    lengths = range(1 * Lcrd, 2 * Lcrd, 9)

    sec     = get_cFSM_section(dimensions, material, load; n_per_segment)
    result  = cFSM.cfsm_analysis(sec.node, sec.elem, sec.prop, lengths; modes=["D"])
    results = Results(result, result.Lcr, result.Rcr)

    return Section(label, material, dimensions, load, results)

end

function calculate_Mcrd_yy_pos_cFSM(dimensions, material, lengths; n_per_segment::Int=4)

    label = "Mcrd_yy_pos_cFSM"
    load  = Load(0.0, 0.0, -1.0, 0.0, 0.0)

    sec     = get_cFSM_section(dimensions, material, load; n_per_segment)
    result  = cFSM.cfsm_analysis(sec.node, sec.elem, sec.prop, lengths; modes=["D"])
    results = Results(result, result.Lcr, result.Rcr)

    return Section(label, material, dimensions, load, results)

end


"""
    calculate_Mcrℓ_yy_neg_cFSM(dimensions, material; n_per_segment=4) -> Section
    calculate_Mcrℓ_yy_neg_cFSM(dimensions, material, lengths; n_per_segment=4) -> Section

Compute the local buckling moment factor for negative weak-axis bending (`Mzz = +1`)
using constrained FSM (cFSM), isolating the local (L) mode.
The search length range mirrors `calculate_Mcrℓ_yy_neg`.
Returns a `Section` with `label = "Mcrℓ_yy_neg_cFSM"`.
"""
function calculate_Mcrℓ_yy_neg_cFSM(dimensions, material; n_per_segment::Int=4)

    label = "Mcrℓ_yy_neg_cFSM"
    load  = Load(0.0, 0.0, 1.0, 0.0, 0.0)

    B = dimensions.B
    H = dimensions.H
    lengths = range(0.5 * maximum([B, H]), 2.0 * maximum([B, H]), 7)

    sec     = get_cFSM_section(dimensions, material, load; n_per_segment)
    result  = cFSM.cfsm_analysis(sec.node, sec.elem, sec.prop, lengths; modes=["L"])
    results = Results(result, result.Lcr, result.Rcr)

    return Section(label, material, dimensions, load, results)

end

function calculate_Mcrℓ_yy_neg_cFSM(dimensions, material, lengths; n_per_segment::Int=4)

    label = "Mcrℓ_yy_neg_cFSM"
    load  = Load(0.0, 0.0, 1.0, 0.0, 0.0)

    sec     = get_cFSM_section(dimensions, material, load; n_per_segment)
    result  = cFSM.cfsm_analysis(sec.node, sec.elem, sec.prop, lengths; modes=["L"])
    results = Results(result, result.Lcr, result.Rcr)

    return Section(label, material, dimensions, load, results)

end


end # module LippedHatSectionBuckling
