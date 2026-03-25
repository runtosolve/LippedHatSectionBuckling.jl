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


using CrossSectionGeometry, CUFSM, AISIS100


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
    n_r = [5, 5, 5, 5, 5, 5];
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
    lengths = [calculate_Lcrd(dimensions, material, load_type)]

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
    
    lengths = range(0.5 * Lcrd, 1.5 * Lcrd, 9)

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
    
    lengths = range(0.5 * Lcrd, 1.5 * Lcrd, 9)

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




end # module LippedHatSectionBuckling
