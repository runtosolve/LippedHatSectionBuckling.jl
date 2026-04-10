include("../src/LippedHatSectionBuckling.jl")
using CUFSM

E = 29500.0
ν = 0.30 

t = 0.0188
D = 0.276
L = 0.551

B = 2.01
H = 1.54
h = H + t*2
r = 0.0790


material = LippedHatSectionBuckling.Material(E, ν)
dimensions = LippedHatSectionBuckling.Dimensions(t, h, B, L, D, r)

section = LippedHatSectionBuckling.calculate_Pcrd(dimensions, material)

Lcrd = LippedHatSectionBuckling.calculate_Lcrd(dimensions, material, "P")
dump(section)

nodes = section.results.model.node
elements = section.results.model.elem


coordinates = LippedHatSectionBuckling.get_section_coordinates(dimensions)

coord = hcat(coordinates.X, coordinates.Y)   # Nx2 matrix

n = length(coordinates.X)
ends = hcat(1:n-1, 2:n, fill(dimensions.t, n-1), ones(n-1))

section_properties = CUFSM.cutwp_prop2(coord, ends)
A = section_properties.A
