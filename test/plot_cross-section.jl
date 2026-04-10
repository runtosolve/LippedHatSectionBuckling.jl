include("../src/LippedHatSectionBuckling.jl")
using CairoMakie
using SectionProperties


E = 29500.0
ν = 0.30 

t = 0.0188
D = 0.276
L = 0.551

B = 2.01
H = 1.54 + t*2
r = 0.0790


material = LippedHatSectionBuckling.Material(E, ν)
dimensions = LippedHatSectionBuckling.Dimensions(t, H, B, L, D, r)

output = LippedHatSectionBuckling.get_section_coordinates(dimensions)

X = output.X
Y = output.Y

fig, ax, _ = scatterlines(X, Y, markersize=4, axis=(; aspect=DataAspect()))

cs = [[output.X[j], output.Y[j]] for j in eachindex(output.X)]
properties = SectionProperties.open_thin_walled(cs, t * ones(Float64, length(cs) - 1))

properties.A
properties.Ixx
properties.Iyy


Lcrd = LippedHatSectionBuckling.calculate_Lcrd(dimensions, material, "P")



coordinates = LippedHatSectionBuckling.get_section_coordinates(dimensions)


function export_coordinates_csv(coordinates, filename)

    open(filename, "w") do f
        println(f, "node,x,y")
        for i in eachindex(coordinates.X)
            println(f, i, ",", coordinates.X[i], ",", coordinates.Y[i])
        end
    end

end

export_coordinates_csv(coordinates, "coordinates.csv")