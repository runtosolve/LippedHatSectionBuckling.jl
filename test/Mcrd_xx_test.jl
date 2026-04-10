include("../src/LippedHatSectionBuckling.jl")

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

Mcrd = LippedHatSectionBuckling.calculate_Mcrd_xx(dimensions, material)



