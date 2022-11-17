module SpectralBoundaryIntegralMethod
using DrWatson
using Revise
using LinearAlgebra, FastGaussQuadrature
using AssociatedLegendrePolynomials
using GLMakie
GLMakie.activate!()

include(srcdir("utilities.jl"))
include(srcdir("geometries.jl"))
include(srcdir("visualize.jl"))

# see srcdir("geometries.jl")
export getSphereGeometry
## see srcdir("utilities.jl")
export hat, gridOnSphere
export associatedLegendreFun
export sphericalHarmonicBasisFun
export sphericalHarmonicAnalysis
## see srcdir("visualize.jl")
export visualizeGeometry

end
