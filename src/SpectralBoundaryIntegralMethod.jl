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

export getSphereGeometry        # see srcdir("geometries.jl")
export hat, gridOnSphere        # see srcdir("utilities.jl")
export associatedLegendreFun    # see srcdir("utilities.jl")
export visualizeGeometry        # see srcdir("visualize.jl")

end
