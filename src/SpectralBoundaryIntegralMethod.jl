module SpectralBoundaryIntegralMethod
using DrWatson
using Revise
using LinearAlgebra, FastGaussQuadrature
using AssociatedLegendrePolynomials
using GLMakie
GLMakie.activate!()

include(srcdir("thinShellMechanics.jl"))
include(srcdir("utilities.jl"))
include(srcdir("geometries.jl"))
include(srcdir("visualize.jl"))

## see srcdir("thinShellMechanics.jl")
export firstFundamentalFormCoeff, secondFundamentalFormCoeff
export unitNormalVector, derivativesUnitNormalVector
export coefficientsOfFundamentalForm
# see srcdir("geometries.jl")
export sphereGeometry
## see srcdir("utilities.jl")
export hat, upSampling
export gridOnSphere, integrationGridOnSphere
export associatedLegendreFun, derivativesAssociatedLegendreFun
export sphericalHarmonicBasisFun
export sphericalHarmonicAnalysis, sphericalHarmonicAnalysis!
export sphericalHarmonicSynthesis!
export gradient
## see srcdir("visualize.jl")
export visualizeGeometry

end
