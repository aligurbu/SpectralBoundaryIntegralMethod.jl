using DrWatson
quickactivate(@__DIR__, "SpectralBoundaryIntegralMethod")
using SpectralBoundaryIntegralMethod
using Test, SafeTestsets

@testset "SpectralBoundaryIntegralMethod.jl" begin
    @test 2+3 == 5
end

@safetestset "Testing the Gauss-Legendre quadrature" begin
    include("FastGaussQuadrature.jl")
end

@safetestset "Testing the associated Legendre functions" begin
    include("AssociatedLegendrePolynomials.jl")
end