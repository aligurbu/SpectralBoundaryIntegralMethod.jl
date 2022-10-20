using Test, SafeTestsets

@testset "SpectralBoundaryIntegralMethod.jl" begin
    # Write your tests here.
end

@safetestset "This is only a test" begin
    include("ThisIsOnlyATest.jl")
end