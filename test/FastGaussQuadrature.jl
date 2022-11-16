using LinearAlgebra, FastGaussQuadrature

@testset "FastGaussQuadrature.jl" begin
    x, w = gausslegendre(3)
    @test sum(x.^2 .* w) ≈ 2/3
    @test dot(x.^2, w) ≈ 2/3
end