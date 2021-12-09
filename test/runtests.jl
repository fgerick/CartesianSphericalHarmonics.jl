using CartesianSphericalHarmonics
using Test

@polyvar x y z

∂ = differentiate
Δ(f) = ∂(f,x,2) + ∂(f,y,2) + ∂(f,z,2)

@testset "CartesianSphericalHarmonics.jl" begin
 
    @test rlm(1,0,x,y,z) == z
    @test rlm(1,1,x,y,z) == -x

    #https://en.wikipedia.org/wiki/Table_of_spherical_harmonics
    @test rlm(0,0,x,y,z; real=false, norm=Laplace{Float64}) ≈ 1/2*sqrt(1/π)
    @test rlm(1,1,x,y,z; real=false, norm=Laplace{Float64}) ≈ -1/2*sqrt(3/2π)*(x+im*y)
    @test rlm(1,1,x,y,z; real=false, norm=Laplace{Float64}) ≈ -1/2*sqrt(3/2π)*(x+im*y)
    @test rlm(2,0,x,y,z; real=false, norm=Laplace{Float64} ) ≈  1/4*sqrt(5/pi)*(3z^2-x^2-y^2-z^2)
    @test rlm(2,2,x,y,z; real=false, norm=Laplace{Float64} ) ≈  1/4*sqrt(15/2pi)*(x+im*y)^2
    @test rlm(2,1,x,y,z; real=false, norm=Laplace{Float64} ) ≈  -1/2*sqrt(15/2pi)*(x+im*y)*z
    
    @test_throws DomainError rlm(1,2,x,y,z)

    @test Δ(rlm(5,2,x,y,z)) == 0
end
