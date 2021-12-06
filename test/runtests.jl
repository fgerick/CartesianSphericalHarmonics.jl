using CartesianSphericalHarmonics
using Test

@polyvar x y z

@testset "CartesianSphericalHarmonics.jl" begin
 
    @test rlm(1,0,x,y,z) == z
end
