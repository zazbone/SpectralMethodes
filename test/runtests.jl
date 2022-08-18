using Spectral
using Test

@testset "distributions" begin
    x = collect(-10:0.1:10)
    @test Spectral.gaussian(0.) ≈ 1. / sqrt(2. * pi)
    @test Spectral.cossignal([0.], [3., 5., 7.], [0.2 4. 10.]) ≈ [15.]
    @test findmax(Spectral.cossignal(x, [3., 5., 7.], [0.2, 4., 10.], true))[1] ≤ 1. 
    @test findmin(Spectral.cossignal(x, [3., 5., 7.], [0.2, 4., 10.], true))[1] ≥ -1.
    @test_throws DimensionMismatch Spectral.cossignal(x, [3., 5., 7.], [0.2])
end

@test Spectral.greet("name") == "Hello name !"
@test Spectral.hello_world() == "Hello World !"
