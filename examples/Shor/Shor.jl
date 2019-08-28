using Yao, YaoExtensions, BitBasis
using Test, Random
include("number_theory.jl")
include("shorlib.jl")

@testset "shor_classical" begin
    Random.seed!(129)
    L = 35
    f = shor(L, Val(:classical))
    @test f == 5 || f == 7

    L = 25
    f = shor(L, Val(:classical))
    @test f == 5

    L = 7*11
    f = shor(L, Val(:classical))
    @test f == 7 || f == 11

    L = 14
    f = shor(L, Val(:classical))
    @test f == 2 || f == 7

    @test NumberTheory.factor_a_power_b(25) == (5, 2)
    @test NumberTheory.factor_a_power_b(15) == nothing
end

@testset "shor quantum" begin
    Random.seed!(129)
    L = 15
    f = shor(L, Val(:quantum))
    @test f == 5 || f == 3
end
