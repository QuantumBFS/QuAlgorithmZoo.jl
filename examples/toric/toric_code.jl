using Yao, LinearAlgebra

@enum Direction RIGHT UP

function toric_code_hamiltonian(nx::Int, ny::Int, Kx::Real, Kz::Real; periodic::Bool)
    indices = vcat([(i, j, RIGHT) for i=1:nx, j=1:ny], [(i, j, UP) for i=1:nx, j=1:ny])
    n = length(indices)
    indexmap = Dict(zip(indices, 1:n))
    # Pauli-X on crosses
    return -Kx * sum(vec([repeat(n, X,
        (indexmap[(i, mod1(j-1, ny), UP)], indexmap[(i, j, RIGHT)],
        indexmap[(i, j, UP)], indexmap[(mod1(i-1, nx), j, RIGHT)]),
    ) for i=(periodic ? 1 : 2):nx, j=(periodic ? 1 : 2):ny])) -
    # Pauli-Z on squares
    Kz * sum(vec([repeat(n, Z,
            (indexmap[(i, j, RIGHT)], indexmap[(i, j, UP)],
            indexmap[(mod1(i+1, nx), j, UP)], indexmap[(i, mod1(j+1, ny), RIGHT)])
        ) for i=1:(periodic ? nx : nx-1), j=1:(periodic ? ny : ny-1)]))
end

using Test
@testset "commutativity" begin
    h = toric_code_hamiltonian(3, 3, 1, 1; periodic=true)
    for x in h[1].content, z in h[2].content
        @test iscommute(x, z)
    end
end

@testset "four fold" begin
    h2 = toric_code_hamiltonian(2, 2, 1, 1; periodic=true)
    E = eigvals(Matrix(h2))
    for i=1:4
        @test E[i] ≈ -8
    end
    @test !(E[5] ≈ -8)
end