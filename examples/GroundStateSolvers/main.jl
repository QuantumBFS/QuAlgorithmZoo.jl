using Yao
using Yao.EasyBuild: heisenberg, variational_circuit
using LinearAlgebra
using Test

include("hamiltonian_solvers.jl")

# the iterative solver (power method)
nbit = 8
h = heisenberg(nbit) |> cache
@test ishermitian(h)
reg = rand_state(nqubits(h))
E0 = expect(h, reg)/nbit
reg |> iter_groundstate!(h, niter=1000)
EG = expect(h, reg)/nbit/4
@test isapprox(EG, -0.4564, atol=1e-4)

# using imaginary time Evolution
reg = rand_state(nqubits(h))
reg |> itime_groundstate!(h, τ=20)
EG = expect(h, reg)/nbit/4
@test isapprox(EG, -0.4564, atol=1e-4)

# using VQE
N = 4
h = heisenberg(N)
E = eigen(h |> mat |> Matrix).values[1]
c = variational_circuit(N, 5)
dispatch!(c, :random)
vqe_solve!(c, h)
E2 = expect(h, zero_state(N) |> c)
@test isapprox(E, E2, rtol=1e-1)