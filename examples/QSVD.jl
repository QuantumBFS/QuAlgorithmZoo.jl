# Quantum SVD
using Yao
using Random
using QuAlgorithmZoo
using LinearAlgebra

"""
Quantum singular value decomposition algorithm.

    * `reg`, input register (A, B) as the target matrix to decompose,
    * `circuit_a`, U matrix applied on register A,
    * `circuit_b`, V matrix applied on register B,
"""
function train!(reg, circuit_a::AbstractBlock{Na}, circuit_b::AbstractBlock{Nb}, optimizer; maxiter::Int=100) where {Na, Nb}
    nbit = Na+Nb
    cnots = chain(control(nbit, i+Na, i=>X) for i=1:min(Na, Nb))
    c = chain(concentrate(nbit, circuit_a, 1:Na), concentrate(nbit, circuit_b, Na+1:nbit), cnots)
    c = c |> autodiff(:QC)   # construct a differentiable circuit for training
    println(c)

    obs = -mapreduce(i->put(nbit, i=>Z), +, 1:Na)
    params = parameters(c)
    for i = 1:maxiter
        grad = opdiff.(() -> copy(reg) |> c, collect_blocks(AbstractDiff, c), Ref(obs))
        QuAlgorithmZoo.update!(params, grad, optimizer)
        @show expect(obs, copy(reg) |> c)
        dispatch!(c, params)
    end
end

# define a matrix of size (2^Na, 2^Nb)
Na = 2
Nb = 2
nbit = Na + Nb
reg = rand_state(nbit)

# the exact result
M = reshape(reg.state, 1<<Na, 1<<Nb)
U_exact, S_exact, V_exact = svd(M)

circuit_a = random_diff_circuit(Na, 5, pair_ring(Na))
circuit_b = random_diff_circuit(Nb, 5, pair_ring(Nb))

Random.seed!(2)
dispatch!(circuit_a, :random)
dispatch!(circuit_b, :random)
train!(reg, circuit_a, circuit_b, Adam(lr=0.1))

order = [1,4,3,2]  # different up to an order
display(mat(circuit_a)[order,:]*U_exact .|> abs2)
display(conj(mat(circuit_b)[order,:])*V_exact .|> abs2)
