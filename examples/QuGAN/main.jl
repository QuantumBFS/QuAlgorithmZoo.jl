using Yao
using Yao.EasyBuild: variational_circuit
include("../common/Adam.jl")
using .SimpleOptimizers: Adam, update!
import Yao: tracedist

"""
Quantum GAN.

Reference:
    Benedetti, M., Grant, E., Wossnig, L., & Severini, S. (2018). Adversarial quantum circuit learning for pure state approximation, 1–14.
"""
struct QuGAN
    nqubits::Int
    target::ArrayReg
    generator::AbstractBlock
    discriminator::AbstractBlock
    reg0::ArrayReg
    witness_op::AbstractBlock
    circuit::AbstractBlock

    function QuGAN(target::ArrayReg, gen::AbstractBlock, dis::AbstractBlock)
        N = nqubits(target)
        c = chain(subroutine(N+1, gen, 1:N), dis)
        witness_op = put(N+1, (N+1)=>ConstGate.P0)
        new(N+1, join(zero_state(1), target), subroutine(N+1, gen, 1:N), dis, zero_state(N+1), witness_op, c)
    end
end

# INTERFACES
circuit(qg::QuGAN) = qg.circuit
loss(qg::QuGAN) = p0t(qg) - p0g(qg)

function gradient(qg::QuGAN)
    grad_gen = expect'(qg.witness_op, qg.reg0 => qg.circuit).second
    grad_tar = expect'(qg.witness_op, qg.target => qg.circuit[2]).second
    ngen = nparameters(qg.generator)
    [-grad_gen[1:ngen]; grad_tar - grad_gen[ngen+1:end]]
end

"""probability to get evidense qubit 0 on generation set."""
p0g(qg::QuGAN) = expect(qg.witness_op, qg.reg0 => qg.circuit) |> real
"""probability to get evidense qubit 0 on target set."""
p0t(qg::QuGAN) = expect(qg.witness_op, qg.target => qg.circuit[2]) |> real
"""generated wave function"""
outputψ(qg::QuGAN) = copy(qg.reg0) |> qg.generator

"""tracedistance between target and generated wave function"""
tracedist(qg::QuGAN) = tracedist(qg.target, outputψ(qg))[]

using Test, Random
Random.seed!(2)

nbit = 3
depth_gen = 4
depth_dis = 4

# define a QuGAN
target = rand_state(nbit)
generator = dispatch!(variational_circuit(nbit, depth_gen), :random)
discriminator = dispatch!(variational_circuit(nbit+1, depth_dis), :random)
qg = QuGAN(target, generator, discriminator)

# check the gradient
grad = gradient(qg)

# learning rates for the generator and discriminator
g_lr = 0.2
d_lr = 0.5
for i=1:300
    ng = nparameters(qg.generator)
    grad = gradient(qg)
    dispatch!(-, qg.generator, grad[1:ng]*g_lr)
    dispatch!(-, qg.discriminator, -grad[ng+1:end]*d_lr)
    println("Step $i, trace distance = $(tracedist(qg))")
end

@test qg |> loss < 0.1
