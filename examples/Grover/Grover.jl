# # [`Grover Search`](@id Grover)
using Yao
using LinearAlgebra

"""
A simple inference oracle, e.g. inference([-1, -8, 5]) is a control block that flip the bit if values of bits on position [1, 8, 5] match [0, 0, 1].
"""
inference_oracle(nbit::Int, locs::Vector{Int}) = control(nbit, locs[1:end-1], abs(locs[end]) => (locs[end]>0 ? Z : chain(phase(π), Z)))

# Compute the propotion of target states to estimate the number of iterations.
reg = ArrayReg(ones(ComplexF64, 1<<nbit))
ratio = norm(statevec(reg)[real(statevec(reg |> oracle)) .< 0])^2
num_grover_step = Int(round(pi/4/sqrt(ratio)))-1

"""
    GroverIter{N}

    GroverIter(oracle, ref::ReflectBlock{N}, psi::ArrayReg, niter::Int)

an iterator that perform Grover operations step by step.
An Grover operation consists of applying oracle and Reflection.
"""
struct GroverIter{N}
    psi::ArrayReg
    oracle
    refl::ReflectBlock{N}
    niter::Int
end

groveriter(psi::ArrayReg, oracle, refl::ReflectBlock{N}, niter::Int) where {N} = GroverIter{N}(psi, oracle, ref, niter)

function Base.iterate(it::GroverIter, st=1)
    if it.niter + 1 == st
        nothing
    else
        apply!(it.psi, it.oracle)
        apply!(it.psi, it.ref), st+1
    end
end

"""
    groverblock(oracle, ref::ReflectBlock{N}, niter::Int=-1)
    groverblock(oracle, psi::ArrayReg, niter::Int=-1)

Return a ChainBlock/Sequential as Grover Iteration, the default `niter` will stop at the first optimal step.
"""
function groverblock(oracle::AbstractBlock{N}, ref::ReflectBlock{N}, niter::Int=-1) where {N}
    if niter == -1 niter = num_grover_step(ref.psi, oracle) end
    chain(N, chain(oracle, ref) for i = 1:niter)
end

function groverblock(oracle, ref::ReflectBlock{N}, niter::Int=-1) where {N}
    if niter == -1 niter = num_grover_step(ref.psi, oracle) end
    sequence(sequence(oracle, ref) for i = 1:niter)
end

groverblock(oracle, psi::ArrayReg, niter::Int=-1) = groverblock(oracle, ReflectBlock(psi |> copy), niter)
# ## Target Space and Evidense
num_bit = 12
oracle = matblock(Diagonal((v = ones(ComplexF64, 1<<num_bit); v[100:101]*=-1; v)))
target_state = zeros(1<<num_bit); target_state[100:101] .= sqrt(0.5)

# now we want to search the subspace with [1,3,5,8,9,11,12] fixed to 1 and [4,6] fixed to 0.
evidense = [1, 3, -4, 5, -6, 8, 9, 11, 12]

# ## Search
# then solve the above problem
it = groveriter(uniform_state(num_bit), oracle)
for (i, reg) in enumerate(it)
    overlap = abs(statevec(reg)'*target_state)
    println("step $(i-1), overlap = $overlap")
end

# ## Inference Example
# we have a state psi0, we know how to prepair it
reg = rand_state(num_bit)

"""
Inference: reg is the initial state, we want to search target space with specific evidense.
e.g. evidense [1, -3, 6] means the [1, 3, 6]-th bits take value [1, 0, 1].
"""
oracle_infer = inference_oracle(evidense)(nqubits(reg))
it = groveriter(reg, oracle_infer)
for (i, reg) in enumerate(it)
    p_target = prob_match_oracle(reg, oracle_infer)
    println("step $(i-1), overlap^2 = $p_target")
end
