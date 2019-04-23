module QuAlgorithmZoo

using LuxurySparse, LinearAlgebra
using MacroTools: @forward
using Yao, Yao.ConstGate, BitBasis
using YaoArrayRegister: u1rows!
import Yao: mat, dispatch!, niparams, getiparams, setiparams!, cache_key, print_block, apply!, PrimitiveBlock, ishermitian, isunitary, isreflexive
import Base: ==, copy, hash

export openbox

"""
    openbox(block::AbstractBlock) -> AbstractBlock

For a black box, like QFTBlock, you can get its white box (loyal simulation) using this function.
"""
function openbox end

#=
export timeevolve, join, DefaultRegister, MatrixBlock, Sequential, ReflectBlock, GeneralMatrixBlock, AddBlock, sequence, matrixgate
export ⊗
import Base: join
const DefaultRegister = ArrayReg
const MatrixBlock = AbstractBlock
const Sequential = ChainBlock
const ReflectBlock = ReflectGate
const GeneralMatrixGate = GeneralMatrixBlock
const AddBlock = Sum
sequence(args...) = chain(args...)
matrixgate(args...) = matblock(args...)
Yao.shift(θ::Real) = shift(Float64(θ))
Yao.phase(θ::Real) = phase(Float64(θ))
timeevolve(args...) = TimeEvolution(args...)
=#
#=
function join(reg1::DefaultRegister{B, T1}, reg2::DefaultRegister{B, T2}) where {B, T1, T2}
    s1 = reg1 |> rank3
    s2 = reg2 |> rank3
    T = promote_type(T1, T2)
    state = Array{T,3}(undef, size(s1, 1)*size(s2, 1), size(s1, 2)*size(s2, 2), B)
    for b = 1:B
        @inbounds @views state[:,:,b] = kron(s1[:,:,b], s2[:,:,b])
    end
    DefaultRegister{B}(reshape(state, size(state, 1), :))
end
join(reg1::DefaultRegister{1}, reg2::DefaultRegister{1}) = DefaultRegister{1}(kron(reg1.state, reg2.state))
join(A::AbstractRegister, B::AbstractRegister) = cat(A, B)
# joining two registers
⊗(reg::AbstractRegister, reg2::AbstractRegister) = join(reg, reg2)
⊗(A::AbstractArray, B::AbstractArray) = kron(A, B)
=#

include("Miscellaneous.jl")
include("Diff.jl")
include("Adam.jl")
include("QFT.jl")
include("CircuitBuild.jl")
include("QCOptProblem.jl")
include("RotBasis.jl")
include("Grover.jl")
include("PhaseEstimation.jl")
include("HHL.jl")
include("hamiltonian_solvers.jl")
include("HadamardTest.jl")
include("lin_diffEq_HHL.jl")


end # module
