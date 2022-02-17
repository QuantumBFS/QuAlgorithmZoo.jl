# several kernel functions
export rbf_kernel, brbf_kernel, rbf_functional, brbf_functional
rbf_kernel(x, y, σ::Real) = exp(-1/2σ * abs2(x-y))
rbf_kernel(x::BitStr, y::BitStr, σ::Real) = exp(-1/2σ * abs2(buffer(x-y)))
rbf_kernel(σ::Real) = (x, y) -> rbf_kernel(x, y, σ)
brbf_kernel(x, y, σ::Real) = exp(-1/2σ * count_ones(x⊻y))
brbf_kernel(σ::Real) = (x, y) -> brbf_kernel(x, y, σ)

rbf_functional(σ::Real) = StatFunctional{2}(rbf_kernel(σ))
brbf_functional(σ::Real) = StatFunctional{2}(brbf_kernel(σ))

export MMD, rbf_mmd_loss, brbf_mmd_loss
"""
    MMD{T,FT,VT<:AbstractVector{T}}
    MMD(f, probs) -> MMD

MMD loss, `VT` is the typeof probability vector, `FT` is the type of kernel function.
"""
struct MMD{T,FT,VT<:AbstractVector{T}}
    f::FT
    probs::VT
end

function Yao.expect(mmd::MMD, px::NDWeights, py::NDWeights=px)
    px_ = NDWeights(px.values .- mmd.probs)
    py_ = NDWeights(py.values .- mmd.probs)
    expect(StatFunctional{2}(mmd.f), px_, py_)
end

function (::Adjoint{Any,typeof(expect)})(mmd::MMD, circuit::Pair{<:ArrayReg, <:AbstractBlock})
    stat = StatFunctional{2}(mmd.f)
    reg, c = circuit
    out = copy(reg) |> c
    outδ = ArrayReg(witness_vec(stat, probs(out).-mmd.probs).*statevec(out))
    (in, inδ), paramsδ = YaoBlocks.AD.apply_back((out, outδ), c)
    return outδ => paramsδ.*2
end

@inline function faithful_grad(mmd::MMD, pair::Pair{<:ArrayReg,<:AbstractBlock})
    stat = StatFunctional{2}(mmd.f)
    initial = probs(copy(pair.first) |> pair.second) .- mmd.probs
    wvec = witness_vec(stat, initial)
    map(get_diffblocks(pair.second)) do diffblock
        r1, r2 = _perturb(()->_dropdims(sum(probs(copy(pair.first) |> pair.second) .* wvec, dims=1), dims=1)/2, diffblock, π/2)
        (r2 - r1)*ndims(stat)/2
    end
end

rbf_mmd_loss(σ::Real, probs) = MMD(rbf_kernel(σ), probs)
brbf_mmd_loss(σ::Real, probs) = MMD(brbf_kernel(σ), probs)
