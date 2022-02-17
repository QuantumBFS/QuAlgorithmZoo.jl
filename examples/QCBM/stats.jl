using LinearAlgebra: Adjoint
export as_weights
export StatFunctional

"""
    StatFunctional{N, F}
    StatFunctional{N}(f::F) -> StatFunctional{N}

`f` is a function, `f(xᵢ,xⱼ,xₖ...)`, this functional is `1/C(r,n)... ∑ᵢⱼₖ...f(xᵢ,xⱼ,xₖ...)`, see U-statistics for detail.

References:
    U-statistics, http://personal.psu.edu/drh20/asymp/fall2006/lectures/ANGELchpt10.pdf
"""
struct StatFunctional{N, F}
    f::F
    StatFunctional{N}(f::F) where {N, F} = new{N, F}(f)
end
Base.ndims(stat::StatFunctional{N}) where N = N

"""
    NDWeights{N, AT<:AbstractArray{T,N} where T}

Weights not limited to 1 dimension.
"""
struct NDWeights{N, AT<:AbstractArray{T,N} where T}
    values::AT
end
Base.getindex(w::NDWeights, args...) = getindex(w.values, args...)
Base.length(w::NDWeights) = length(w.values)
Base.size(w::NDWeights, args...) = size(w.values, args...)

"U-statistics of order 2."
function Yao.expect(stat::StatFunctional{2}, xs::AbstractVecOrMat{T}) where T<:Integer
    N = size(xs, 1)
    res = map(1:size(xs, 2)) do b
        res = zero(stat.f(xs[1], xs[1]))
        for i = 2:N
            for j = 1:i-1
                @inbounds res += stat.f(xs[i,b], xs[j,b])
            end
        end
        res/binomial(N,2)
    end
    length(res) == 1 ? res[] : res
end

function Yao.expect(stat::StatFunctional{2}, xs::AbstractVecOrMat{T}, ys::AbstractVecOrMat{T}) where T<:Integer
    ci = CartesianIndices((size(xs,1), size(ys,1)))
    res = map(1:size(xs, 2)) do b
        @inbounds mapreduce(ind->stat.f(xs[ind[1], b], ys[ind[2], b]), +, ci)/length(ci)
    end
    length(res) == 1 ? res[] : res
end

function Yao.expect(stat::StatFunctional{2}, px::NDWeights, py::NDWeights=px)
    Tx = BitStr64{log2dim1(px)}
    Ty = BitStr64{log2dim1(py)}
    ci = CartesianIndices((size(px,1), size(py,1)))
    res = [@inbounds mapreduce(ind->px[ind[1], b]*stat.f(Tx(ind[1]-1), Ty(ind[2]-1))*py[ind[2], b], +, ci) for b in 1:size(px, 2)]
    length(res) == 1 ? res[] : res
end

function Yao.expect(stat::StatFunctional{1}, xs::AbstractVecOrMat{<:Integer})
    res = mean(stat.f.(xs), dims=1)
    _dropdims(res, dims=1)
end

function Yao.expect(stat::StatFunctional{1}, px::NDWeights)
    T = BitStr64{log2dim1(px)}
    res = [mapreduce(i->stat.f(T(i-1)) * px[i,b], +, 1:size(px,1)) for b in 1:size(px, 2)]
    length(res) == 1 ? res[] : res
end

as_weights(probs::AbstractArray) = NDWeights(probs)
as_weights(reg::ArrayReg) = reg |> probs |> as_weights

@inline function faithful_grad(stat::StatFunctional{2}, pair::Pair{<:ArrayReg,<:AbstractBlock})
    initial = copy(pair.first) |> pair.second |> as_weights
    map(get_diffblocks(pair.second)) do diffblock
        r1, r2 = _perturb(()->expect(stat, copy(pair.first) |> pair.second |> as_weights, initial), diffblock, π/2)
        (r2 - r1)*ndims(stat)/2
    end
end

@inline function faithful_grad(stat::StatFunctional{1}, pair::Pair{<:ArrayReg,<:AbstractBlock})
    map(get_diffblocks(pair.second)) do diffblock
        r1, r2 = _perturb(()->expect(stat, copy(pair.first) |> pair.second |> as_weights), diffblock, π/2)
        (r2 - r1)*ndims(stat)/2
    end
end

function (::Adjoint{Any,typeof(expect)})(stat::StatFunctional, circuit::Pair{<:ArrayReg, <:AbstractBlock})
    reg, c = circuit
    out = copy(reg) |> c
    outδ = ArrayReg(witness_vec(stat, out |> probs).*statevec(out))
    (in, inδ), paramsδ = apply_back((out, outδ), c)
    return outδ => paramsδ.*2
end

function witness_vec(stat::StatFunctional{2}, probs::AbstractVecOrMat)
    T = BitStr64{log2dim1(probs)}
    Nb = size(probs,2)
    Tv = typeof(stat.f(T(0), T(0)))
    res = zeros(Tv, size(probs)...)
    for b = 1:Nb
        res[:,b] = map(i->2*mapreduce(j->stat.f(i, T(j-1))*probs[j,b], +, 1:size(probs,1)), basis(T))
    end
    return res
end

function witness_vec(stat::StatFunctional{1}, probs::AbstractVecOrMat)
    res = map(stat.f, basis(BitStr64{log2dim1(probs)}))
    if ndims(probs) == 2
        return repeat(res, 1, size(probs,2))
    end
    return res
end
