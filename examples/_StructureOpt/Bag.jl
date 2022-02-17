export AbstractBag
"""
    AbstractBag{BT, D}<:TagBlock{BT, D}

Abstract `Bag` is a wrapper of a block that conserves all properties.
Including `mat`, `apply!`, `ishermitian`, `isreflexive`, `isunitary`,
`occupied_locs`, `apply_back!` and `mat_back!`.
"""
abstract type AbstractBag{BT, D}<:TagBlock{BT, D} end

Yao.mat(::Type{T}, bag::AbstractBag) where {T} = mat(T, content(bag))
_apply!(reg::AbstractRegister, bag::AbstractBag) = _apply!(reg, content(bag))
Yao.ishermitian(bag::AbstractBag) = ishermitian(content(bag))
Yao.isreflexive(bag::AbstractBag) = isreflexive(content(bag))
Yao.isunitary(bag::AbstractBag) = isunitary(content(bag))
Yao.occupied_locs(bag::AbstractBag) = occupied_locs(content(bag))

function Yao.AD.apply_back!(state, b::AbstractBag, collector)
    Yao.AD.apply_back!(state, content(b), collector)
end
function Yao.AD.mat_back!(::Type{T}, b::AbstractBag, adjy, collector) where T
    Yao.AD.mat_back!(T, content(b), adjy, collector)
end

export Bag, enable_block!, disable_block!, setcontent!, isenabled
"""
    Bag{D}<:TagBlock{AbstractBlock, D}

A bag is a trivil container, but can
    * `setcontent!(bag, content)`
    * `disable_block!(bag)`
    * `enable_block!(bag)`
"""
mutable struct Bag{D}<:AbstractBag{AbstractBlock, D}
    content::AbstractBlock{D}
    mask::Bool
end
Bag(b::AbstractBlock) = Bag(b,true)

Yao.content(bag::Bag{D}) where {D} = bag.mask ? bag.content : put(nqudits(bag), 1=>IdentityGate{D}(nqudits(bag)))
Yao.chcontent(bag::Bag, content) = Bag(content)
setcontent!(bag::Bag, content) = (bag.content = content; bag)
disable_block!(b::Bag) = (b.mask = false; b)
enable_block!(b::Bag) = (b.mask = true; b)
isenabled(b::Bag) = b.mask

function Yao.print_annotation(io::IO, bag::Bag)
    printstyled(io, isenabled(bag) ? "[⊙] " : "[⊗] "; bold=true, color=isenabled(bag) ? :green : :red)
end

function Base.show(io::IO, ::MIME"plain/text", blk::Bag)
    return print_tree(io, blk; title=false, compact=false)
end
