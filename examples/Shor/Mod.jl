# TODO
# compile Mod and KMod to elementary gates.
export Mod, KMod

"""
    mod_inverse(x::Int, N::Int) -> Int

Return `y` that `(x*y)%N == 1`, notice the `(x*y)%N` operations in Z* forms a group and this is the definition of inverse.
"""
function mod_inverse(x::Int, N::Int)
    for i=1:N
        (x*i)%N == 1 && return i
    end
    throw(ArgumentError("Can not find the inverse, $x is probably not in Z*($N)!"))
end

"""
    Mod <: PrimitiveBlock

calculates `mod(a*x, L)`, notice `gcd(a, L)` should be 1.
"""
struct Mod <: PrimitiveBlock{2}
    n::Int
    a::Int
    L::Int
    function Mod(n, a, L)
        @assert gcd(a, L) == 1 && L<=1<<n
        new(n, a, L)
    end
end
Yao.nqudits(m::Mod) = m.n

function YaoBlocks.unsafe_apply!(reg::AbstractArrayReg, m::Mod)
    nstate = zero(reg.state)
    for i in basis(reg)
        _i = buffer(i) >= m.L ? buffer(i)+1 : mod(buffer(i)*m.a, m.L)+1
        for j in 1:size(nstate, 2)
            @inbounds nstate[_i,j] = reg.state[buffer(i)+1,j]
        end
    end
    reg.state = nstate
    reg
end

function Yao.mat(::Type{T}, m::Mod) where {T}
    perm = Vector{Int}(undef, 1<<m.n)
    for i in Yao.basis(m.n)
        @inbounds perm[i >= m.L ? i+1 : mod(i*m.a, m.L)+1] = i+1
    end
    PermMatrix(perm, ones(T, 1<<m.n))
end

Base.adjoint(m::Mod) = Mod(nqubits(m), mod_inverse(m.a, m.L), m.L)
Yao.print_block(io::IO, m::Mod) = print(io, "Mod: $(m.a)*x % $(m.L) (nqubits = $(nqubits(m)))")

Yao.isunitary(::Mod) = true
# Yao.ishermitian(::Mod) = true  # this is not true for L = 1.
# Yao.isreflexive(::Mod) = true  # this is not true for L = 1.

"""
    KMod <: PrimitiveBlock{2}

The first `k` qubits are exponent, and the rest `n-k` are base `a`,
it calculates `mod(a^k*x, L)`, notice `gcd(a, L)` should be 1.
"""
struct KMod <: PrimitiveBlock{2}
    n::Int
    k::Int
    a::Int
    L::Int
    function KMod(n, k, a, L)
        @assert gcd(a, L) == 1 && L<=1<<(n-k)
        new(n, k, a, L)
    end
end

Yao.nqudits(m::KMod) = m.n
nqubits_data(m::KMod) = nqubits(m)-m.k

function bint2_reader(T, k::Int)
    mask = bmask(T, 1:k)
    return b -> (b&mask, b>>k)
end

function YaoBlocks.unsafe_apply!(reg::ArrayReg, m::KMod)
    nstate = zero(reg.state)

    reader = bint2_reader(Int, m.k)
    for b in Yao.basis(reg)
        k, i = reader(buffer(b))
        _i = i >= m.L ? i : mod(i*powermod(m.a, k, m.L), m.L)
        _b = k + _i<<m.k + 1
        for j in 1:size(nstate, 2)
            @inbounds nstate[_b,j] = reg.state[buffer(b)+1,j]
        end
    end
    reg.state = nstate
    reg
end

function Yao.mat(::Type{T}, m::KMod) where {T}
    perm = Vector{Int}(undef, 1<<m.n)
    reader = bint2_reader(Int, m.k)
    for b in Yao.basis(m.n)
        k, i = reader(buffer(b))
        _i = i >= m.L ? i : mod(i*powermod(m.a, k, m.L), m.L)
        _b = k + _i<<m.k + 1
        @inbounds perm[_b] = buffer(b)+1
    end
    YaoBlocks.LuxurySparse.PermMatrix(perm, ones(T, 1<<m.n))
end

Base.adjoint(m::KMod) = KMod(m.n, m.k, mod_inverse(m.a, m.L), m.L)
Yao.print_block(io::IO, m::KMod) = print(io, "Mod: $(m.a)^k*x % $(m.L) (nqubits = $(nqudits(m)), number of control bits = $(m.k))")

Yao.isunitary(::KMod) = true
# Yao.ishermitian(::Mod) = true  # this is not true for L = 1.
# Yao.isreflexive(::Mod) = true  # this is not true for L = 1.
