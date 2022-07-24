# Bernstein-Vazirani algorithm
# https://en.wikipedia.org/wiki/Bernstein%E2%80%93Vazirani_algorithm
using Yao
using YaoPlots

# build a oracle circuit that implements inner product on GF(2): ``f(x) = (-1)^{1+\sum_i x_i s_i}``
# `T` is the matrix data type,
# `s` is an unknown bit string that the program wants to find.
oracle_circuit(::Type{T}, s::BitStr{N,TI}) where {N,T,TI} = matblock(Diagonal(T[(-1)^count_ones(s & x) for x in basis(BitStr{N,TI})]); tag="oracle")

# this oracle circuit is equivalent to applying Z gates at locations i that ``s_i = 1``.
operator_fidelity(oracle_circuit(ComplexF64, bit"1101"), kron(4, Z, I2, Z, Z))

# The main circuit contains: a Hadamard transformation, oracle followed by another Hadamard transformation.
function bernstein_vazirani_circuit(::Type{T}, s::BitStr) where T
    n = length(s)
    return chain(repeat(n, H), oracle_circuit(T, s), repeat(n, H))
end
bernstein_vazirani_circuit(s::BitStr) = bernstein_vazirani_circuit(ComplexF64, s)

vizcircuit(bernstein_vazirani_circuit(bit"111"))

function main(s)
    circuit = bernstein_vazirani_circuit(s)
    return measure!(apply!(zero_state(nqubits(circuit)), circuit))
end

main(bit"0011101")
main(bit"1010101")

# TODO: implement the Hidden linear function problem
# https://en.wikipedia.org/wiki/Hidden_linear_function_problem