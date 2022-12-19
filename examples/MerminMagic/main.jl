# Mermin-Peres Magic Square Game Algorithm
# https://www.gregegan.net/SCIENCE/MPMS/MPMS.html
using Yao

# measurement operator table
# see Table S1 in https://arxiv.org/abs/2206.12042
op_mtx = [kron(I2, Z) kron(Z, I2) kron(Z, Z); kron(X, I2) kron(I2, X) kron(X, X); kron(-X, Z) kron(-Z, X) kron(Y, Y)]

"""
    bell_state(which::Int)::AbstractRegister

Initialize bell state on a quantum register.

# Arguments
- 'which::Int': which of the four bell states to initialize
  1 corresponds to ``(|00> + |11>)/√2``
  2 corresponds to ``(|00> - |11>)/√2``
  3 corresponds to ``(|01> + |10>)/√2``
  4 corresponds to ``(|01> - |10>)/√2``
"""
function bell_state(which::Int)

    @assert 1 <= which <= 4 "Input $(which) is not in range 1..4"
    reg = zero_state(2)
    apply!(reg, chain(put(2, 1 => H), cnot(2, 1, 2)))
    (which == 2 || which == 4) && apply!(reg, put(2, 1 => Z))
    (which == 3 || which == 4) && apply!(reg, put(2, 1 => X))
    return reg
end

"""
    init_double_bell_state()::AbstractRegister

Initialize two pairs of bell state ``(|00> + |11>) ⊗ (|00> + |11>)/2``.

First qubit of each bell state pair belongs to Alice and
the others belong to Bob.
"""
function init_double_bell_state()
    return join(bell_state(1), bell_state(1))
end


"""
    meas_by_idx(reg::AbstractRegister, i::Int, j::Int)::Vector{Float64}

Perform Measurement based on index assigned by dealer.

# Arguments
- `reg`: is the quantum register to measure
- 'i::Int': the row idx assigned
- 'j::Int': the column idx assigned

"""
function meas_by_idx(reg::AbstractRegister, i::Int, j::Int)

    alice_measure!(reg::AbstractRegister, row::Int) = [real(measure!(op, reg, (1, 3))) for op in op_mtx[row, :]]
    bob_measure!(reg::AbstractRegister, col::Int) = [real(measure!(op, reg, (2, 4))) for op in op_mtx[:, col]]
    return [alice_measure!(reg, i)'..., bob_measure!(reg, j)'...]
end

"""
    input_idx(which::String)::Int

# Arguments
- 'which::String': whether user is providing index of the row or column
"""
function input_idx(which::String)
    idx = 0
    while idx > 3 || idx < 1
        println("Please provide a $which number in the range [1,3]")
        try
            idx = parse(Int64, readline())
        catch e
            println("$e ,Please enter a NUMBER!")
        end
    end
    return idx
end

function main()
    # Dealer randomly generate a pair of numbers (i,j)
    # i is the row number given to Alice
    # j is the column number given to Bob
    println("Dealer: Let's Play Mermin's Magic Square Game!")

    i = input_idx("row")
    j = input_idx("column")

    println("Alice will provide three numbers to be written on row: $(i)")
    println("Bob will provide three numbers to be written on column: $(j)")

    # initialize register to two pairs of bell state
    reg1 = init_double_bell_state()

    # do measurement and get answer string
    a1, a2, a3, b1, b2, b3 = meas_by_idx(reg1, i, j)

    println("Alice should give answer: $a1 , $a2 , $a3")
    println("Bob should give answer: $b1 , $b2 , $b3")
    println("We should see $([a1,a2,a3][j]) is equal to $([b1,b2,b3][i])")
    @assert [a1, a2, 1 * a1 * a2][j] == [b1, b2, -1 * b1 * b2][i]

    println("To reveal some magic behind the scenes, you may check that measurements on each row and each columns commute with each other and all are hermitian")
    @assert all(k -> iscommute(op_mtx[k, :]...), 1:3) && all(k -> iscommute(op_mtx[:, k]...), 1:3) && all(ishermitian, op_mtx)
end

main()
