# Mermin-Peres Magic Square Game Algorithm
# https://www.gregegan.net/SCIENCE/MPMS/MPMS.html
using Yao

# measurement operator table
# see Table S1 in https://arxiv.org/abs/2206.12042
const magic_square = [kron(I2, Z) kron(Z, I2) kron(Z, Z)
    kron(X, I2) kron(I2, X) kron(X, X)
    kron(-X, Z) kron(-Z, X) kron(Y, Y)]

"""
    bell_state(which::Int) -> AbstractRegister

Initialize bell state on a quantum register.

# Arguments
- `which::Int`: which of the four bell states to initialize
  1 corresponds to ``(|00> + |11>)/√2``
  2 corresponds to ``(|00> - |11>)/√2``
  3 corresponds to ``(|01> + |10>)/√2``
  4 corresponds to ``(|01> - |10>)/√2``

# Returns
- `AbstractRegister`: the register on which we have initialized the bell state
"""
function bell_state(which::Int)

    @assert 1 <= which <= 4 "Input $(which) is not in range 1..4"
    reg = zero_state(2)
    apply!(reg, chain(put(2, 1 => H), cnot(2, 1, 2)))
    (which == 2 || which == 4) && apply!(reg, put(2, 1 => Z))
    (which == 3 || which == 4) && apply!(reg, put(2, 1 => X))
    return reg
end

# Alice perform Measurement using Operator specified by row idx
alice_measure!(reg::AbstractRegister, row::Int) = [real(measure!(op, reg, (1, 3))) for op in op_mtx[row, :]]
# Bob perform Measurement using Operator specified by row idx
bob_measure!(reg::AbstractRegister, col::Int) = [real(measure!(op, reg, (2, 4))) for op in op_mtx[:, col]]

"""
    input_idx(which::String) -> Int

# Arguments
- `which::String`: whether user is asked to provide the index for the row or column

# Returns
- `Int`: index determined by the user
"""
function input_idx(which::String)
    idx = 0
    while idx > 3 || idx < 1
        println("Please provide a $which number in the range [1,3]")
        try
            idx = parse(Int64, readline())
        catch e
            showerror(stdout, e)
            println("Please enter a NUMBER!")
        end
    end
    return idx
end

function main()
    println("Dealer: Let's Play Mermin's Magic Square Game!")

    # ask the user to input row and column index
    i = input_idx("row")
    j = input_idx("column")

    println("Alice will provide three numbers to be written on row: $(i)")
    println("Bob will provide three numbers to be written on column: $(j)")

    # Initialize two pairs of bell state ``(|00> + |11>) ⊗ (|00> + |11>)/2``.
    # First qubit of each bell state pair belongs to Alice and the others belong to Bob.
    reg1 = join(bell_state(1), bell_state(1))

    # do measurement and get answer string
    a1, a2, a3 = alice_measure!(reg1, i)
    b1, b2, b3 = bob_measure!(reg1, j)

    println("Alice should give answer: $a1 , $a2 , $a3")
    println("Bob should give answer: $b1 , $b2 , $b3")
    #We should see j_th element  in [a1,a2,a3] is equal to i_th element in [b1,b2,b3]
    @assert [a1, a2, a3][j] == [b1, b2, b3][i] ("Opps, something is wrong, the j_th item in
       Alice's answer does not match the i_th item in Bob's answer")

    # To reveal some magic behind the scenes, you may check that
    # measurements on each row and each columns commute with
    # each other and all are hermitian
    @assert (all(k -> iscommute(op_mtx[k, :]...), 1:3) &&
             all(k -> iscommute(op_mtx[:, k]...), 1:3) &&
             all(ishermitian, op_mtx))
end


main()
