# Mermin-Peres Magic Square Game Algorithm
# https://www.gregegan.net/SCIENCE/MPMS/MPMS.html
using Yao
using Yao.EasyBuild, YaoPlots


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
    (which == 3 || which == 4) && apply!(reg, put(2, 2 => X))
    return reg
end

"""
    init_double_bell_state(reg::ArrayReg)

Initialize two pairs of bell state ``(|00> + |11>) ⊗ (|00> + |11>)/2``.

First qubit of each bell state pair belongs to Alice and
the others belong to Bob.
"""
function init_double_bell_state()
    return join(bell_state(1), bell_state(1))
end


"""
    meas_by_idx(reg::ArrayReg, i::Int, j::Int)

Perform Measurement based on index assigned by dealer.

# Arguments
- `reg`: is the quantum register to measure
- 'i::Int': the row idx assigned
- 'j::Int': the column idx assigned

"""
function meas_by_idx(reg::ArrayReg, i::Int, j::Int)

    alice_measure!(reg::AbstractRegister, row::Int) = [real(measure!(op, reg, (1, 3))) for op in op_mtx[row, :]]
    bob_measure!(reg::AbstractRegister, col::Int) = [real(measure!(op, reg, (2, 4))) for op in op_mtx[:, col]]
    ans = zeros(Float64, 6)

    ans[1:3] = alice_measure!(reg, i)'
    ans[4:6] = bob_measure!(reg, j)'
    return ans
end

"""

        end
    end

end

function main()
    # Dealer randomly generate a pair of numbers (i,j)
    # i is the row number given to Alice
    # j is the column number given to Bob
    println("Dealer: Let's Play Mermin's Magic Square Game!")

    println("Please enter a number in the range [1,3] as row index")
    i = parse(Int64,readline())
    println("Please enter a number in the range [1,3] as col index")
    j = parse(Int64,readline())

    # test to make sure in bound
    while i > 3 || i < 1
        println("Please be nice and enter a number in the range [1,3]")
        i = parse(Int64,readline())
    end

    while j > 3 || j < 1
        println("Please be nice and enter a number in the range [1,3]")
        j = parse(Int64,readline())
    end

    println("Alice will provide three numbers to be written on row: $(i)")
    println("Bob will provide three numbers to be written on column: $(j)")

    # initialize quantum registero
    reg1 = zero_state(4)

    init_double_bell_state(reg1)

    # do measurement and get answer string
    # I sort of cheated here, since we know Alice Bob
    # need to adhere to the rules, the third measurement
    # won't be necessary hence I just used the fact
    # that a3 = 1 * a1 * a2 and b3 = -1 * b1 * b2
    a1, a2, b1, b2 = meas_by_idx(reg1, i, j)

    println("Alice should give answer: $a1 , $a2 , $(1*a1*a2)")
    println("Bob should give answer: $b1 , $b2 , $(-1*b1*b2)")
    println("We should see $([a1,a2,1*a1*a2][j]) is equal to $([b1,b2,-1*b1*b2][i])")
    @assert [a1, a2, 1 * a1 * a2][j] == [b1, b2, -1 * b1 * b2][i]

    println("To reveal some magic behind the scenes, you may check that measurements on each row and each columns commute with each other and all are hermitian")
    @assert all(k -> iscommute(op_mtx[k, :]...), 1:3) && all(k -> iscommute(op_mtx[:, k]...), 1:3) && all(ishermitian, op_mtx)
end

main()
