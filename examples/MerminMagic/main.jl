# Mermin-Peres Magic Square Game Algorithm
# https://www.gregegan.net/SCIENCE/MPMS/MPMS.html
using Yao
using Yao.EasyBuild, YaoPlots
using Random

# initialize four qubits to be two pairs of Bell Pairs (00)
bell_state_00(ctr_idx::Int) = chain(4, put(ctr_idx => H), control(ctr_idx, (ctr_idx + 1) => X))

function play()
    # measurement operator table
    op_mtx1 = Matrix{ConstantGate{1,2}}(undef, 3, 3)
    op_mtx1 = [I2 Z Z; X I2 X; -X -Z Y]
    op_mtx2 = Matrix{ConstantGate{1,2}}(undef, 3, 3)
    op_mtx2 = [Z I2 Z; I2 X X; Z X Y]

    # Dealer randomly generate a pair of numbers (i,j)
    # i is the row number given to Alice
    # j is the column number given to Bob
    println("Dealer: Let's Play Mermin's Magic Square Game!")
    (i, j) = rand([1, 2, 3], (1, 2))
    println("Alice will provide three numbers to be written on row: $(i)")
    println("Bob will provide three numbers to be written on column: $(j)")

    # initialize quantum registero
    reg1 = zero_state(4)

    # initialize circuit
    circ = chain(bell_state_00(1), bell_state_00(3))

    reg1 |> circ

    # do measurement for Alice depending on the index given
    a1 = measure!(kron(op_mtx1[i, 1], I2, op_mtx2[i, 1], I2), reg1)[1]
    a2 = measure!(kron(op_mtx1[i, 2], I2, op_mtx2[i, 2], I2), reg1)[1]
    # do measurement for Bob depending on the index given
    b1 = measure!(kron(I2, op_mtx1[1, j], I2, op_mtx2[1, j]), reg1)[1]
    b2 = measure!(kron(I2, op_mtx1[2, j], I2, op_mtx2[2, j]), reg1)[1]

    println("Alice should give answer: $a1 , $a2 , $(1*a1*a2)")
    println("Bob should give answer: $b1 , $b2 , $(-1*b1*b2)")
    println("We should see $([a1,a2,1*a1*a2][j]) is equal to $([b1,b2,-1*b1*b2][i])")
    @assert [a1, a2, 1 * a1 * a2][j] == [b1, b2, -1 * b1 * b2][i]
end

function main()
    play()
end

main()
