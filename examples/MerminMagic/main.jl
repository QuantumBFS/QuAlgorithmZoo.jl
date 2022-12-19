# Mermin-Peres Magic Square Game Algorithm
# https://www.gregegan.net/SCIENCE/MPMS/MPMS.html
using Yao
using Yao.EasyBuild, YaoPlots
using Random

# initialize four qubits to be two pairs of Bell Pairs (00)
function init_double_bell_state(reg::ArrayReg)
    for k in [1,3]
        apply!(reg,put(4,k=>H))
        apply!(reg,cnot(4,k,k+1))
    end
end

function meas_by_idx(reg::ArrayReg, i::Int, j::Int)
    # measurement operator table
    op_mtx1 = Matrix{ConstantGate{1,2}}(undef, 3, 3)
    op_mtx1 = [I2 Z Z; X I2 X; -X -Z Y]
    op_mtx2 = Matrix{ConstantGate{1,2}}(undef, 3, 3)
    op_mtx2 = [Z I2 Z; I2 X X; Z X Y]

    ans = zeros(Float64, 4)
    for k in 1:2
        # do measurement for Alice depending on the index given
        ans[k] = measure!(kron(op_mtx1[i, k], I2, op_mtx2[i, k], I2), reg)[1]
        # do measurement for Bob depending on the index given
        ans[k+2] = measure!(kron(I2, op_mtx1[k, j], I2, op_mtx2[k, j]), reg)[1]
    end

    return real.(ans)
end

function play()

    # Dealer randomly generate a pair of numbers (i,j)
    # i is the row number given to Alice
    # j is the column number given to Bob
    println("Dealer: Let's Play Mermin's Magic Square Game!")
    (i, j) = rand([1, 2, 3], (1, 2))
    println("Alice will provide three numbers to be written on row: $(i)")
    println("Bob will provide three numbers to be written on column: $(j)")

    # initialize quantum registero
    reg1 = zero_state(4)

    init_double_bell_state(reg1)

    # do measurement and get answer string
    a1, a2, b1, b2 = meas_by_idx(reg1, i, j)

    println("Alice should give answer: $a1 , $a2 , $(1*a1*a2)")
    println("Bob should give answer: $b1 , $b2 , $(-1*b1*b2)")
    println("We should see $([a1,a2,1*a1*a2][j]) is equal to $([b1,b2,-1*b1*b2][i])")
    @assert [a1, a2, 1 * a1 * a2][j] == [b1, b2, -1 * b1 * b2][i]
end

function main()
    play()
end

# main()
