# Mermin-Peres Magic Square Game Algorithm
# https://www.gregegan.net/SCIENCE/MPMS/MPMS.html
using Yao
using Yao.EasyBuild, YaoPlots
using Random

"""
    init_double_bell_state(reg::ArrayReg)

Initialize two pairs of bell state.

First qubit of each bell state pair belongs to Alice and
the others belong to Bob.

# Arguments
- 'reg::ArrayReg': the array of quantum register to initize state upon
"""
function init_double_bell_state(reg::ArrayReg)
    for k in [1, 3]
        apply!(reg, put(4, k => H))
        apply!(reg, cnot(4, k, k + 1))
    end
end


"""
    meas_by_idx(reg::ArrayReg, i::Int, j::Int)

Perform Measurement based on index assigned by dealer.

# Arguments
- 'reg::ArrayReg': the array of quantum register to measure
- 'i::Int': the row idx assigned
- 'j::Int': the column idx assigned

"""
function meas_by_idx(reg::ArrayReg, i::Int, j::Int)
    # measurement operator table
    # see Table S1 in https://arxiv.org/abs/2206.12042
    op_mtx = Matrix{KronBlock}(undef, 3, 3)
    op_mtx = [kron(I2, Z) kron(Z, I2) kron(Z, Z); kron(X, I2) kron(I2, X) kron(X, X); kron(-X, Z) kron(-Z, X) kron(Y, Y)]

    ans = zeros(Float64, 4)
    for k in 1:2
        # do measurement for Alice depending on the index given
        @inbounds ans[k] = measure!(op_mtx[i, k], reg, (1, 3))[1]
        # do measurement for Bob depending on the index given
        @inbounds ans[k+2] = measure!(op_mtx[k, j], reg, (2, 4))[1]
    end

    return real.(ans)
end


"""
    Show the measurement in each row and column commutes and are Hermitian.
"""
function show_hermitian_and_commute()
    op_mtx = Matrix{KronBlock}(undef, 3, 3)
    op_mtx = [kron(I2, Z) kron(Z, I2) kron(Z, Z); kron(X, I2) kron(I2, X) kron(X, X); kron(-X, Z) kron(-Z, X) kron(Y, Y)]

    # show all are hermitian
    for j in 1:3
        for i in 1:3
            @inbounds println("$(op_mtx[i,j].blocks) is Hermitian: $(ishermitian(op_mtx[i,j]))")
        end
    end

    # show elements in each row commutes
    for j in 1:3
        for i in 1:2
            @inbounds println("$(op_mtx[i,j].blocks) commutes with $(op_mtx[i+1,j].blocks): $(iscommute(op_mtx[i,j],op_mtx[i+1,j]))")
        end
    end

    # show elements in each column commutes
    for j in 1:2
        for i in 1:3
            @inbounds println("$(op_mtx[i,j].blocks) commutes with $(op_mtx[i,j+1].blocks): $(iscommute(op_mtx[i,j],op_mtx[i,j+1]))")
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
    show_hermitian_and_commute()
end

main()
