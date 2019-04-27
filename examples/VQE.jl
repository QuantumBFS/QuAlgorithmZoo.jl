using Yao, Yao.Blocks
using QuAlgorithmZoo
using KrylovKit

function ed_groundstate(h::MatrixBlock)
    E, V = eigsolve(h |> mat, 1, :SR, ishermitian=true)
    println("Ground State Energy is $(E[1])")
    register(V[1])
end

N = 5
c = random_diff_circuit(N, N, [i=>mod(i,N)+1 for i=1:N], mode=:Merged) |> autodiff(:QC)
dispatch!(c, :random)
hami = heisenberg(N)

dbs = collect(AbstractDiff, c)

function scan(db, nparam::Int)
    x = LinRange(0, 2pi, nparam)
    y = zeros(nparam)
    for (i,x) in enumerate(x)
        dispatch!(db, x)
        y[i] = expect(hami, zero_state(N) |> c) |> real
    end
    x, y
end

using Plots
xs, ys = scan(dbs[10], 100)

# see a sine curve
plot(xs, ys)

# get the exact ground state
ed_groundstate(hami)

# vqe ground state
vqe_solve(c, hami)
