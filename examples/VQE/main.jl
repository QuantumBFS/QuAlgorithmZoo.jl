# # [Variational Quantum Eigensolver](@id VQE)

# ## References
# * [A variational eigenvalue solver on a quantum processor](https://arxiv.org/abs/1304.3061)
# * [Variational Quantum Eigensolver with Fewer Qubits](https://arxiv.org/abs/1902.02663)

# ## Define a hamiltonian

# construct a 5-site heisenberg hamiltonian

using Yao, Yao.EasyBuild
using KrylovKit: eigsolve
import Optimisers

N = 5
hami = heisenberg(N)

# The ground state can be obtained by a sparse matrix ground state solver.
# The high performance `mat` function in `Yao.jl` makes computation time lower than `10s`
# to construct a `20` site Heisenberg hamltonian

function ed_groundstate(h::AbstractBlock)
    E, V = eigsolve(h |> mat, 1, :SR, ishermitian=true)
    E[1], V[1]
end

ed_groundstate(hami)

# Here we use the `heisenberg` hamiltonian that defined in `Yao.EasyBuild`,
# for tutorial purpose, we pasted the code for construction here.
# ```julia
# function heisenberg(nbit::Int; periodic::Bool=true)
#    sx = i->put(nbit, i=>X)
#    sy = i->put(nbit, i=>Y)
#    sz = i->put(nbit, i=>Z)
#    map(1:(periodic ? nbit : nbit-1)) do i
#        j=i%nbit+1
#        sx(i)*sx(j)+sy(i)*sy(j)+sz(i)*sz(j)
#     end |> sum
# end
# ```

# ## Define an ansatz
# As an ansatz, we use the canonical circuit for demonstration [`variational_circuit`](@ref Yao.EasyBuild.variational_circuit)
# defined in [`Yao.EasyBuild.jl`](https://github.com/QuantumBFS/Yao.jl).
c = variational_circuit(N)
dispatch!(c, :random)

# ## Run
# Use the `ADAM` optimizer for parameter optimization,
# we provide a poorman's implementation in `QuAlgorithmZoo`
params = parameters(c)
optimizer = Optimisers.setup(Optimisers.ADAM(0.01), params)
niter = 100

for i = 1:niter
    ## `expect'` gives the gradient of an observable.
    grad_input, grad_params = expect'(hami, zero_state(N) => c)

    ## feed the gradients into the circuit.
    dispatch!(c, Optimisers.update!(optimizer, params, grad_params))
    println("Step $i, Energy = $(expect(hami, zero_state(N) |> c))")
end

# ## Hydrogen atoms
# The hamiltonian can be found in arXiv: 1704.05018, table S2
function hydrogen_hamiltonian()
    Z1 = put(2,1=>Z)
    Z2 = put(2,2=>Z)
    X1 = put(2,1=>X)
    X2 = put(2,2=>X)
    0.011280*Z1*Z2 + 0.397936*Z1 + 0.397936*Z2 + 0.180931*X1*X2
end

hami = hydrogen_hamiltonian()
emin = eigvals(Matrix(mat(h)))[1]

c = variational_circuit(2)
dispatch!(c, :random)


params = parameters(c)
optimizer = Optimsers.setup(Optimisers.ADAM(0.01), params)
niter = 100

for i = 1:niter
    ## `expect'` gives the gradient of an observable.
    grad_input, grad_params = expect'(hami, zero_state(2) => c)

    ## feed the gradients into the circuit.
    dispatch!(c, Optimisers.update!(optimizer, params, grad_params))
    println("Step $i, Energy = $(expect(hami, zero_state(2) |> c))")
end

