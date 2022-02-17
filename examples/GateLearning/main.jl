using Yao
using Yao.EasyBuild: general_U4
using Random
using Optim: LBFGS, optimize
using Optim

"""
    learn_u4(u::AbstractMatrix; niter=100)

Learn a general U4 gate. The optimizer is LBFGS.
"""
function learn_u4(u::AbstractBlock; niter=100)
    ansatz = dispatch(general_U4(), :random)
    params = parameters(ansatz)
    println("initial fidelity = $(operator_fidelity(u,ansatz))")
    res = optimize(x->-operator_fidelity(u, dispatch(ansatz, x)),
            (G, x) -> (G .= -operator_fidelity'(u, dispatch(ansatz, x))[2]),
            parameters(ansatz),
            LBFGS(),
            Optim.Options(iterations=niter))
    println("final fidelity = $(operator_fidelity(u,dispatch(ansatz, res.minimizer)))")
    return dispatch!(ansatz, res.minimizer)
end

using Random
Random.seed!(2)
u = matblock(rand_unitary(4))
c = learn_u4(u; niter=150)
