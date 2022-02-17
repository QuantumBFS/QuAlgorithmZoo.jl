using Yao
using Test, Random
using Optim: LBFGS, optimize
using Optim

# port the `Matrix` function to Yao's AD.
using Zygote

ansatz = general_U4() * put(2, 1=>phase(0.0))  # initial values are 0, here, we attach a global phase.

u = rand_unitary(4)

function loss(params)
    m = Matrix(dispatch(ansatz, params))
    sum(abs.(u .- m))
end

"""
    learn_u4(u::AbstractMatrix; niter=100)

Learn a general U4 gate. The optimizer is LBFGS.
"""
function learn_u4(u::AbstractMatrix; niter=100)
    params = randn(nparameters(ansatz))
    g!(G, x) = (G .= Zygote.gradient(loss, x)[1])
    res = optimize(loss, g!, params,
                    LBFGS(), Optim.Options(iterations=niter))
    println("final loss = $(loss(res.minimizer))")
    return dispatch(ansatz, res.minimizer)
end

using Random
Random.seed!(3)
c = learn_u4(u; niter=150)
