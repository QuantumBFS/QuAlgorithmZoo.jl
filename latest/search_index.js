var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": "CurrentModule = QuAlgorithmZoo"
},

{
    "location": "#Quantum-Algorithm-Zoo-1",
    "page": "Home",
    "title": "Quantum Algorithm Zoo",
    "category": "section",
    "text": "A curated implementation of quantum algorithms with Yao.jl"
},

{
    "location": "tutorial/QCBM/#",
    "page": "Quantum Circuit Born Machine",
    "title": "Quantum Circuit Born Machine",
    "category": "page",
    "text": "EditURL = \"https://github.com/QuantumBFS/QuAlgorithmZoo.jl/blob/master/examples/QCBM.jl\""
},

{
    "location": "tutorial/QCBM/#Quantum-Circuit-Born-Machine-1",
    "page": "Quantum Circuit Born Machine",
    "title": "Quantum Circuit Born Machine",
    "category": "section",
    "text": "Quantum circuit born machine is a fresh approach to quantum machine learning. It use a parameterized quantum circuit to learning machine learning tasks with gradient based optimization. In this tutorial, we will show how to implement it with Yao (幺) framework.using Yao, Plots@doc 幺"
},

{
    "location": "tutorial/QCBM/#Training-target-1",
    "page": "Quantum Circuit Born Machine",
    "title": "Training target",
    "category": "section",
    "text": "a gaussian distributionfunction gaussian_pdf(n, μ, σ)\n    x = collect(1:1<<n)\n    pl = @. 1 / sqrt(2pi * σ^2) * exp(-(x - μ)^2 / (2 * σ^2))\n    pl / sum(pl)\nendf(x left mu sigma^2right) = frac1sqrt2pisigma^2 e^-frac(x-mu)^22sigma^2const n = 6\nconst maxiter = 20\npg = gaussian_pdf(n, 2^5-0.5, 2^4)\nfig = plot(0:1<<n-1, pg)"
},

{
    "location": "tutorial/QCBM/#Build-Circuits-1",
    "page": "Quantum Circuit Born Machine",
    "title": "Build Circuits",
    "category": "section",
    "text": ""
},

{
    "location": "tutorial/QCBM/#Building-Blocks-1",
    "page": "Quantum Circuit Born Machine",
    "title": "Building Blocks",
    "category": "section",
    "text": "Gates are grouped to become a layer in a circuit, this layer can be Arbitrary Rotation or CNOT entangler. Which are used as our basic building blocks of Born Machines."
},

{
    "location": "tutorial/QCBM/#Arbitrary-Rotation-1",
    "page": "Quantum Circuit Born Machine",
    "title": "Arbitrary Rotation",
    "category": "section",
    "text": "Arbitrary Rotation is built with Rotation Gate on Z, Rotation Gate on X and Rotation Gate on Z:Rz(theta) cdot Rx(theta) cdot Rz(theta)Since our input will be a 0dots 0rangle state. The first layer of arbitrary rotation can just use Rx(theta) cdot Rz(theta) and the last layer of arbitrary rotation could just use Rz(theta)cdot Rx(theta)In 幺, every Hilbert operator is a block type, this includes all quantum gates and quantum oracles. In general, operators appears in a quantum circuit can be divided into Composite Blocks and Primitive Blocks.We follow the low abstraction principle and thus each block represents a certain approach of calculation. The simplest Composite Block is a Chain Block, which chains other blocks (oracles) with the same number of qubits together. It is just a simple mathematical composition of operators with same size. e.g.textchain(X Y Z) iff X cdot Y cdot ZWe can construct an arbitrary rotation block by chain Rz, Rx, Rz together.chain(Rz(0), Rx(0), Rz(0))chain(X)Rx, Ry and Rz will construct new rotation gate, which are just shorthands for rot(X, 0.0), etc.Then, let\'s pile them up vertically with another method called rollrepeatlayer(x::Symbol) = layer(Val(x))\nlayer(::Val{:first}) = rollrepeat(chain(Rx(0), Rz(0)))\nrollrepeat(chain(Rx(0), Rz(0)))(4)In 幺, the factory method rollrepeat will construct a block called Roller. It is mathematically equivalent to the kronecker product of all operators in this layer:rollrepeat(n U) iff roll(n texti=U for i = 1n) iff kron(n texti=U for i=1n) iff U otimes dots otimes U@doc Valroll(4, i=>X for i = 1:4)rollrepeat(4, X)kron(4, i=>X for i = 1:4)However, kron is calculated differently comparing to roll. In principal, Roller will be able to calculate small blocks with same size with higher efficiency. But for large blocks Roller may be slower. In 幺, we offer you this freedom to choose the most suitable solution.all factory methods will lazy evaluate the first arguements, which is the number of qubits. It will return a lambda function that requires a single interger input. The instance of desired block will only be constructed until all the information is filled.rollrepeat(X)rollrepeat(X)(4)When you filled all the information in somewhere of the declaration, 幺 will be able to infer the others.chain(4, rollrepeat(X), rollrepeat(Y))We will now define the rest of rotation layerslayer(::Val{:last}) = rollrepeat(chain(Rz(0), Rx(0)))\nlayer(::Val{:mid}) = rollrepeat(chain(Rz(0), Rx(0), Rz(0)))"
},

{
    "location": "tutorial/QCBM/#CNOT-Entangler-1",
    "page": "Quantum Circuit Born Machine",
    "title": "CNOT Entangler",
    "category": "section",
    "text": "Another component of quantum circuit born machine is several CNOT operators applied on different qubits.entangler(pairs) = chain(control([ctrl, ], target=>X) for (ctrl, target) in pairs)We can then define such a born machinefunction QCBM(n, nlayer, pairs)\n    circuit = chain(n)\n    push!(circuit, layer(:first))\n\n    for i = 1:(nlayer - 1)\n        push!(circuit, cache(entangler(pairs)))\n        push!(circuit, layer(:mid))\n    end\n\n    push!(circuit, cache(entangler(pairs)))\n    push!(circuit, layer(:last))\n\n    circuit\nendWe use the method cache here to tag the entangler block that it should be cached after its first run, because it is actually a constant oracle. Let\'s see what will be constructedQCBM(4, 1, [1=>2, 2=>3, 3=>4, 4=>1])"
},

{
    "location": "tutorial/QCBM/#MMD-Loss-and-Gradients-1",
    "page": "Quantum Circuit Born Machine",
    "title": "MMD Loss & Gradients",
    "category": "section",
    "text": "The MMD loss is describe below:beginaligned\nmathcalL = left sum_x p theta(x) phi(x) - sum_x pi(x) phi(x) right^2\n            = langle K(x y) rangle_x sim p_theta ysim p_theta - 2 langle K(x y)\n                  rangle_xsim p_theta ysim pi + langle K(x y) rangle_xsimpi ysimpi\nendalignedWe will use a squared exponential kernel here.struct Kernel\n    sigma::Float64\n    matrix::Matrix{Float64}\nend\n\nfunction Kernel(nqubits, sigma)\n    basis = collect(0:(1<<nqubits - 1))\n    Kernel(sigma, kernel_matrix(basis, basis, sigma))\nend\n\nexpect(kernel::Kernel, px::Vector{Float64}, py::Vector{Float64}) = px\' * kernel.matrix * py\nloss(qcbm, kernel::Kernel, ptrain) = (p = get_prob(qcbm) - ptrain; expect(kernel, p, p))Next, let\'s define the kernel matrixfunction kernel_matrix(x, y, sigma)\n    dx2 = (x .- y\').^2\n    gamma = 1.0 / (2 * sigma)\n    K = exp.(-gamma * dx2)\n    K\nend"
},

{
    "location": "tutorial/QCBM/#Gradients-1",
    "page": "Quantum Circuit Born Machine",
    "title": "Gradients",
    "category": "section",
    "text": "the gradient of MMD loss isbeginaligned\nfracpartial mathcalLpartial theta^i_l = langle K(x y) rangle_xsim p_theta^+ ysim p_theta - langle K(x y) rangle_xsim p_theta^- ysim p_theta\n- langle K(x y) rangle _xsim p_theta^+ ysimpi + langle K(x y) rangle_xsim p_theta^- ysimpi\nendalignedWe have to update one parameter of each rotation gate each time, and calculate its gradient then collect them. Since we will need to calculate the probability from the state vector frequently, let\'s define a shorthand first.Firstly, you have to define a quantum register. Each run of a QCBM\'s input is a simple 00cdots 0rangle state. We provide string literal bit to help you define one-hot state vectors like thisr = register(bit\"0000\")circuit = QCBM(6, 10, [1=>2, 3=>4, 5=>6, 2=>3, 4=>5, 6=>1]);Now, we define its shorthandget_prob(qcbm) = apply!(register(bit\"0\"^nqubits(qcbm)), qcbm) |> statevec .|> abs2We will first iterate through each layer contains rotation gates and allocate an array to store our gradientfunction gradient(n, nlayers, qcbm, kernel, ptrain)\n    prob = get_prob(qcbm)\n    grad = zeros(real(datatype(qcbm)), nparameters(qcbm))\n    idx = 0\n    for ilayer = 1:2:(2 * nlayers + 1)\n        idx = grad_layer!(grad, idx, prob, qcbm, qcbm[ilayer], kernel, ptrain)\n    end\n    grad\nendThen we iterate through each rotation gate.function grad_layer!(grad, idx, prob, qcbm, layer, kernel, ptrain)\n    count = idx\n    for each_line in blocks(layer)\n        for each in blocks(each_line)\n            gradient!(grad, count+1, prob, qcbm, each, kernel, ptrain)\n            count += 1\n        end\n    end\n    count\nendWe update each parameter by rotate it -pi2 and pi2function gradient!(grad, idx, prob, qcbm, gate, kernel, ptrain)\n    dispatch!(+, gate, pi / 2)\n    prob_pos = get_prob(qcbm)\n\n    dispatch!(-, gate, pi)\n    prob_neg = get_prob(qcbm)\n\n    dispatch!(+, gate, pi / 2) # set back\n\n    grad_pos = expect(kernel, prob, prob_pos) - expect(kernel, prob, prob_neg)\n    grad_neg = expect(kernel, ptrain, prob_pos) - expect(kernel, ptrain, prob_neg)\n    grad[idx] = grad_pos - grad_neg\n    grad\nend"
},

{
    "location": "tutorial/QCBM/#Optimizer-1",
    "page": "Quantum Circuit Born Machine",
    "title": "Optimizer",
    "category": "section",
    "text": "We will use the Adam optimizer. Since we don\'t want you to install another package for this, the following code for this optimizer is copied from Knet.jlReference: Kingma, D. P., & Ba, J. L. (2015). Adam: a Method for Stochastic Optimization. International Conference on Learning Representations, 1–13.using LinearAlgebra\n\nmutable struct Adam\n    lr::AbstractFloat\n    gclip::AbstractFloat\n    beta1::AbstractFloat\n    beta2::AbstractFloat\n    eps::AbstractFloat\n    t::Int\n    fstm\n    scndm\nend\n\nAdam(; lr=0.001, gclip=0, beta1=0.9, beta2=0.999, eps=1e-8)=Adam(lr, gclip, beta1, beta2, eps, 0, nothing, nothing)\n\nfunction update!(w, g, p::Adam)\n    gclip!(g, p.gclip)\n    if p.fstm===nothing; p.fstm=zero(w); p.scndm=zero(w); end\n    p.t += 1\n    lmul!(p.beta1, p.fstm)\n    BLAS.axpy!(1-p.beta1, g, p.fstm)\n    lmul!(p.beta2, p.scndm)\n    BLAS.axpy!(1-p.beta2, g .* g, p.scndm)\n    fstm_corrected = p.fstm / (1 - p.beta1 ^ p.t)\n    scndm_corrected = p.scndm / (1 - p.beta2 ^ p.t)\n    BLAS.axpy!(-p.lr, @.(fstm_corrected / (sqrt(scndm_corrected) + p.eps)), w)\nend\n\nfunction gclip!(g, gclip)\n    if gclip == 0\n        g\n    else\n        gnorm = vecnorm(g)\n        if gnorm <= gclip\n            g\n        else\n            BLAS.scale!(gclip/gnorm, g)\n        end\n    end\nendThe training of the quantum circuit is simple, just iterate through the steps.function train!(qcbm, ptrain, optim; learning_rate=0.1, niter=50)initialize the parameters    params = 2pi * rand(nparameters(qcbm))\n    dispatch!(qcbm, params)\n    kernel = Kernel(nqubits(qcbm), 0.25)\n\n    n, nlayers = nqubits(qcbm), (length(qcbm)-1)÷2\n    history = Float64[]\n\n    for i = 1:niter\n        grad = gradient(n, nlayers, qcbm, kernel, ptrain)\n        curr_loss = loss(qcbm, kernel, ptrain)\n        push!(history, curr_loss)\n        params = collect(parameters(qcbm))\n        update!(params, grad, optim)\n        dispatch!(qcbm, params)\n    end\n    history\nendoptim = Adam(lr=0.1)\nhis = train!(circuit, pg, optim, niter=50, learning_rate=0.1)\nplot(1:50, his, xlabel=\"iteration\", ylabel=\"loss\")p = get_prob(circuit)\np = plot(0:1<<n-1, p, xlabel=\"x\", ylabel=\"p\")\nplot!(p, 0:1<<n-1, pg)This page was generated using Literate.jl."
},

]}
