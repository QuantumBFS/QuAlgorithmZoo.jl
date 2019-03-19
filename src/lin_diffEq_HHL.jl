export array_qudiff, prepare_init_state, solve_qudiff, AlgQuDiff, bval, aval
export QuEuler, QuLeapfrog, QuAB2, QuAB3, QuAB4

"""
    Based on : arxiv.org/abs/1010.2745v2

    * array_qudiff(N_t,N,h,A) - generates matrix for k-step solver
    * prepare_init_state(b,x,h,N_t) - generates inital states
    * solve_qudiff - solver

    x' = Ax + b

    * A - input matrix.
    * b - input vector.
    * x - inital vector
    * N - dimension of b (as a power of 2).
    * h - step size.
    * tspan - time span.
"""

"""
    AlgQuDiff
    * step - step for multistep method
    * α - coefficients for xₙ
    * β - coefficent for xₙ'
"""

struct AlgQuDiff
    step::Int
    α::Array
    β::Array
    AlgQuDiff(step,α,β) = new(step,α,β)
end

function bval(alg::AlgQuDiff,t,h,g::Function)
    b = zero(g(1))
    for i in 1:(alg.step)
        b += alg.β[i]*g(t-(i-1)*h)
    end
    b
end

function aval(alg::AlgQuDiff,t,h,g::Function)
    sz, = size(g(1))
    A = Array{ComplexF64}(undef,sz,(alg.step + 1)*sz)
    i_mat = Matrix{Float64}(I, size(g(1)))
    A[1:sz,sz*(alg.step) + 1:sz*(alg.step + 1)] = i_mat
    for i in 1:alg.step
        A[1:sz,sz*(i - 1) + 1: sz*i] = -1*(alg.α[alg.step - i + 1]*i_mat + h*alg.β[alg.step - i + 1]*g(t - (alg.step - i)*h))
    end
    A
end
"""
    Explicit Linear Multistep Methods
"""
QuEuler = AlgQuDiff(1,[1],[1]);
QuLeapfrog = AlgQuDiff(2,[0 1],[2 0])
QuAB2 = AlgQuDiff(2,[1 0], [1.5 -0.5])
QuAB3 = AlgQuDiff(3,[1 0 0], [23/12 -16/12 5/12])
QuAB4 = AlgQuDiff(4,[1 0 0 0], [55/24 -59/24 37/24 -9/24])


function prepare_init_state(tspan::Tuple,x::Vector,h::Float64,g::Function,alg::AlgQuDiff)
    N_t = Int(round(2*(tspan[2] - tspan[1])/h + 3))
    sz, = size(g(1))
    init_state = zeros(ComplexF64,2*(N_t + 1)*sz)
    #inital value
    init_state[1:sz] = x
    for i in 2:(N_t+1)/2 - 1
        b = bval(alg,h*(i - 1) + tspan[1],h,g)
        init_state[Int(sz*(i - 1) + 1):Int(sz*(i))] = h*b
    end
    init_state
end

function array_qudiff(tspan::Tuple,h::Float64,g::Function,alg::AlgQuDiff)
    sz, = size(g(1))
    i_mat = Matrix{Float64}(I, size(g(1)))
    N_t = Int(round(2*(tspan[2] - tspan[1])/h + 3))
    A_ = zeros(ComplexF64,(N_t + 1)*sz,(N_t + 1)*sz)
    # Generates First two rows
    A_[1:sz, 1:sz] = i_mat
    A_[sz + 1:2*sz,1:sz*2] = [-1*(i_mat + h*g(tspan[1])) i_mat]
    #Generates additional rows based on k - step
    for i in 3:alg.step
        A_[Int(sz*(i - 1) + 1):Int(sz*(i)), Int(sz*(i - 3) + 1):Int(sz*i)] = aval(QuAB2,(i-2)*h + tspan[1],h,g)
    end
    for i in alg.step + 1:(N_t + 1)/2 - 1
        A_[Int(sz*(i - 1) + 1):Int(sz*(i)),Int(sz*(i - alg.step - 1) + 1):Int(sz*i)]= aval(alg,(i - 2)*h + tspan[1],h,g)
    end
    #Generates half mirroring matrix
    for i in (N_t + 1)/2:N_t + 1
        A_[Int(sz*(i - 1) + 1):Int(sz*(i)),Int(sz*(i - 2) + 1):Int(sz*i)] = [ -1*i_mat i_mat]
    end
    A_ = [zero(A_) A_;A_' zero(A_)]
    A_
end

function solve_qudiff(alg::AlgQuDiff ,A::Function, b::Function, x::Vector,tspan::Tuple, h::Float64, n_reg::Int)
    mat = array_qudiff(tspan, h, A, alg)
    state = prepare_init_state(tspan, x, h, b, alg)
    λ = maximum(eigvals(mat))
    C_value = minimum(eigvals(mat) .|> abs)*0.01;
    mat = 1/(λ*2)*mat
    state = state*1/(2*λ) |> normalize!
    res = hhlsolve(mat,state, n_reg, C_value)
    res = res/λ
    res
end;
