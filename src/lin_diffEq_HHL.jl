export Array_QuDiff, prepare_init_state, solve_QuDiff, alg_lin_diffeq, bval, Aval
export QuEuler, QuLeapfrog, QuAB2, QuAB3, QuAB4

"""
    Based on : arxiv.org/abs/1010.2745v2

    * Array_QuDiff(N_t,N,h,A) - generates matrix for k-step solver
    * prepare_init_state(b,x,h,N_t) - generates inital states
    * solve_QuDiff - solver

    x' = Ax + b

    * A - input matrix.
    * b - input vector.
    * x - inital vector
    * N - dimension of b (as a power of 2).
    * h - step size.
    * tspan - time span.
"""

"""
    alg_lin_diffeq
    * step - step for multistep method
    * α - coefficients for xₙ
    * β - coefficent for xₙ'
"""
struct alg_lin_diffeq
  step::Int
  α::Array
  β::Array
  alg_lin_diffeq(step,α,β) = new(step,α,β)
end

function bval(alg::alg_lin_diffeq,t,h,g::Function)
  b = zero(g(1))
  for i = 1:(alg.step)
    b += alg.β[i]*g(t-(i-1)*h)
  end
  b
end

function Aval(alg::alg_lin_diffeq,t,h,g::Function)
  I_mat = Matrix{Float64}(I, size(g(1)))
  A = Matrix{Float64}(I, size(g(1)))
  for i = 1:(alg.step)
    A = [-1*(alg.α[i]*I_mat + h*alg.β[i]*g(t-(i-1)*h)) A]
  end
  A
end
"""
    Explicit Linear Multistep Methods
"""
QuEuler = alg_lin_diffeq(1,[1],[1]);
QuLeapfrog = alg_lin_diffeq(2,[0 1],[2 0])
QuAB2 = alg_lin_diffeq(2,[1 0], [1.5 -0.5])
QuAB3 = alg_lin_diffeq(3,[1 0 0], [23/12 -16/12 5/12])
QuAB4 = alg_lin_diffeq(4,[1 0 0 0], [55/24 -59/24 37/24 -9/24])


function prepare_init_state(tspan::Tuple,x::Vector,h::Float64,g::Function,alg::alg_lin_diffeq)
  init_state = x;
  N_t = Int(round(2*(tspan[2] - tspan[1])/h + 3))
  b = similar(g(1))
  for i = 1:N_t
   if i < (N_t+1)/2 -1
      b = bval(alg,h*(i-1) + tspan[1],h,g)
      init_state = [init_state;h*b]
    else
      init_state = [init_state;zero(b)]
   end

  end
  init_state = [init_state;zero(init_state)]
  init_state
end

function Array_QuDiff(tspan::Tuple,h::Float64,g::Function,alg::alg_lin_diffeq)
  N_t = Int(round(2*(tspan[2] - tspan[1])/h + 3))
  sz = size(g(1))
  I_mat = Matrix{Float64}(I, size(g(1)));
  A_ = I_mat;
  # Generates First two rows
  A_ = [A_ zeros(sz[1],sz[1]*(N_t));-1*(I_mat + h*g(tspan[1])) I_mat zeros(sz[1],sz[1]*(N_t-1))]
  #Generates additional rows based on k - step
  if alg.step > 1
    A_ = [A_; Aval(QuAB2,2*h + tspan[1],h,g) zeros(sz[1],sz[1]*(N_t-2))]
  end

  for i = 3:alg.step-1
    tA_ = Aval(QuAB2,(i-1)*h + tspan[1],h,g)
    tA_ = [zeros(sz[1],sz[1]*(i-2)) tA_]
    A_ = [A_; tA_ zeros(sz[1],sz[1]*(N_t-i))]
  end

  if alg.step > 2
    A_ = [A_;Aval(alg,(alg.step)*h + tspan[1],h,g) zeros(sz[1],sz[1]*(N_t-alg.step)) ]
  end

  for i = alg.step+1: (N_t + 1)/2 -2
    tA_ = Aval(alg,(i-1)*h + tspan[1],h,g)
    A_= [A_;zeros(sz[1],Int(sz[1]*(i-alg.step))) tA_ zeros(sz[1],Int(sz[1]*(N_t-i)))]
   end
  #Generates half mirroring matrix
  for i = (N_t + 1)/2  : N_t
    A_ = [A_; zeros(sz[1],Int(sz[1]*(i-2))) -1*I_mat I_mat zeros(sz[1],Int(sz[1]*(N_t + 1 -i)))]
  end
  A_ = [A_; zeros(sz[1],sz[1]*(N_t-1)) -1*I_mat I_mat];

  A_ = [zero(A_) A_;A_' zero(A_)]
  A_
end

function solve_QuDiff(alg::alg_lin_diffeq ,A::Function, b::Function, x::Vector,tspan::Tuple, h::Float64)

  mat = Array_QuDiff(tspan,h,A,alg)
  state = prepare_init_state(tspan,x,h,b,alg)
  λ = maximum(eigvals(mat))
  C_value = minimum(eigvals(mat) .|> abs)*0.01;
  n_reg = 12;
  mat = 1/(λ*2)*mat
  state = state*1/(2*λ) |> normalize!
  res = hhlsolve(mat,state, n_reg, C_value)
  res = res/λ
  res
end
