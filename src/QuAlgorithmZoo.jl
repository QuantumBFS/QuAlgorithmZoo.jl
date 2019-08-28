module QuAlgorithmZoo

using LinearAlgebra
using Yao, BitBasis
using YaoExtensions

include("Grover.jl")
include("Adam.jl")
include("PhaseEstimation.jl")
include("hamiltonian_solvers.jl")
include("HadamardTest.jl")
include("QSVD.jl")

@deprecate random_diff_circuit variational_circuit

end # module
