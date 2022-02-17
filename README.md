# QuAlgorithmZoo

A curated implementation of quantum algorithms with [Yao.jl](https://github.com/QuantumBFS/Yao.jl)

## Installation

QuAlgorithmZoo.jl is not registered, please use the following command:

```julia
pkg> dev https://github.com/QuantumBFS/QuAlgorithmZoo.jl.git
```

Then open directory `.julia/dev/QuAlgorithmZoo/examples` to find algorithms.

## Contents

- QFT (`Yao.EasyBuild.qft_circuit`)
- Phase Estimation (`Yao.EasyBuild.phase_estimation_circuit`)
- Hadamard Test (`Yao.EasyBuild.hadamard_test_circuit`)
- State Overlap Algorithms (`Yao.EasyBuild.swap_test_circuit`)

In examples folder, you will find

- [Variational Quantum Eigensolver](examples/QSVD)
- [Quantum SVD](examples/QSVD)
- [Imaginary Time Evolution Quantum Eigensolver](examples/GroundStateSolvers)
- [HHL algorithm](examples/HHL)
- [QAOA](examples/QAOA)
- [Quantum Circuit Born Machine](examples/QCBM)
- [QuGAN](examples/QuGAN)
- [Shor](examples/Shor)
- [Grover search](examples/Grover)

Examples of using Yao in other projects
- [QuODE](https://github.com/QuantumBFS/QuDiffEq.jl)
- [TensorNetwork Inspired Circuits](https://github.com/GiggleLiu/QuantumPEPS.jl)

## License

QuAlgorithmZoo.jl is released under Apache License 2.0.
