# QuAlgorithmZoo

A curated implementation of quantum algorithms with [Yao.jl](https://github.com/QuantumBFS/Yao.jl)

## Installation

`QuAlgorithmZoo.jl` is not registered, please use the following command to download the examples:

```bash
$ git clone https://github.com/QuantumBFS/QuAlgorithmZoo.jl.git
```

Then open directory `.julia/dev/QuAlgorithmZoo/examples` to find algorithms.

## Contents

- QFT (`Yao.EasyBuild.qft_circuit`)
- Phase Estimation (`Yao.EasyBuild.phase_estimation_circuit`)
- Hadamard Test (`Yao.EasyBuild.hadamard_test_circuit`)
- State Overlap Algorithms (`Yao.EasyBuild.swap_test_circuit`)

In examples folder, you will find

- [Variational Quantum Eigensolver](examples/VQE)
- [Variational Quantum Eigensolver (with OpenFermion)](examples/VQE_openfermion)
- [Quantum SVD](examples/QSVD)
- [Imaginary Time Evolution Quantum Eigensolver](examples/GroundStateSolvers)
- [HHL algorithm](examples/HHL)
- [QAOA](examples/QAOA)
- [Quantum Circuit Born Machine](examples/QCBM)
- [Quantum GAN](examples/QuGAN)
- [Shor](examples/Shor)
- [Grover search](examples/Grover)
- [Learning a general two qubit unitary gate](examples/GateLearning)

Examples of using Yao in other projects
- [QuODE](https://github.com/QuantumBFS/QuDiffEq.jl)
- [TensorNetwork Inspired Circuits](https://github.com/GiggleLiu/QuantumPEPS.jl)
- [Beta-VQE](https://github.com/wangleiphy/BetaVQE.jl)
- [Quantum Neural Network Classifier](https://github.com/LWKJJONAK/Quantum_Neural_Network_Classifiers)

## License

QuAlgorithmZoo.jl is released under Apache License 2.0.
