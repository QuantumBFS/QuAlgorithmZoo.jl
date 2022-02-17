# Variational Quantum Eigensolver

## Setup your python environment
1. set your Python environment in Julia
```
julia> ENV["PYTHON"] = "... path of the python executable ..."
```
5. install [PyCall](https://github.com/JuliaPy/PyCall.jl)
```
pkg> add PyCall
pkg> build PyCall
```
6. install [OpenFermion](https://github.com/quantumlib/OpenFermion)
```bash
pip install openfermion
```
7. run H2_OpenFermion example 
```bash
julia examples/VQE/H2_OpenFermion.jl
```
