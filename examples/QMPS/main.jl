# # Solving TFI model with only 2 qubits - the 幺 simulation

# Reference: Variational Quantum Eigensolver with Fewer Qubits¶
# Jin-Guo Liu, Yi-Hong Zhang, Yuan Wan, Lei Wang

# https://arxiv.org/abs/1902.02663
# Original post: https://giggleliu.github.io/TwoQubit-VQE.html

using Yao, YaoPlots
using Statistics: mean
using LinearAlgebra
using Plots
using Optimisers

# ## Build a quantum circuit inspired by MPS

rotor(noleading::Bool=false, notrailing::Bool=false) = noleading ? (notrailing ? Rx(0) : chain(Rx(0), Rz(0))) : (notrailing ? chain(Rz(0), Rx(0)) : chain(Rz(0), Rx(0), Rz(0)))

"""
    twoqubit_circuit(nlayer::Int, nrepeat::Int)

Construct the above ansatz, `nrepeat` is the number of measure operations,
`nlayer` is the length of each block.
"""
function twoqubit_circuit(nlayer::Int, nrepeat::Int)
    nbit_measure = nbit_virtual = 1
    nbit_used = nbit_measure + nbit_virtual
    circuit = chain(nbit_used)
    bases, measures = AbstractBlock[], Measure[]

    for i=1:nrepeat
        unit = chain(nbit_used)
        #push!(unit, put(nbit_used, 2=>H))
        for j=1:nlayer
            push!(unit, put(nbit_used, 1=>rotor(true, false)))
            push!(unit, put(nbit_used, 2=>H))
            push!(unit, put(nbit_used, 2=>Rz(0.0)))
            push!(unit, control(nbit_used, 1, 2=>shift(0.0)))
            if j == nlayer
                push!(unit, put(nbit_used, 1=>rotor(true, false)))
                push!(unit, put(nbit_used, 2=>H))
                push!(unit, put(nbit_used, 2=>Rz(0.0)))
            end
        end
        for k=(1:(i==nrepeat ? nbit_used : 1))  # the first repeat requires an extra measure
            push!(unit, chain([put(nbit_used, k=>I2)]))  # rotate basis
            push!(bases, unit[end])
            if i==nrepeat
                push!(unit, Measure(nbit_used; locs=(k,)))  # no reset in the last step
            else
                push!(unit, Measure(nbit_used; locs=(k,), resetto=0))
            end
            push!(measures, unit[end])
        end
        push!(circuit, unit)
    end
    dispatch!(circuit, :random), bases, measures
end

circuit, bases, measures = twoqubit_circuit(1, 3)

vizcircuit(circuit)


"""
    gensample(circuit, bases, measures, operator; nbatch=1024) -> bit strings

Generate samples from MPS-inspired circuit. Here, `nbatch` means nshot.
`operator` is the pauli operator to measure.
This function returns a vector of `Measure` gates, results are stored in `m.results`.
"""
function gensample(circuit, bases, measures, operator; nbatch=1024)
    E, base = YaoBlocks.eigenbasis(operator)
    for m in bases
        m[1] = put(nqubits(circuit), m[1].locs=>base')  # rotate to the operator basis
    end
    reg = zero_state(nqubits(circuit); nbatch=nbatch)
    reg |> circuit
    return [diag(mat(Float64, E))[Int.(m.results) .+ 1] for m in measures]
end

res = gensample(circuit, bases, measures, X; nbatch=1024)

# ## Model Hamiltonians
# ### Transverse field Ising Model
"""
for simplicity, we require an AbstractModel contains `size` and `periodic` members.
"""
abstract type AbstractModel{D} end

nspin(model::AbstractModel) = prod(model.size)

"""
transverse field ising model, `h` is the strength of transverse field.
"""
struct TFI{D} <:AbstractModel{1}
    size::NTuple{D, Int}
    h::Float64
    periodic::Bool
    TFI(size::Int...; h::Real, periodic::Bool) = new{length(size)}(size, Float64(h), periodic)
end



function get_bonds(model::AbstractModel{1})
    nbit, = model.size
    [(i, i%nbit+1) for i in 1:(model.periodic ? nbit : nbit-1)]
end

function get_bonds(model::AbstractModel{2})
    m, n = model.size
    cis = LinearIndices(model.size)
    bonds = Tuple{Int, Int, Float64}[]
    for i=1:m, j=1:n
        (i!=m || model.periodic) && push!(bonds, (cis[i,j], cis[i%m+1,j]))
        (j!=n || model.periodic) && push!(bonds, (cis[i,j], cis[i,j%n+1]))
    end
    bonds
end

function hamiltonian(model::TFI{1})
    nbit = nspin(model)
    sum([repeat(nbit, Z, (i,j)) for (i,j) in get_bonds(model)])*0.25 +
    sum([put(nbit, i=>X) for i=1:nbit])*0.5model.h
end

tfi_model = TFI(4; h=0.5, periodic=false)

tfi_h = hamiltonian(tfi_model)

function ising_energy(circuit, bonds, basis; nbatch=nbatch)
    results = gensample(circuit, bases, measures, basis; nbatch=nbatch)
    local eng = 0.0
    for (a, b) in bonds
        eng += mean(results[a] .* results[b])
    end
    eng/=4
end

function energy(circuit, model::TFI; nbatch=1024)
    # measuring Z
    eng = ising_energy(circuit, get_bonds(model), Z; nbatch=nbatch)
    # measuring X
    results = gensample(circuit, bases, measures, X; nbatch=nbatch)
    engx = sum(mean.(results))
    eng + model.h*engx/2
end

energy(circuit, tfi_model; nbatch=100000)

# ### Heisenberg Model
struct Heisenberg{D} <: AbstractModel{D}
    size::NTuple{D, Int}
    periodic::Bool
    Heisenberg(size::Int...; periodic::Bool) = new{length(size)}(size, periodic)
end

heisenberg_ij(nbit::Int, i::Int, j::Int=i+1) = put(nbit, i=>X)*put(nbit, j=>X) + put(nbit, i=>Y)*put(nbit, j=>Y) + put(nbit, i=>Z)*put(nbit, j=>Z)
function hamiltonian(model::Heisenberg)
    nbit = nspin(model)
    sum(x->heisenberg_ij(nbit, x[1], x[2]), get_bonds(model))*0.25
end

function energy(circuit, model::Heisenberg; nbatch=1024)
    bonds = get_bonds(model)
    sum(basis->ising_energy(circuit, bonds, basis; nbatch=nbatch), [X, Y, Z])
end

heisenberg_ij(nbit::Int, i::Int, j::Int=i+1) = put(nbit, i=>X)*put(nbit, j=>X) + put(nbit, i=>Y)*put(nbit, j=>Y) + put(nbit, i=>Z)*put(nbit, j=>Z)
function hamiltonian(model::Heisenberg)
    nbit = nspin(model)
    sum(x->heisenberg_ij(nbit, x[1], x[2]), get_bonds(model))*0.25
end

function energy(circuit, model::Heisenberg; nbatch=1024)
    bonds = get_bonds(model)
    sum(basis->ising_energy(circuit, bonds, basis; nbatch=nbatch), [X, Y, Z])
end

hei_model = Heisenberg(4; periodic=false)
energy(circuit, hei_model)

# ## Build the expanded view and check the energy
function expand_circuit(circuit)
    nbit = length(collect_blocks(Measure, circuit))
    nm = 1
    nv = 1
    c = chain(nbit)
    for (i, blk) in enumerate(circuit)
        blk = chain([b for b in blk if !(b ∈ measures || b ∈ bases)]...)
        push!(c, concentrate(nbit, blk, [(i-1)*nm+1:i*nm..., nbit-nv+1:nbit...]))
    end
    c
end

c_expand = expand_circuit(circuit)
vizcircuit(c_expand)

function wave_function(circuit)
    ec = expand_circuit(circuit)
    zero_state(nqubits(ec)) |> ec
end

@show expect(tfi_h, wave_function(circuit)) |> real
@show expect(hamiltonian(hei_model), wave_function(circuit)) |> real;

# ## Obtaining the ground state

function fidelity(circuit, VG)
    psi = zero_state(nqubits(circuit)) |> circuit
    abs(statevec(psi)' * VG)
end

nparameters(circuit)

function train_adam(circuit, model; VG=nothing, maxiter=200, α=0.3, nbatch=1024)
    # collect blocks with parameters
    rots = collect_blocks(Union{RotationGate, ControlBlock{<:ShiftGate}}, circuit)
    c_expanded = expand_circuit(circuit)
    loss_history = Float64[]
    params = vcat(parameters.(rots)...)
    grad = zero(params)
    optimizer = Optimisers.setup(Optimisers.ADAM(α), params)
    for i in 0:maxiter
        for (j,r) in enumerate(rots)
            dispatch!(+, r, (π/2,))
            E₊ = energy(circuit, model; nbatch=nbatch)
            dispatch!(-, r, (π,))
            E₋ = energy(circuit, model; nbatch=nbatch)
            dispatch!(+, r, (π/2,))
            g = 0.5(E₊ - E₋)
            grad[j] = g
        end
        Optimisers.update!(optimizer, params, grad)
        dispatch!.(rots, params)
        push!(loss_history, energy(circuit, model, nbatch=nbatch)/nspin(model))
        
        if i%10 == 0
            print("Iter $i, E/N = $(loss_history[end])")
            if !(VG isa Nothing)
                dispatch!(c_expanded, params)
                fid = fidelity(c_expanded, VG)
                println(", fidelity = $fid")
            else
                println()
            end
        end
    end
    loss_history, circuit
end

lattice_size = 4
circuit, bases, measures = twoqubit_circuit(1, lattice_size-1)
model = TFI(lattice_size; h=0.5, periodic=false)
#model = Heisenberg(lattice_size; periodic=false)
h = hamiltonian(model)

# obtain the exact ground state energy
res = eigen(mat(h) |> Matrix)
EG = res.values[1]/nspin(model)
@show EG
VG = res.vectors[:,1];

nparameters(circuit)

dispatch!(circuit, :random)
loss_history, circuit = train_adam(circuit, model; maxiter=100, VG=VG, α=0.3);

M = length(loss_history)
plot(0:M-1, [loss_history, fill(EG, M)], label=["QMPS", "Exact"], lw=3, ylabel="Energy")

circuit

parameters(circuit)