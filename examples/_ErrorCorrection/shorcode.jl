### A Pluto.jl notebook ###
# v0.15.1

using Markdown
using InteractiveUtils

# ╔═╡ 18616886-e1af-460d-8ded-e50c65b9d2aa
begin
	using Pkg
	Pkg.activate()
	using Yao, YaoPlots, Test
end

# ╔═╡ 68efde8a-0127-11ec-3844-c56fbe9cb384
md"# Error correction"

# ╔═╡ 72457507-dce0-4b21-98c7-1fbb8d258fdd
md"## Shor's code"

# ╔═╡ 86b9519a-0c18-4709-ac30-0a1783df07d9
shor_code = chain(9, cnot(1,4), cnot(1, 7), repeat(H, (1,4,7)), cnot(1, 2), cnot(1,3), cnot(4,5), cnot(4,6), cnot(7,8), cnot(7,9))

# ╔═╡ 85f00eab-8fbe-4116-ab6c-022f00390f45
vizcircuit(shor_code, scale=0.5)

# ╔═╡ 385c7b9c-77b5-47b4-803b-ffab77b79047
@testset "zero state and one state" begin
	pos = (product_state(bit"000") + product_state(bit"111"))/sqrt(2)
	neg = (product_state(bit"000") - product_state(bit"111"))/sqrt(2)
	@test zero_state(9) |> shor_code ≈ join(pos, pos, pos)
	@test zero_state(9) |> put(9, 1=>X) |> shor_code ≈ join(neg, neg, neg)
end

# ╔═╡ ed4d93a8-d846-414b-8821-dfb9cc716d6f
function shor_encode(reg)
	join(zero_state(8), reg) |> shor_code
end

# ╔═╡ 8e79ce69-0ea0-4f92-8b16-4867b2907de1
function shor_decode!(reg::AbstractRegister)
	res = measure_remove!(reg |> shor_code', 2:9)
	#@assert Int(res) == 0
	return reg
end

# ╔═╡ a0b3940c-ee00-4ded-9a8f-b084e8e3afa1
@testset "encode and decode" begin
	ψ = rand_state(1)
	cp = shor_encode(ψ)
	ϕ = shor_decode!(cp)
	@test fidelity(ϕ, ψ) ≈ 1
end

# ╔═╡ 0a819cc4-02dd-4a26-bc9a-13e6cef8065a
md"## Error syndrome detection"

# ╔═╡ 994349b6-2da5-44d4-a3e0-f88852593ad5
md"Should not change the information that one wants to encode. i.e. The measurement result does not tell any infromation about the original qubit state."

# ╔═╡ fe7232b3-8b48-4c6b-bec3-49548c005de2
function shor_syndrome!(reg)
	# Flip error (X)
	xerror, zerror = 0, 0
	
	for j=1:3
		offset = (j-1)*3
		z1 = Float64(measure!(repeat(9, Z, (1+offset,2+offset)), reg)) < 0
		z2 = Float64(measure!(repeat(9, Z, (2+offset,3+offset)), reg)) < 0
		xerror += (z1 && !z2) << offset      # 1st qubit errors
		xerror += (z1 && z2) << (offset+1)      # 2nd qubit errors
		xerror += (!z1 && z2) << (offset+2)     # 3rd qubit errors
	end
	# Phase error (Z)
	x1 = Float64(measure!(repeat(9, X, 1:6), reg)) < 0
	x2 = Float64(measure!(repeat(9, X, 4:9), reg)) < 0
	
	zerror += (x1 && !x2)      			# a qubit in 1st block errors
	zerror += (x1 & x2) << 3     		# a qubit in 2nd block errors
	zerror += (!x1 & x2) << 6    		# a qubit in 3rd block errors
	return xerror, zerror
end

# ╔═╡ 24124d5c-ed5a-43a3-8cf1-209ca374917c
@testset "syndrome" begin
	reg1 = zero_state(1)
	reg = shor_encode(reg1)
	for i=1:9
		regi = copy(reg) |> put(9, i=>X)
		sx, sz = shor_syndrome!(regi)
		@test sx == 1<<(i-1)
		@test sz == 0
	end
	for i=1:9
		regi = copy(reg) |> put(9, i=>Z)
		sx, sz = shor_syndrome!(regi)
		@test sx == 0
		@test sz == 1<<((i-1)÷3*3)
	end
end

# ╔═╡ ecda6bcd-04b4-40ba-a1ef-79bbbe8f60b8
md"## Recovery"

# ╔═╡ 7e13f518-9089-4efa-a093-ca3e2dd5df7e
function shor_recover!(reg::AbstractRegister, syndrome_x, syndrome_z)
	n = nqubits(reg)
	for i in 1:n
		YaoArrayRegister.readbit(syndrome_x, i) == 1 && apply!(reg, put(n, i=>X))
		YaoArrayRegister.readbit(syndrome_z, i) == 1 && apply!(reg, put(n, i=>Z))
	end
	return reg
end

# ╔═╡ e391b164-6470-4fc3-85da-fddd17278080
@testset "encode and decode" begin
	ψ2 = rand_state(1)
	# encode
	cp2 = shor_encode(ψ2)
	# ZX error on 4
	cp2 |> put(9, 4=>Z*X)
	syndrome_x, syndrome_z = shor_syndrome!(cp2)
	@test syndrome_x==1<<3 && syndrome_z==1<<3
	shor_recover!(cp2, syndrome_x, syndrome_z)
	ϕ2 = shor_decode!(cp2)
	@test fidelity(ϕ2, ψ2) ≈ 1
end

# ╔═╡ df4f7012-8ce4-46bc-9eec-acdd681ce8fa
md"### Stablizer"

# ╔═╡ 8942e4be-dd1f-41a4-9d36-9186dfe063f1
shor_stablizers = [[repeat(9,Z, (i,i+1)) for i in [1,2,4,5,7,8]]..., repeat(9, X, 1:6), repeat(9, X,4:9)]

# ╔═╡ 1dc1c287-3c91-4e59-a011-cb8ccf95385b
@testset "shor stablizer" begin
	zgate = repeat(9, X,1:9)
	xgate = repeat(9, Z,1:9)
	for reg in [zero_state(1), (zero_state(1) |> X)]
		REG = shor_encode(reg)
		for g in shor_stablizers
			@test copy(REG) |> g ≈ REG
		end
		@test fidelity(copy(REG) |> zgate, shor_encode(reg |> Z)) ≈ 1
		@test fidelity(copy(REG) |> xgate, shor_encode(reg |> X)) ≈ 1
	end
end

# ╔═╡ e15dd67a-5a7c-436a-8110-67bcca59867e
md"## Steane code"

# ╔═╡ 49e4aaa8-f0d4-42e8-aaf0-307596d2e921
gs = [kron(I2,I2,I2,X,X,X,X), kron(I2,X,X,I2,I2,X,X), kron(X,I2,X,I2,X,I2,X), kron(I2,I2,I2,Z,Z,Z,Z), kron(I2,Z,Z,I2,I2,Z,Z), kron(Z,I2,Z,I2,Z,I2,Z)]

# ╔═╡ ab31a2dd-cb1b-4540-b478-d5a5fe011058
steane_zero = sum(x->product_state(x), [bit"0000000", bit"1010101", bit"0110011", bit"1100110", bit"0001111", bit"1011010", bit"0111100", bit"1101001"])

# ╔═╡ 00ad35e8-2b6f-46bd-9c5b-0f763d58785e
steane_one = sum(x->product_state(x), [bit"1111111", bit"0101010", bit"1001100", bit"0011001", bit"1110000", bit"0100101", bit"1000011", bit"0010110"]) |> normalize!

# ╔═╡ 7a102ede-e829-4599-9465-c4dc6b16cdb2
@testset "steane code" begin
	for g in gs
		@test steane_zero |> g ≈ steane_zero
		@test steane_one |> g ≈ steane_one
	end
end

# ╔═╡ Cell order:
# ╟─68efde8a-0127-11ec-3844-c56fbe9cb384
# ╟─72457507-dce0-4b21-98c7-1fbb8d258fdd
# ╠═18616886-e1af-460d-8ded-e50c65b9d2aa
# ╠═86b9519a-0c18-4709-ac30-0a1783df07d9
# ╠═85f00eab-8fbe-4116-ab6c-022f00390f45
# ╠═385c7b9c-77b5-47b4-803b-ffab77b79047
# ╠═ed4d93a8-d846-414b-8821-dfb9cc716d6f
# ╠═8e79ce69-0ea0-4f92-8b16-4867b2907de1
# ╠═a0b3940c-ee00-4ded-9a8f-b084e8e3afa1
# ╟─0a819cc4-02dd-4a26-bc9a-13e6cef8065a
# ╟─994349b6-2da5-44d4-a3e0-f88852593ad5
# ╠═fe7232b3-8b48-4c6b-bec3-49548c005de2
# ╠═24124d5c-ed5a-43a3-8cf1-209ca374917c
# ╟─ecda6bcd-04b4-40ba-a1ef-79bbbe8f60b8
# ╠═7e13f518-9089-4efa-a093-ca3e2dd5df7e
# ╠═e391b164-6470-4fc3-85da-fddd17278080
# ╟─df4f7012-8ce4-46bc-9eec-acdd681ce8fa
# ╠═8942e4be-dd1f-41a4-9d36-9186dfe063f1
# ╠═1dc1c287-3c91-4e59-a011-cb8ccf95385b
# ╟─e15dd67a-5a7c-436a-8110-67bcca59867e
# ╠═49e4aaa8-f0d4-42e8-aaf0-307596d2e921
# ╠═ab31a2dd-cb1b-4540-b478-d5a5fe011058
# ╠═00ad35e8-2b6f-46bd-9c5b-0f763d58785e
# ╠═7a102ede-e829-4599-9465-c4dc6b16cdb2
