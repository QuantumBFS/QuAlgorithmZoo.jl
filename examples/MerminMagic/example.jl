### A Pluto.jl notebook ###
# v0.19.18

using Markdown
using InteractiveUtils

# ╔═╡ a1d7da56-1810-4ae0-ad52-10206807d287
begin
	import Pkg;
	Pkg.add("Yao")
	Pkg.add("YaoPlots")
end

# ╔═╡ 667cd31a-7e79-11ed-2f5c-9be185e115e2
include("./main.jl")

# ╔═╡ 385a00a2-fcab-4ec1-b22f-0680ae10893a
for i in 1:10
    play()
end

# ╔═╡ Cell order:
# ╠═a1d7da56-1810-4ae0-ad52-10206807d287
# ╠═667cd31a-7e79-11ed-2f5c-9be185e115e2
# ╠═385a00a2-fcab-4ec1-b22f-0680ae10893a
