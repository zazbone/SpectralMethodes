### A Pluto.jl notebook ###
# v0.19.11

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 9a33fdc0-1be6-11ed-280d-73b9d4da9a14
begin
	using Pkg
	Pkg.activate(".")
	using Plots
	using PlutoUI
	using Spectral
	using FFTW: fft, fftfreq, ifft
	Pkg.status()
end

# ╔═╡ 36ed3bc6-b314-44d3-8ac5-2741d9a0984e
begin
	MINX = -1.
	MAXX = 1.
	NX = 1024
	T = 1.
	NT = 1024
	ν = -0.1
end;

# ╔═╡ 934db463-d01a-4339-8dcb-cf296fcb4737
md"""
### Choose the initial condition
Initial condition function define by:

$(@bind u0 Select(["guassian", "signal"]))(x, $(@bind arg1 Slider(-1:0.1:1, show_value = true, default=0)), $(@bind arg2 Slider(-1:0.1:1, show_value = true, default=0)))
"""

# ╔═╡ 0330c70e-6605-41ce-9823-dbb67aac8f85
begin
	u_full = zeros((NX, NT))
	u_full[:, 1] = distrib_switch(u0).(range(MINX, stop=MAXX, length=NX), arg1, arg2)
	du = zeros(ComplexF64, NX)
	ddu = zeros(ComplexF64, NX)
end;

# ╔═╡ 6624e383-3e1c-4961-8f32-0d82c2a52011
dt = T / NT

# ╔═╡ 24350072-8997-400e-a281-29a0ad8f2219
ik = 1im * fftfreq(NX, NX) / (MAXX - MINX)

# ╔═╡ e46c3825-3bac-44fb-91e2-1a856de42c1a
k2 = ik .^ 2

# ╔═╡ c9e82d1b-47de-499c-ac79-b4b19bee0823
begin
	plot(range(-2π, 2π, NX), -sin.(range(-2π, 2π, NX)))
	plot!(range(-2π, 2π, NX), real.(ifft(k2 .* fft(sin.(range(-2π, 2π, NX))))))
end

# ╔═╡ 8f145e47-cc73-499b-956f-f43141e1c252
@time for i in 1:(NX - 1)
	du[:] = fft(u_full[:, i] .+ .0im)
	ddu[:] = fft(du)
	du[:] = du .* ik
	ddu[:] = ddu .* k2
	u_full[:, i + 1] = u_full[:, i] + dt * u_full[:, i] .* real.(ifft(du)) + ν * dt * real.(ifft(ddu))
end

# ╔═╡ 8473a4f9-3e07-4dcc-86d4-62c89ed5effa
plot(range(MINX, MAXX, NX), u_full[:, 2])

# ╔═╡ 3d78b8c6-fc15-40b6-80dc-ae5b6ec6efce
heatmap(u_full[:, 1:1])

# ╔═╡ Cell order:
# ╠═9a33fdc0-1be6-11ed-280d-73b9d4da9a14
# ╠═36ed3bc6-b314-44d3-8ac5-2741d9a0984e
# ╟─934db463-d01a-4339-8dcb-cf296fcb4737
# ╠═0330c70e-6605-41ce-9823-dbb67aac8f85
# ╠═6624e383-3e1c-4961-8f32-0d82c2a52011
# ╠═24350072-8997-400e-a281-29a0ad8f2219
# ╠═e46c3825-3bac-44fb-91e2-1a856de42c1a
# ╠═c9e82d1b-47de-499c-ac79-b4b19bee0823
# ╠═8f145e47-cc73-499b-956f-f43141e1c252
# ╠═8473a4f9-3e07-4dcc-86d4-62c89ed5effa
# ╠═3d78b8c6-fc15-40b6-80dc-ae5b6ec6efce
