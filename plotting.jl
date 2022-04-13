### A Pluto.jl notebook ###
# v0.18.4

using Markdown
using InteractiveUtils

# ╔═╡ b42336e6-58ac-42d9-a6d9-b1dd07b5c3d4
begin
	using Revise
	using Pkg
	#Pkg.activate("/home/brooks/github/Gradients")
	Pkg.activate()
	#Pkg.develop("Gradients")
end

# ╔═╡ db488b20-9896-4ceb-82f3-cf49b39d8d92
begin
	using HDF5
	using Plots
end

# ╔═╡ 615a0698-e463-43c4-a425-46deed41355a
using Gradients

# ╔═╡ 5561e533-5fb2-46c5-ab7c-dc357c8c1f53
Gradients.test2()

# ╔═╡ 74e73576-59c8-4375-97bf-6890199a9fb8
function save_png(path)
	save = "/home/brooks/fractional-step/paper/figs/" * path
	png(save)
end

# ╔═╡ a30c733d-9cf9-4f90-9926-e4efba1ce0a0
test2()

# ╔═╡ 5a329021-320c-41ed-a28d-8b98bfbbb745
#data_path = "/home/brooks/fractional-step/results/flowfield.h5"
data_path = "/home/brooks/fractional-step/paper_results/re_1000_40_sec.h5"

# ╔═╡ 70a8481c-3373-4623-a207-47759a24e9f9
h5file = h5open(data_path, "r")

# ╔═╡ b0cf1d82-2d69-4d5a-bd37-151eec44ab26
begin
	velocity = read(h5file["velocity"])
	velocityT = permutedims(velocity, (1,2,4,3))
	time = read(h5file["times"])
end

# ╔═╡ a4bbd6a3-f09b-473f-a758-82bbb9fbfb5b
num_writes, _, nx, ny = size(velocity)

# ╔═╡ 4d723360-c404-4eb8-b980-35ed65d0e808
size(velocityT[1, 2, :, :])

# ╔═╡ 6ae07f94-4495-4b55-85dd-84b42608f945
begin
	vorticity = zeros(num_writes, nx, ny)
	dx = 1 / nx
	dy = 1/ ny

	for write = 1:num_writes
		dV_dX = df_dx(velocity[write, 2, :, :], dx) 
		dX_dY = df_dy(velocity[write, 1, :, :], dy)
		# 
		vorticity[write, :, :] = dV_dX - dX_dY
	end

	vorticity = permutedims(vorticity, (1,2,3))
	# vorticity = reverse(vorticity; dims=2)
	vorticity = reverse(vorticity; dims=3)
	vorticity = reverse(vorticity; dims=2)

end

# ╔═╡ 360aadb3-2090-4eb7-98c9-ed6fe9bdce39
typeof(h5file["velocity"])

# ╔═╡ 19abf3ed-85de-467a-8bb5-09a63fe0c555
velocity[end, 1, :, nx]

# ╔═╡ c3a7aa30-63a4-4bef-96b2-04addc6087f2
begin
	local min_poss = minimum(velocity[:, 1, :, :])
	local max_poss = maximum(velocity[:, 1, :, :])
	bound_xvelo = max(abs(min_poss), abs(max_poss))
	
	local anim = @animate for write = 1:num_writes
		heatmap(
			velocityT[write, 1, :, :],
			title = "u velocity",
			clim = (-bound_xvelo, bound_xvelo),
			c = :bluesreds
		)
	end
	gif(anim)
end

# ╔═╡ 9c653171-5be3-4d06-a098-3242e44ecb99
begin
	local min_poss = minimum(velocityT[:, 2, :, :])
	local max_poss = maximum(velocityT[:, 2, :, :])
	bound_yvelo = max(abs(min_poss), abs(max_poss))
	# bound = .3
	local anim = @animate for write = 1:num_writes
		heatmap(
			velocityT[write, 2, :, :],
			title = "v velocity",
			clims = (-bound_yvelo, bound_yvelo),
			c = :bluesreds,
		)
	end
	gif(anim)
end

# ╔═╡ beb8c5db-5ce2-4c1a-a061-34df027e3da1
begin
	idx = num_writes

	local plt = contour(
		1:nx,
		1:ny,
		vorticity[idx, :, :],
		levels=100,
		# levels = -2:0.05:2,
		title = "Vorticity for Re = 1000 at t = 40 sec",
		dpi = 200,
		margin=5Plots.mm,
		size = (600,600)
	)

	save_png("vorticity_brooks.png")

	plt
end

# ╔═╡ c42099e1-e5e0-4ce1-8b6f-f9dc42ea0fa1
begin
	local anim = @animate for write = 1:num_writes
		local hm1 = heatmap(
			velocityT[write,1, :, :],
			title = "x velocity",
			clim = (-bound_xvelo, bound_xvelo),
			c = :bluesreds
		)

		local hm2 = heatmap(
			velocityT[write,2, :, :],
			title = "y velocity",
			clim = (-bound_yvelo, bound_yvelo),
			c = :bluesreds
		)
		plot(hm1, hm2, layout = (1,2), size = (1000, 400))
	end
	gif(anim)
end

# ╔═╡ 497f564f-eac2-4a1d-a9b7-011d1d9ad107
begin
	quiver_slice = 1:5:(ny*nx)
	x_grid = repeat(1:nx, 1, ny)[quiver_slice] |> vec
	y_grid = transpose(repeat(1:nx, 1,ny))[quiver_slice] |> vec
end

# ╔═╡ c6ffcd1a-af89-417e-9f44-59a846c0d22e
begin
	local step = num_writes
	local vx = velocity[step, 1, :, :]
	local vy = velocity[step, 2, :, :]
	
	local vx = reshape(vx, nx*ny)[quiver_slice]
	local vy = reshape(vy, nx*ny)[quiver_slice]
	
	local plt = quiver(
		x_grid,
		y_grid,
		quiver=(vx, vy),
		dpi = 200,
		title="Velocity Field at t = 40 sec, Re = 1000",
		margin = 5Plots.mm,
		xlabel = "x",
		ylabel = "y",
		aspect_ratio=:equal,
		size = (600, 600)
	)

	save_png("quiver_brooks.png")

	plt
end

# ╔═╡ 0f2c7918-5e1e-4969-a1d9-7289fae089f2
x_grid

# ╔═╡ 2269826d-0800-461a-bdaa-cb484c82e189
y_grid

# ╔═╡ 5cbc7c56-b36c-4a36-818c-bf86dba7b8f3
begin
	local anim = @animate for write = 1:num_writes
		local vx = velocity[write, 1, :, :]
		local vy = velocity[write, 2, :, :]

		slice = 1:5:(ny*nx)
	
		local vx = reshape(vx, nx*ny)[quiver_slice]
		local vy = reshape(vy, nx*ny)[quiver_slice]

		# local gridx = x_grid[slice]
		# local gridy = y_grid[slice]
		
		quiver(
			x_grid,
			y_grid,
			quiver=(vx, vy),
			size=(500,500),
			dpi=200,
			xlabel = "x",
			ylabel = "y",
			title = "velocity evolution over time",
			margin=5Plots.mm
		)
	end

	gif(anim)
end

# ╔═╡ Cell order:
# ╠═b42336e6-58ac-42d9-a6d9-b1dd07b5c3d4
# ╠═db488b20-9896-4ceb-82f3-cf49b39d8d92
# ╠═615a0698-e463-43c4-a425-46deed41355a
# ╠═5561e533-5fb2-46c5-ab7c-dc357c8c1f53
# ╠═74e73576-59c8-4375-97bf-6890199a9fb8
# ╠═a30c733d-9cf9-4f90-9926-e4efba1ce0a0
# ╠═5a329021-320c-41ed-a28d-8b98bfbbb745
# ╠═70a8481c-3373-4623-a207-47759a24e9f9
# ╠═b0cf1d82-2d69-4d5a-bd37-151eec44ab26
# ╠═a4bbd6a3-f09b-473f-a758-82bbb9fbfb5b
# ╠═4d723360-c404-4eb8-b980-35ed65d0e808
# ╠═6ae07f94-4495-4b55-85dd-84b42608f945
# ╠═360aadb3-2090-4eb7-98c9-ed6fe9bdce39
# ╠═19abf3ed-85de-467a-8bb5-09a63fe0c555
# ╠═c3a7aa30-63a4-4bef-96b2-04addc6087f2
# ╠═9c653171-5be3-4d06-a098-3242e44ecb99
# ╠═beb8c5db-5ce2-4c1a-a061-34df027e3da1
# ╠═c42099e1-e5e0-4ce1-8b6f-f9dc42ea0fa1
# ╠═497f564f-eac2-4a1d-a9b7-011d1d9ad107
# ╠═c6ffcd1a-af89-417e-9f44-59a846c0d22e
# ╠═0f2c7918-5e1e-4969-a1d9-7289fae089f2
# ╠═2269826d-0800-461a-bdaa-cb484c82e189
# ╠═5cbc7c56-b36c-4a36-818c-bf86dba7b8f3
