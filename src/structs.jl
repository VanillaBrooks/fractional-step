module Structs
	export BoundaryConditions
	export Dims
	export create_dims

	struct BoundaryConditions
		u_t::Vector{Float64}
		u_b::Vector{Float64}
		u_r::Vector{Float64}
		u_l::Vector{Float64}
		# v direction
		v_t::Vector{Float64}
		v_b::Vector{Float64}
		v_r::Vector{Float64}
		v_l::Vector{Float64}
	end

	Base.maximum(bcs::BoundaryConditions)::Float64 = 
		max(
			maximum(bcs.u_t),
			maximum(bcs.u_b),
			maximum(bcs.u_r),
			maximum(bcs.u_l),
			maximum(bcs.v_t),
			maximum(bcs.v_b),
			maximum(bcs.v_r),
			maximum(bcs.v_l)
		)

	struct Dims
		nx::Int
		ny::Int
		np::Int
		nu::Int
		dx::Float64
		dy::Float64
	end

	function create_dims(length::Float64, height::Float64, nx::Int, ny::Int)::Dims 
		np = nx * ny
		nu = (nx - 1) * ny + (ny - 1) * nx

		dx = length / nx
		dy = height / ny
		return Dims(nx, ny, np, nu, dx, dy)
	end

end
