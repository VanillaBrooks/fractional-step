module Structs
	export BoundaryConditions
	export Dims

	struct BoundaryConditions
		u_t::Matrix{Float64}
		u_b::Matrix{Float64}
		u_r::Matrix{Float64}
		u_l::Matrix{Float64}
		# v direction
		v_t::Matrix{Float64}
		v_b::Matrix{Float64}
		v_r::Matrix{Float64}
		v_l::Matrix{Float64}
	end

	struct Dims
		nx::Int
		ny::Int
		np::Int
		nu::Int
		dx::Float64
		dy::Float64
	end
end
