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
end
