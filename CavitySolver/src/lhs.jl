module lhs
	using .Main.Structs: Dims, BoundaryConditions
	using .Main.laplacian: lap
	using .Main.gradient: grad
	using .Main.divergence: div_

	export AxCalculator, calculate_ax, FirstStepAx, SecondStepAx, ConstantMatrix

	# all AxCalculator subtypes must implement calculate_ax function so they may be used in the
	# conjugate gradient formula
	abstract type AxCalculator end

	#default impl
	function calculate_ax(context::T, x::Vector{Float64})::Vector{Float64} where T <: AxCalculator
		error("calculate_ax not implemented for this type")
	end

	struct FirstStepAx <: AxCalculator
		dims::Dims
		zero_bcs::BoundaryConditions
		iu::Matrix{Int64}
		iv::Matrix{Int64}
		dt::Float64
		re::Float64
	end

	# calculation of full LHS of the _first step_
	function calculate_ax(ctx::FirstStepAx, x::Vector{Float64})::Vector{Float64}
		# equivalent to R * u_f
		# because R = I - dt * nu * L / 2
		# so we change to R = I - L / (2 re) 
		# after multiplying by u_f then we have 
		# R * u_f = u_f - L*u_f / (2 * re)
		#
		# the whole term is effectivly divided by dt
		lap_nobc_n = lap(ctx.dims, ctx.zero_bcs, ctx.iu, ctx.iv, x)
		q = (
			x / ctx.dt
			-
			(1 / (2 * ctx.re)) * lap_nobc_n
		)
		return q
	end

	struct SecondStepAx <: AxCalculator
		dims::Dims
		zero_bcs::BoundaryConditions
		iu::Matrix{Int64}
		iv::Matrix{Int64}
		ip::Matrix{Int64}
		dt::Float64
		re::Float64
	end

	function calculate_ax(ctx::SecondStepAx, p::Vector{Float64})::Vector{Float64}
		# take the gradient of the pressure
		G = grad(ctx.dims, p, ctx.iu, ctx.iv, ctx.ip)

		# then take the divergence of R^-1
		div = div_(ctx.dims, ctx.zero_bcs, G, ctx.iu, ctx.iv, ctx.ip)

		return div
	end

	struct ConstantMatrix <: AxCalculator
		A::Matrix{Float64}
	end

	function calculate_ax(ctx::ConstantMatrix, x::Vector{Float64})::Vector{Float64}
		ctx.A * x
	end
end
