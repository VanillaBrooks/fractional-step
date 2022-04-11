module lhs
	using .Main.Structs: Dims, BoundaryConditions
	using .Main.laplacian: lap
	using .Main.gradient: grad
	using .Main.divergence: div_
	using .Main.indexing: Indexable

	export AxCalculator, calculate_ax, FirstStepAx, SecondStepAx, ConstantMatrix

	# all AxCalculator subtypes must implement calculate_ax function so they may be used in the
	# conjugate gradient formula
	abstract type AxCalculator end

	#default impl
	function calculate_ax(context::T, x::Vector{Float64})::Vector{Float64} where T <: AxCalculator
		error("calculate_ax not implemented for this type")
	end

	struct FirstStepAx{IU, IV} <: AxCalculator
		dims::Dims
		zero_bcs::BoundaryConditions
		iu::IU
		iv::IV
		dt::Float64
		re::Float64
		lap_nobc_n::Vector{Float64}
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
		#ctx.lap_nobc_n = 
		lap(ctx.dims, ctx.zero_bcs, ctx.iu, ctx.iv, x, ctx.lap_nobc_n)
		q = (
			x / ctx.dt
			-
			(1 / (2 * ctx.re)) * ctx.lap_nobc_n
		)
		return q
	end

	struct SecondStepAx{IU, IV, IP} <: AxCalculator
		dims::Dims
		zero_bcs::BoundaryConditions
		iu::IU
		iv::IV
		ip::IP
		dt::Float64
		re::Float64
		div_buffer::Vector{Float64}
	end

	function calculate_ax(ctx::SecondStepAx, p::Vector{Float64})::Vector{Float64}
		# take the gradient of the pressure
		G = grad(ctx.dims, p, ctx.iu, ctx.iv, ctx.ip)

		# then take the divergence of R^-1
		#ctx.div_buffer = 
		div_(ctx.dims, ctx.zero_bcs, G, ctx.iu, ctx.iv, ctx.ip, ctx.div_buffer)

		return ctx.div_buffer
	end

	struct ConstantMatrix <: AxCalculator
		A::Matrix{Float64}
	end

	function calculate_ax(ctx::ConstantMatrix, x::Vector{Float64})::Vector{Float64}
		ctx.A * x
	end
end
