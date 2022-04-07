module conjugate_gradient
	using .Main.Structs: Dims, BoundaryConditions
	using .Main.laplacian: lap
	using Printf

	export conj_grad, calculate_ax, FirstStepAx, debug_conjugate_gradient

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
		q = 
			x / ctx.dt
			-
			(1 / (2 * ctx.re)) * lap_nobc_n
		return q
	end

	struct ConstantMatrix <: AxCalculator
		A::Matrix{Float64}
	end

	function calculate_ax(ctx::ConstantMatrix, x::Vector{Float64})::Vector{Float64}
		ctx.A * x
	end

	#jfunction debug_conjugate_gradient(A::Matrix{Float64}, x::Matrix{Float64}, b::Matrix{Float64})
	function debug_conjugate_gradient()::Vector{Float64}
		A = Float64.([
			5 -2 0;
			-2 5 1;
			0 1 5
		])

		b = Float64.([
			 20;10; -10;
		])

		println("A matrix is ")
		display(A)
		#println("\nx matrix is ")
		#display(x)
		println("\nb matrix is ")
		display(b)
		#println("\nAx (should be equivalen to b) is")
		#display(A*x)
		println("")

		x_guess= zeros(3)

		step_context = ConstantMatrix(A)
		#Ax = A*x
		Ax_guess = A*x_guess

		estimated_x = conj_grad(step_context, Ax_guess, x_guess, b, 0.00001)

		println("estimated matrix x is:")
		display(estimated_x)

		println("\nA*estimated x is:")
		display(A * estimated_x)

		return estimated_x
	end

	function conj_grad(
		step_context::T,
		Ax::Vector{Float64}, 
		x_guess::Vector{Float64},
		b::Vector{Float64}, 
		tol::Float64,
	)::Vector{Float64} where T <: AxCalculator
		i = 0
		# r should be a effectively vector 
		r = b - Ax

		d = r

		delta_new::Float64 = transpose(r) * r
		delta_0::Float64 = delta_new

		imax = length(b)

		x = x_guess

		while i < imax && (delta_new > tol^2 * delta_0 || i == 0)
			# q = Ad
			q = calculate_ax(step_context, d)

			# alpha
			alpha = delta_new / (transpose(r) * q)[1,1]

			x = x + (alpha * d) #TODO: This might need to be d instead of r

			if i % 50 == 0 && i != 0
				# calculate a new Ax 
				Ax_new = calculate_ax(step_context, x)

				# update r using this new calculation to avoid floating points errors
				r = b - Ax_new
			else
				r = r - alpha * q
			end

			delta_old = delta_new
			delta_new = transpose(r) * r

			beta = delta_new / delta_old

			d = r + beta * d

			i += 1
		end

		println("conjugate gradient ran ", i, " steps")

		return x
	end
end
