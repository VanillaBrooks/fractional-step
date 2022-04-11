module conjugate_gradient
	using .Main.Structs: Dims, BoundaryConditions
	using .Main.laplacian: lap
	using .Main.lhs: calculate_ax, AxCalculator, ConstantMatrix
	using Printf

	export conj_grad, debug_conjugate_gradient

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
			alpha = delta_new / (transpose(r) * q)

			x = x .+ alpha .* d

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

			d[:] = r .+ beta .* d

			i += 1
		end

		#println("\nconjugate gradient ran ", i, " steps")

		return x
	end
end
