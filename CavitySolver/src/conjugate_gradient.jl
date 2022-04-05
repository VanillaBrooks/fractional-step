module conjugate_gradient
	using .Main.Structs: Dims, BoundaryConditions
	using .Main.laplacian: lap
	using Printf

	struct FirstStep
	end

	function conj_grad(
		step::FirstStep, 
		dims::Dims, 
		zero_bcs::BoundaryConditions, 
		iu::Matrix{Int64}, 
		iv::Matrix{Int64}, 
		dt::Float64,
		Ax::Matrix{Float64}, 
		x_guess::Matrix{Float64}
		b::Matrix{Float64}, 
		tol::Float64
	)::Matrix{Float64}
		i = 0
		# r should be a effectively vector 
		r = b - Ax

		d = r

		delta_new = transpose(r) * r
		delta_0 = delta_new

		imax = length(b)

		x = x_guess

		while i < imax && delta_new > tol^2 * delta_0 
			#
			# q = Ad
			#
			lap_nobc_n = lap(dims, zero_bcs, iu, iv, d)
			q = 
				d / dt
				-
				(1 / (2 * re)) * lap_nobc_n

			# alpha
			alpha = delta / (transpose(r) * q)

			# calculate new x
			x = x + alpha * r

			if i % 50 == 0
				#
				# calculate a new Ax 
				# 
				lap_nobc_n = lap(dims, zero_bcs, iu, iv, x)
				Ax_new = 
					x / dt
					-
					(1 / (2 * re)) * lap_nobc_n

				# update r using this new calculation to avoid floating points errors
				r = b - Ax_new
			else
				r = r - alpha * q
			end

			delta_old = delta_new
			delta_new = transpose(r) * r

			beta = delta_new / delta_old

			d = r + beta * d

			delta = transpose(r) * r

			i += 1
		end

		return x
	end
end
