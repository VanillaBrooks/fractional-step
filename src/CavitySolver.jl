include("structs.jl")
include("indexing.jl")
include("adv.jl")
include("lap.jl")
include("gradient.jl")
include("divergence.jl")
include("lhs.jl")
include("conjugate_gradient.jl")

module CavitySolver
	using Printf
	using HDF5

	using .Main.Structs:BoundaryConditions
	using .Main.Structs:Dims
	using .Main.Structs: create_dims
	using .Main.advection: adv!
	using .Main.laplacian: lap!
	using .Main.gradient: grad!
	using .Main.divergence: div!
	using .Main.conjugate_gradient: conj_grad, debug_conjugate_gradient
	using .Main.lhs: FirstStepAx, SecondStepAx, calculate_ax
	using .Main.indexing: create_ip, create_iu_iv

	export create_dims, create_iu_iv, test, create_boundary_conditions, create_ip
	export adv, lap, grad, div_, conj_grad, BoundaryConditions, Dims
	export execute_solver, SolverResults, gradient_test

	# temp
	export debug_conjugate_gradient

	function print_vector(vec::Vector{<:AbstractFloat}, name::String)
		println("for ", name, ":")
		for i = 1:length(vec)
			println(i, " ", vec[i])
		end
	end

	function print_vector(
			vec1::Vector{<:AbstractFloat}, name1::String,
			vec2::Vector{<:AbstractFloat}, name2::String 
	)
		println(" \t", name1, "\t", name2)
		for i = 1:length(vec1)
			#println(i, "\t", vec1[i], "\t", vec2[i])
			@printf "%d \t %.5f \t %.5f\n" i vec1[i] vec2[i]
		end
	end

	function create_boundary_conditions(dims::Dims)::BoundaryConditions
		nx = dims.nx
		ny = dims.ny

		u_t = ones(nx-1);
		u_b = zeros(nx-1);
		u_r = zeros(ny);
		u_l = zeros(ny);

		v_t = zeros(nx);
		v_b = zeros(nx);
		v_r = zeros(ny-1);
		v_l = zeros(ny-1);

		return BoundaryConditions( u_t, u_b, u_r, u_l, v_t, v_b, v_r, v_l)
	end

	function zero_boundary_conditions(dims::Dims)::BoundaryConditions
		nx = dims.nx
		ny = dims.ny

		u_t = zeros(nx-1);
		u_b = zeros(nx-1);
		u_r = zeros(ny);
		u_l = zeros(ny);

		v_t = zeros(nx);
		v_b = zeros(nx);
		v_r = zeros(ny-1);
		v_l = zeros(ny-1);

		return BoundaryConditions( u_t, u_b, u_r, u_l, v_t, v_b, v_r, v_l)
	end



	struct SolverResults
		#velocity::Matrix{Float64}
	end

	function execute_solver(dims::Dims, bcs::BoundaryConditions, T::Float64, re::Float64)::SolverResults

		#iu, iv = create_iu_iv(dims)
		#ip = create_ip(dims)
		iu = Main.indexing.Iu(dims)
		iv = Main.indexing.Iv(dims)
		ip = Main.indexing.Ip(dims)

		cfl_target = 1.0

		tol::Float64 = 1E-6

		dt = calculate_dt(dims, bcs, cfl_target, re)
		#dt = 0.0391

		n_step = Int64(floor(T / dt))

		#
		# Setup Hdf5 file stuff
		#

		io_flowfield = 10
		filename = "./results/flowfield.h5"
		flowfield_file_id = h5open(filename, "w")
		num_flowfields = Int64(floor(n_step / io_flowfield))
		velocity_dataset = create_dataset(flowfield_file_id, "velocity", zeros(num_flowfields, dims.nx, dims.ny))

		zero_bcs = zero_boundary_conditions(dims)

		##
		## Define the working vectors
		##

		# velocity matricies
		uf = zeros(dims.nu)
		q_nm1 = copy(uf)
		q_n = copy(uf)
		q_np1 = copy(uf)

		# pressure matricies
		p_n = zeros(dims.np)
		# P^(n+1)
		p_np1 = copy(p_n)

		# advection matricies
		adv_n = zeros(Float64, dims.nu)
		adv_nm1 = zeros(Float64, dims.nu)

		#
		# laplace
		#

		# used in the main routine
		laplace_bc_n = zeros(Float64, dims.nu)
		# used in the LHS calculation
		laplace_nobc_n = zeros(Float64, dims.nu)

		#
		# divergence buffers
		#

		div_loop = zeros(dims.np)
		div_lhs = zeros(dims.np)

		#
		# gradients
		#
		grad_buffer_loop = zeros(dims.nu)
		grad_buffer_lhs = zeros(dims.nu)

		#
		# Steppers for LHS
		#

		first_step_lhs_calculator = FirstStepAx(
			dims, zero_bcs, iu, iv, dt, re, laplace_nobc_n
		)

		second_step_lhs_calculator = SecondStepAx(
			dims, zero_bcs, iu, iv, ip, dt, re, div_lhs, grad_buffer_lhs
		)

		println("executing ", n_step, " steps for a solver time of ", T, " with a dt = ", dt, "re = ", re)

		for step = 1:n_step
			# step variables forward once
			q_nm1 = q_n
			q_n = q_np1
			p_n = p_np1
			adv_nm1 = adv_n

			println("executing step")

			#
			# First step
			#

			lap!(dims, bcs, iu, iv, q_n, laplace_bc_n)

			adv!(dims, bcs, iu, iv, ip, q_n, adv_n)

			#print_vector(laplace_bc_n, "laplace bc", adv_n, "advection")

			# here S =  I + dt nu / 2 L
			# which becomes S = (I + dt L / (2 Re))
			#
			# the whole term is divided by dt
			sq_rhs = (
				# start of S * u_n term
				q_n / dt
				+
				(1 / (2 .* re)) .* laplace_bc_n
				# end of S * u_n term
				# start of advection term
				+
				(3/2) .* adv_n
				-
				(1/2) .* adv_nm1
			)

			rq_lhs = calculate_ax(first_step_lhs_calculator, q_n)

			#print_vector(rq_lhs, "RHS initial guess")

			uf[:] = conj_grad(first_step_lhs_calculator, rq_lhs, q_n, sq_rhs, tol)

			#if step == n_step
			#	print_vector(uf, "q_star / uf")
			#end

			#
			# Second step
			#
			
			div_loop = div!(dims, bcs, uf, iu, iv, ip, div_loop)

			second_step_rhs = div_loop ./ dt

			lhs_estimate = calculate_ax(second_step_lhs_calculator, p_n)

			p_np1 = conj_grad(second_step_lhs_calculator, lhs_estimate, p_n, second_step_rhs, tol)
			lhs_final = calculate_ax(second_step_lhs_calculator, p_np1)

			#if step == n_step 
			#	println("\n\n\nLHS final --- RHS actual")
			#	for i in 1:dims.np
			#		println(lhs_final[i], "   ", second_step_rhs[i])
			#	end
			#end

			#
			# Third fracitonal step
			#

			# u^{n+1} calculation - this is the third step
			grad!(dims, p_np1, iu, iv, ip, grad_buffer_loop)
			q_np1 = uf .- (dt .* grad_buffer_loop)

			#if step == n_step 
			#	println("\n\nfinal velocity")
			#	display(q_np1)
			#	println("")
			#end

		end

		close(flowfield_file_id)

		return SolverResults()
	end

	# find the largest dt that we are allowed to use while upholding the CFL target
	function calculate_dt(dims::Dims, bcs::BoundaryConditions, cfl_target::Float64, re::Float64)::Float64

		max_velo = maximum(bcs)

		# advective CFL = U dt / dx
		dt_advective = cfl_target * min(dims.dx, dims.dy) / max_velo

		# diffusive CFL = dt / (Re * dx^2)
		dt_diffusive = re * min(dims.nx, dims.ny)^2 * cfl_target

		dt = min(dt_diffusive, dt_advective)

		return dt
	end

	function gradient_test(dims::Dims)
		pressure = zeros(dims.nx * dims.ny)
		expected = zeros(dims.nu)

		iu, iv = create_iu_iv(dims)
		ip = create_ip(dims)

		for i = 1:dims.nx
			for j = 1:dims.ny
				pressure[ip[i,j]] = sin(i) + cos(j)
			end
		end

		no_bcs = zero_boundary_conditions(dims)

		out = grad(dims, pressure, iu, iv, ip)

		println("gradient is ")
		display(transpose(out))

		println("expected is ")
		display(transpose(expected))


	end

end

#function debug_solver(dims::CavitySolver.Dims, bcs::BoundaryConditions)
#
#	iu, iv = create_iu_iv(dims)
#	ip = create_ip(dims)
#
#	q = zeros(dims.nu)
#	p = ones(dims.np)
#
#	a, b = size(iu)
#	c, d = size(iv)
#	e = size(q)
#	show(a*b + c*d)
#	println("")
#	show(e)
#	println("")
#
#	adv(dims, bcs, iu, iv, ip, q)
#	lap(dims, bcs, iu, iv, q)
#	grad(dims, p, iu, iv, ip)
#
#	div_(dims, bcs, q, iu, iv, ip)
#end

#n = 16

#using .CavitySolver
#const dims = create_dims(1.0, 1.0, n,n)
#bcs = create_boundary_conditions(dims)

#debug_solver(dims, bcs)
#Re = 10.0
#time_end = 0.01
time_end = 0.1
# Re = 100
# T = 1
# steps = 128
# runtime = 355 seconds
#@time 
#execute_solver(dims, bcs, time_end, Re)

#gradient_test(dims)

#debug_conjugate_gradient()

#using .indexing
#check_indexing()
