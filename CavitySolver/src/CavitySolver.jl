include("structs.jl")
include("adv.jl")
include("lap.jl")
include("gradient.jl")
include("divergence.jl")
include("conjugate_gradient.jl")

module CavitySolver
	using Printf

	using .Main.Structs:BoundaryConditions
	using .Main.Structs:Dims
	using .Main.advection: adv
	using .Main.laplacian: lap
	using .Main.gradient: grad
	using .Main.divergence: div_
	using .Main.conjugate_gradient: conj_grad, FirstStepAx, calculate_ax, debug_conjugate_gradient

	export create_dims, create_iu_iv, test, create_boundary_conditions, create_ip
	export adv, lap, grad, div_, conj_grad, BoundaryConditions, Dims
	export execute_solver, SolverResults

	# temp
	export debug_conjugate_gradient

	function create_dims(length::Float64, height::Float64, nx::Int, ny::Int)::Dims 
		np = nx * ny
		nu = (nx - 1) * ny + (ny - 1) * nx

		dx = length / nx
		dy = height / ny
		return Dims(nx, ny, np, nu, dx, dy)
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

		u_t = zeros(nx-1, 1);
		u_b = zeros(nx-1, 1);
		u_r = zeros(ny, 1);
		u_l = zeros(ny, 1);

		v_t = zeros(nx, 1);
		v_b = zeros(nx, 1);
		v_r = zeros(ny-1, 1);
		v_l = zeros(ny-1, 1);

		return BoundaryConditions( u_t, u_b, u_r, u_l, v_t, v_b, v_r, v_l)
	end

	function create_iu_iv(dims::Dims)::Tuple{Matrix{Int}, Matrix{Int}}
		iu = zeros(dims.nx-1, dims.ny)
		iv = zeros(dims.nx, dims.ny-1)

		k = 1
		# init iu values
		for i = 1:dims.nx -1
			for j = 1:dims.ny
				iu[i,j] = k
				k += 1
			end
		end

		# init iv values
		for i = 1:dims.nx
			for j = 1:dims.ny -1
				iv[i,j] = k
				k += 1
			end
		end

		return iu,iv
	end

	function create_ip(dims::Dims)::Matrix{Int} 
		ip = zeros(dims.nx, dims.ny)

		k = 1
		for i = 1:dims.nx
			for j = 1:dims.ny
				ip[i,j] = k
				k += 1
			end
		end
		
		return ip
	end


	struct SolverResults
		#velocity::Matrix{Float64}
	end

	function execute_solver(dims::Dims, bcs::BoundaryConditions, T::Float64, re::Float64)::SolverResults

		iu, iv = create_iu_iv(dims)
		ip = create_ip(dims)

		cfl_target = 1.0

		dt = calculate_dt(dims, bcs, cfl_target, re)

		n_step = floor(T / dt)

		zero_bcs = zero_boundary_conditions(dims)

		##
		## Define the working vectors
		##

		# velocity matricies
		q_star = zeros(dims.nu)
		q_nm1 = copy(q_star)
		q_n = copy(q_star)
		q_np1 = copy(q_star)

		# pressure matricies
		p_n = zeros(dims.np)
		# P^(n+1)
		p_np1 = copy(p_n)

		first_step_lhs_calculator = FirstStepAx(
			dims, zero_bcs, iu, iv, dt, re
		)

		for step = 1:n_step
			# step variables forward once
			q_nm1 = q_n
			q_n = q_np1
			p_n = p_np1

			#
			# run calculations
			# 
			lap_bc_n = lap(dims, bcs, iu, iv, q_n)
			adv_n = adv(dims, bcs, iu, iv, ip, q_n)
			# TODO: cache this result from the last run - it should be the same
			adv_nm1 = adv(dims, bcs, iu, iv, ip, q_nm1)
			lap_nobc_n = lap(dims, zero_bcs, iu, iv, q_n)

			#
			# assemble things
			#

			# here S =  I + dt nu / 2 L
			# which becomes S = (I + dt L / (2 Re))
			#
			# the whole term is divided by dt
			sq_rhs = 
				# start of S * u_n term
				q_n / dt
				+
				(1 / (2 * re)) * lap_bc_n
				# end of S * u_n term
				# start of advection term
				+
				(3/2) * adv_n
				-
				(1/2) * adv_nm1

			rq_lhs = calculate_ax(first_step_lhs_calculator, q_n)

			# TODO: determine if the laplacian in R^-1 acts 
			# with or without boundary conditions
			# the other L (laplace) operators in R have no boundary conditions
			# but the L (laplace) in S _does_ have boundary conditions
			
			println("calling into conjugate gradient")
			conj_grad(first_step_lhs_calculator, rq_lhs, q_n, sq_rhs, 0.01)

			break
			
		end

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

end

using .CavitySolver
using Printf

function debug_solver(dims::CavitySolver.Dims, bcs::BoundaryConditions)

	iu, iv = create_iu_iv(dims)
	ip = create_ip(dims)

	q = zeros(dims.nu)
	p = ones(dims.np)

	a, b = size(iu)
	c, d = size(iv)
	e = size(q)
	show(a*b + c*d)
	@printf "\n"
	show(e)
	@printf "\n"

	adv(dims, bcs, iu, iv, ip, q)
	#lap(dims, bcs, iu, iv, q)
	#grad(dims, p, iu, iv, ip)

	#div_(dims, bcs, q, iu, iv, ip)
end

const dims = create_dims(1.0, 1.0, 10,10)
bcs = create_boundary_conditions(dims)

debug_solver(dims, bcs)
#execute_solver(dims, bcs, 10., 1000.0)

#debug_conjugate_gradient()
