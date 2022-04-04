include("structs.jl")
include("adv.jl")
include("lap.jl")
include("gradient.jl")
include("divergence.jl")

module CavitySolver
	using Printf

	using .Main.Structs:BoundaryConditions
	using .Main.Structs:Dims
	using .Main.advection: adv
	using .Main.laplacian: lap
	using .Main.gradient: grad
	using .Main.divergence: div_

	export create_dims, create_iu_iv, test, create_boundary_conditions, create_ip
	export adv, lap, grad, div_, BoundaryConditions, Dims

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

		u_t = ones(nx-1, 1);
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

	function test(dims::Dims)
		@printf "called2"
	end
end

using .CavitySolver
#using .Main.laplacian: lap
#using .Structs
using Printf

const dims = create_dims(1.0, 1.0, 10,10)
bcs = create_boundary_conditions(dims)

iu, iv = create_iu_iv(dims)
ip = create_ip(dims)

q = zeros(dims.nu, 1)
p = ones(dims.np, 1)

a, b = size(iu)
c, d = size(iv)
e, f = size(q)
show(a*b + c*d)
@printf "\n"
show(e*f)
@printf "\n"

#lap(dims, bcs, iu, iv, q)
#adv(dims, bcs, iu, iv, ip, q)
#grad(dims, p, iu, iv, ip)

div_(dims, bcs, q, iu, iv, ip)
