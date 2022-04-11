module laplacian
	using Printf
	using .Main.Structs: BoundaryConditions, Dims
	using .Main.indexing: Indexable

	export lap

	struct Indexer
		left::Function
		right::Function
		top::Function
		bottom::Function
		center::Function
	end

	function lap(
		dims::Dims, 
		bcs::BoundaryConditions, 
		iu::IU,
		iv::IV,
		q::Vector{Float64},
		# dims.nu points long
		Lq::Vector{Float64}
	)::Vector{Float64} where IV <: Indexable where IU <: Indexable 

		#Lq = zeros(dims.nu)

		nx = dims.nx
		ny = dims.ny

		####
		#### Inner domain
		####

		indexer = Indexer(reg_left, reg_right, reg_top, reg_bot, reg_center)

		# u velocity points
		calc_lap(dims, Lq, bcs, 2:nx-2, 2:ny-1, q, iu, indexer)
		# v velocity points
		calc_lap(dims, Lq, bcs, 2:nx-1, 2:ny-2, q, iv, indexer)

		####
		#### Bottom Inner
		####

		# u velocity points
		indexer = Indexer(reg_left, reg_right, reg_top, u_bottom_bc, reg_center)
		calc_lap(dims, Lq, bcs, 2:nx-2, 1:1, q, iu, indexer)
		# v velocity points
		indexer = Indexer(reg_left, reg_right, reg_top, v_bottom_bc, reg_center)
		calc_lap(dims, Lq, bcs, 2:nx-1, 1:1, q, iv, indexer)

		####
		#### Bottom Left
		####

		# u velocity points
		indexer = Indexer(u_left_bc, reg_right, reg_top, u_bottom_bc, reg_center)
		calc_lap(dims, Lq, bcs, 1:1, 1:1, q, iu, indexer)
		# v velocity points
		indexer = Indexer(v_left_bc, reg_right, reg_top, v_bottom_bc, reg_center)
		calc_lap(dims, Lq, bcs, 1:1, 1:1, q, iv, indexer)

		####
		#### Bottom Right
		####

		# u velocity points
		indexer = Indexer(reg_left, u_right_bc, reg_top, u_bottom_bc, reg_center)
		calc_lap(dims, Lq, bcs, nx-1:nx-1, 1:1, q, iu, indexer)
		# v velocity points
		indexer = Indexer(reg_left, v_right_bc, reg_top, v_bottom_bc, reg_center)
		calc_lap(dims, Lq, bcs, nx:nx, 1:1, q, iv, indexer)

		####
		#### Left Inner
		####

		# u velocity points
		indexer = Indexer(u_left_bc, reg_right, reg_top, reg_bot, reg_center)
		calc_lap(dims, Lq, bcs, 1:1, 2:ny-1, q, iu, indexer)
		# v velocity points
		indexer = Indexer(v_left_bc, reg_right, reg_top, reg_bot, reg_center)
		calc_lap(dims, Lq, bcs, 1:1, 2:ny-2, q, iv, indexer)

		####
		#### Right Inner
		####

		# u velocity points
		indexer = Indexer(reg_left, u_right_bc, reg_top, reg_bot, reg_center)
		calc_lap(dims, Lq, bcs, nx-1:nx-1, 2:ny-1, q, iu, indexer)
		# v velocity points
		indexer = Indexer(reg_left, v_right_bc, reg_top, reg_bot, reg_center)
		calc_lap(dims, Lq, bcs, nx:nx, 2:ny-2, q, iv, indexer)

		#println("\nlaplacian before top boundary conditions\n")
		#display(transpose(Lq))

		####
		#### Top Inner
		####

		# u velocity points
		indexer = Indexer(reg_left, reg_right, u_top_bc, reg_bot, reg_center)
		calc_lap(dims, Lq, bcs, 2:nx-2, ny:ny, q, iu, indexer)
		# v velocity points
		indexer = Indexer(reg_left, reg_right, v_top_bc, reg_bot, reg_center)
		calc_lap(dims, Lq, bcs, 2:nx-1, ny-1:ny-1, q, iv, indexer)

		####
		#### Top Left
		####

		# u velocity points
		indexer = Indexer(u_left_bc, reg_right, u_top_bc, reg_bot, reg_center)
		calc_lap(dims, Lq, bcs, 1:1, ny:ny, q, iu, indexer)
		# v velocity points
		indexer = Indexer(v_left_bc, reg_right, v_top_bc, reg_bot, reg_center)
		calc_lap(dims, Lq, bcs, 1:1, ny-1:ny-1, q, iv, indexer)

		####
		#### Top Right
		####

		# u velocity points
		indexer = Indexer(reg_left, u_right_bc, u_top_bc, reg_bot, reg_center)
		calc_lap(dims, Lq, bcs, nx-1:nx-1, ny:ny, q, iu, indexer)
		# v velocity points
		indexer = Indexer(reg_left, v_right_bc, v_top_bc, reg_bot, reg_center)
		calc_lap(dims, Lq, bcs, nx:nx, ny-1:ny-1, q, iv, indexer)

		#println("\nlaplacian after top boundary conditions\n")
		#display(transpose(Lq))

		####
		#### Return
		####

		return Lq
	end

	@inline function calc_lap(
		dims::Dims,
		Lq::Vector{Float64},
		bcs::BoundaryConditions,
		x_range::UnitRange{Int64},
		y_range::UnitRange{Int64},
		q::Vector{Float64},
		points::T,
		indexer::Indexer,
	) where T <: Indexable
		for i = x_range
			for j = y_range
				left = indexer.left(q, bcs, points, i, j)
				right = indexer.right(q, bcs, points, i, j)

				top = indexer.top(q, bcs, points, i, j)
				bottom = indexer.bottom(q, bcs, points, i, j)

				center = indexer.center(q, bcs, points, i, j)

				d_dx2 = (right - 2*center + left) / dims.dx^2
				d_dy2 = (top - 2*center + bottom) / dims.dy^2

				Lq[points[i,j]] = d_dx2 + d_dy2
			end
		end
	end

	function reg_left(q::Vector{Float64}, bcs::BoundaryConditions, indexer::T, i::Int64, j::Int64)::Float64 where T <: Indexable
		return q[indexer[i-1, j]]
	end

	function reg_right(q::Vector{Float64}, bcs::BoundaryConditions, indexer::T, i::Int64, j::Int64)::Float64 where T <: Indexable
		return q[indexer[i+1, j]]
	end

	function reg_top(q::Vector{Float64}, bcs::BoundaryConditions, indexer::T, i::Int64, j::Int64)::Float64 where T <: Indexable 
		return q[indexer[i, j+1]]
	end

	function reg_bot(q::Vector{Float64}, bcs::BoundaryConditions, indexer::T, i::Int64, j::Int64)::Float64 where T <: Indexable 
		return q[indexer[i, j-1]]
	end
	function reg_center(q::Vector{Float64}, bcs::BoundaryConditions, indexer::T, i::Int64, j::Int64)::Float64 where T <: Indexable 
		return q[indexer[i, j]]
	end

	#####
	##### u velocity specific boundaries
	#####

	function u_bottom_bc(q::Vector{Float64}, bcs::BoundaryConditions, indexer::T, i::Int64, j::Int64)::Float64 where T <: Indexable
		return 2 * bcs.u_b[i] - q[indexer[i,j]]
	end

	function u_top_bc(q::Vector{Float64}, bcs::BoundaryConditions, indexer::T, i::Int64, j::Int64)::Float64 where T <: Indexable
		return 2 * bcs.u_t[i] - q[indexer[i,j]]
	end

	function u_right_bc(q::Vector{Float64}, bcs::BoundaryConditions, indexer::T, i::Int64, j::Int64)::Float64 where T <: Indexable
		return bcs.u_r[j]
	end

	function u_left_bc(q::Vector{Float64}, bcs::BoundaryConditions, indexer::T, i::Int64, j::Int64)::Float64 where T <: Indexable
		return bcs.u_l[j]
	end

	#####
	##### v velocity specific boundaries
	#####

	function v_right_bc(q::Vector{Float64}, bcs::BoundaryConditions, indexer::T, i::Int64, j::Int64)::Float64 where T <: Indexable
		return 2 * bcs.v_r[j] - q[indexer[i,j]]
	end

	function v_left_bc(q::Vector{Float64}, bcs::BoundaryConditions, indexer::T, i::Int64, j::Int64)::Float64 where T <: Indexable
		return 2 * bcs.v_l[j] - q[indexer[i,j]]
	end

	function v_bottom_bc(q::Vector{Float64}, bcs::BoundaryConditions, indexer::T, i::Int64, j::Int64)::Float64 where T <: Indexable
		return bcs.v_b[i]
	end

	function v_top_bc(q::Vector{Float64}, bcs::BoundaryConditions, indexer::T, i::Int64, j::Int64)::Float64 where T <: Indexable
		return bcs.v_t[i]
	end

end
