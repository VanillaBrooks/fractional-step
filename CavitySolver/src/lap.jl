module laplacian
	using Printf
	using .Main.Structs: BoundaryConditions, Dims

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
		iu::Matrix{Int64},
		iv::Matrix{Int64},
		q::Matrix{Float64}
	)::Matrix{Float64}

		Lq = zeros((dims.nu, 1))

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
		#### Top Inner
		####

		# u velocity points
		indexer = Indexer(reg_left, u_right_bc, u_top_bc, reg_bot, reg_center)
		calc_lap(dims, Lq, bcs, nx-1:nx-1, ny:ny, q, iu, indexer)
		# v velocity points
		indexer = Indexer(reg_left, v_right_bc, v_top_bc, reg_bot, reg_center)
		calc_lap(dims, Lq, bcs, nx:nx, ny-1:ny-1, q, iv, indexer)

		####
		#### Return
		####

		return Lq
	end


	function calc_lap(
		dims::Dims,
		Lq::Matrix{Float64},
		bcs::BoundaryConditions,
		x_range::UnitRange{Int64},
		y_range::UnitRange{Int64},
		q::Matrix{Float64},
		points::Matrix{Int64},
		indexer::Indexer,
	)
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

	reg_left(q::Matrix{Float64}, bcs::BoundaryConditions, indexer::Matrix{Int64}, i::Int64, j::Int64) = 
		q[indexer[i-1, j]]
	reg_right(q::Matrix{Float64}, bcs::BoundaryConditions, indexer::Matrix{Int64}, i::Int64, j::Int64) = 
		q[indexer[i, j]]
	reg_top(q::Matrix{Float64}, bcs::BoundaryConditions, indexer::Matrix{Int64}, i::Int64, j::Int64) = 
		q[indexer[i, j+1]]
	reg_bot(q::Matrix{Float64}, bcs::BoundaryConditions, indexer::Matrix{Int64}, i::Int64, j::Int64) = 
		q[indexer[i, j-1]]
	reg_center(q::Matrix{Float64}, bcs::BoundaryConditions, indexer::Matrix{Int64}, i::Int64, j::Int64) = 
		q[indexer[i, j]]

	#####
	##### u velocity specific boundaries
	#####

	u_bottom_bc(q::Matrix{Float64}, bcs::BoundaryConditions, indexer::Matrix{Int64}, i::Int64, j::Int64) = 
		2 * bcs.u_b[i] - q[indexer[i,j]]

	u_top_bc(q::Matrix{Float64}, bcs::BoundaryConditions, indexer::Matrix{Int64}, i::Int64, j::Int64) = 
		2 * bcs.u_t[i] - q[indexer[i,j]]

	u_right_bc(q::Matrix{Float64}, bcs::BoundaryConditions, indexer::Matrix{Int64}, i::Int64, j::Int64) = 
		bcs.u_r[j]

	u_left_bc(q::Matrix{Float64}, bcs::BoundaryConditions, indexer::Matrix{Int64}, i::Int64, j::Int64) = 
		bcs.u_l[j]

	#####
	##### v velocity specific boundaries
	#####

	v_right_bc(q::Matrix{Float64}, bcs::BoundaryConditions, indexer::Matrix{Int64}, i::Int64, j::Int64) = 
		2 * bcs.v_r[j] - q[indexer[i,j]]

	v_left_bc(q::Matrix{Float64}, bcs::BoundaryConditions, indexer::Matrix{Int64}, i::Int64, j::Int64) = 
		2 * bcs.v_l[j] - q[indexer[i,j]]

	v_bottom_bc(q::Matrix{Float64}, bcs::BoundaryConditions, indexer::Matrix{Int64}, i::Int64, j::Int64) = 
		bcs.v_b[i]

	v_top_bc(q::Matrix{Float64}, bcs::BoundaryConditions, indexer::Matrix{Int64}, i::Int64, j::Int64) = 
		bcs.v_t[i]

end
