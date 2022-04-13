module advection
	using Printf

	using .Main.Structs: BoundaryConditions, Dims
	using .Main.indexing: Indexable

	export adv

	# represents all the lambda functions required to calculate 
	# the advection in the X direction
	struct IndexerX
		u_bar_x_left::Function
		u_bar_x_right::Function
		u_bar_y_top::Function
		u_bar_y_bot::Function
		v_bar_x_top::Function
		v_bar_x_bot::Function
	end

	# represents all the lambda functions required to calculate 
	# the advection in the Y direction
	struct IndexerY
		v_bar_x_left::Function
		v_bar_x_right::Function
		v_bar_y_top::Function
		v_bar_y_bot::Function
		u_bar_y_left::Function
		u_bar_y_right::Function
	end

	function adv!(
		dims::Dims, 
		bcs::BoundaryConditions, 
		iu::IU,
		iv::IV,
		ip::IP,
		q::Vector{Float64},
		# dims.nu long
		nq::Vector{Float64}
	) where IU <: Indexable where IV <: Indexable where IP <: Indexable
		#nq = zeros(dims.nu)
		nx = dims.nx
		ny = dims.ny

		##
		## Inner Domain
		##

		x_dir_idx = IndexerX(reg_bar_left, reg_bar_right, reg_bar_top, reg_bar_bot, x_dir_v_bar_top, x_dir_v_bar_bot)
		calculate_nq(dims, bcs, iu, iv, q, nq, 2:nx-2, 2:ny-1, x_dir_idx)

		y_dir_idx = IndexerY(reg_bar_left, reg_bar_right, reg_bar_top, reg_bar_bot, y_dir_u_bar_left, y_dir_u_bar_right)
		calculate_nq(dims, bcs, iu, iv, q, nq, 2:nx-1, 2:ny-2, y_dir_idx)

		##
		## Bottom Inner
		##
		
		x_dir_idx = IndexerX(reg_bar_left, reg_bar_right, reg_bar_top, x_dir_u_bar_bot_bc, x_dir_v_bar_top, x_dir_v_bar_bot_bc)
		calculate_nq(dims, bcs, iu, iv, q, nq, 2:nx-2, 1:1, x_dir_idx)

		y_dir_idx = IndexerY(reg_bar_left, reg_bar_right, reg_bar_top, y_dir_v_bar_bot_bc, y_dir_u_bar_left, y_dir_u_bar_right)
		calculate_nq(dims, bcs, iu, iv, q, nq, 2:nx-1, 1:1, y_dir_idx)

		##
		## Bottom Left
		##
		
		x_dir_idx = IndexerX(x_dir_u_bar_left_bc, reg_bar_right, reg_bar_top, x_dir_u_bar_bot_bc, x_dir_v_bar_top, x_dir_v_bar_bot_bc)
		calculate_nq(dims, bcs, iu, iv, q, nq, 1:1, 1:1, x_dir_idx)

		y_dir_idx = IndexerY(y_dir_v_bar_left_bc, reg_bar_right, reg_bar_top, y_dir_v_bar_bot_bc, y_dir_u_bar_left_bc, y_dir_u_bar_right)
		calculate_nq(dims, bcs, iu, iv, q, nq, 1:1, 1:1, y_dir_idx)
		
		##
		## Bottom Right
		##
		
		x_dir_idx = IndexerX(reg_bar_left, x_dir_u_bar_right_bc, reg_bar_top, x_dir_u_bar_bot_bc, x_dir_v_bar_top, x_dir_v_bar_bot_bc)
		calculate_nq(dims, bcs, iu, iv, q, nq, nx-1:nx-1, 1:1, x_dir_idx)

		y_dir_idx = IndexerY(reg_bar_left, y_dir_v_bar_right_bc, reg_bar_top, y_dir_v_bar_bot_bc, y_dir_u_bar_left, y_dir_u_bar_right_bc)
		calculate_nq(dims, bcs, iu, iv, q, nq, nx:nx, 1:1, y_dir_idx)

		##
		## Left Inner
		##
		
		x_dir_idx = IndexerX(x_dir_u_bar_left_bc, reg_bar_right, reg_bar_top, reg_bar_bot, x_dir_v_bar_top, x_dir_v_bar_bot)
		calculate_nq(dims, bcs, iu, iv, q, nq, 1:1, 2:ny-1, x_dir_idx)

		y_dir_idx = IndexerY(y_dir_v_bar_left_bc, reg_bar_right, reg_bar_top, reg_bar_bot, y_dir_u_bar_left_bc, y_dir_u_bar_right)
		calculate_nq(dims, bcs, iu, iv, q, nq, 1:1, 2:ny-2, y_dir_idx)

		##
		## Right Inner
		##
		
		x_dir_idx = IndexerX(reg_bar_left, x_dir_u_bar_right_bc, reg_bar_top, reg_bar_bot, x_dir_v_bar_top, x_dir_v_bar_bot)
		calculate_nq(dims, bcs, iu, iv, q, nq, nx-1:nx-1, 2:ny-1, x_dir_idx)

		y_dir_idx = IndexerY(reg_bar_left, y_dir_v_bar_right_bc, reg_bar_top, reg_bar_bot, y_dir_u_bar_left, y_dir_u_bar_right_bc)
		calculate_nq(dims, bcs, iu, iv, q, nq, nx:nx, 2:ny-2, y_dir_idx)

		##
		## TOP Inner
		##
		
		x_dir_idx = IndexerX(reg_bar_left, reg_bar_right, x_dir_u_bar_top_bc, reg_bar_bot, x_dir_v_bar_top_bc, x_dir_v_bar_bot)
		calculate_nq(dims, bcs, iu, iv, q, nq, 2:nx-2, ny:ny, x_dir_idx)

		y_dir_idx = IndexerY(reg_bar_left, reg_bar_right, y_dir_v_bar_top_bc, reg_bar_bot, y_dir_u_bar_left, y_dir_u_bar_right)
		calculate_nq(dims, bcs, iu, iv, q, nq, 2:nx-1, ny-1:ny-1, y_dir_idx)

		##
		## TOP Left
		##

		x_dir_idx = IndexerX(x_dir_u_bar_left_bc, reg_bar_right, x_dir_u_bar_top_bc, reg_bar_bot, x_dir_v_bar_top_bc, x_dir_v_bar_bot)
		calculate_nq(dims, bcs, iu, iv, q, nq, 1:1, ny:ny, x_dir_idx)

		y_dir_idx = IndexerY(y_dir_v_bar_left_bc, reg_bar_right, y_dir_v_bar_top_bc, reg_bar_bot, y_dir_u_bar_left_bc, y_dir_u_bar_right)
		calculate_nq(dims, bcs, iu, iv, q, nq, 1:1, ny-1:ny-1, y_dir_idx)

		##
		## TOP Right
		##

		x_dir_idx = IndexerX(reg_bar_left, x_dir_u_bar_right_bc, x_dir_u_bar_top_bc, reg_bar_bot, x_dir_v_bar_top_bc, x_dir_v_bar_bot)
		calculate_nq(dims, bcs, iu, iv, q, nq, nx-1:nx-1, ny:ny, x_dir_idx)

		y_dir_idx = IndexerY(reg_bar_left, y_dir_v_bar_right_bc, y_dir_v_bar_top_bc, reg_bar_bot, y_dir_u_bar_left, y_dir_u_bar_right_bc)
		calculate_nq(dims, bcs, iu, iv, q, nq, nx:nx, ny-1:ny-1, y_dir_idx)
	end

	# calculate entire nq for loop for X direction values
	@inline function calculate_nq(
		dims::Dims,
		bcs::BoundaryConditions,
		iu::IU,
		iv::IV,
		q::Vector{Float64},
		nq::Vector{Float64},
		x_range::UnitRange{Int},
		y_range::UnitRange{Int},
		indexer::IndexerX
	) where IU <: Indexable where IV <: Indexable

		for i = x_range
			for j = y_range
				u_bar_x_left = indexer.u_bar_x_left(q, iu, bcs, i,j)
				u_bar_x_right = indexer.u_bar_x_right(q, iu, bcs, i,j)

				u_bar_y_top = indexer.u_bar_y_top(q, iu, bcs, i,j)
				u_bar_y_bot = indexer.u_bar_y_bot(q, iu, bcs, i,j)

				v_bar_x_top = indexer.v_bar_x_top(q, iv, bcs, i, j)
				v_bar_x_bot = indexer.v_bar_x_bot(q, iv, bcs, i, j)

				# calculate the final value
				@inbounds nq[iu[i,j]] = (u_bar_x_right^2 - u_bar_x_left^2)/dims.dx + 
					(u_bar_y_top * v_bar_x_top - u_bar_y_bot * v_bar_x_bot) / dims.dy
			end
		end
	end

	# calculate entire nq for loop for y direction values
	@inline function calculate_nq(
		dims::Dims,
		bcs::BoundaryConditions,
		iu::IU,
		iv::IV,
		q::Vector{Float64},
		nq::Vector{Float64},
		x_range::UnitRange{Int},
		y_range::UnitRange{Int},
		indexer::IndexerY
	) where IU <: Indexable where IV <: Indexable
		for i = x_range
			for j = y_range
				v_bar_x_left = indexer.v_bar_x_left(q, iv, bcs, i,j)
				v_bar_x_right = indexer.v_bar_x_right(q, iv, bcs, i,j)

				v_bar_y_top = indexer.v_bar_y_top(q, iv, bcs, i,j)
				v_bar_y_bot = indexer.v_bar_y_bot(q, iv, bcs, i,j)

				u_bar_y_left = indexer.u_bar_y_left(q, iu, bcs, i, j)
				u_bar_y_right = indexer.u_bar_y_right(q, iu, bcs, i, j)

				@inbounds nq[iv[i,j]] = (u_bar_y_right * v_bar_x_right - u_bar_y_left * v_bar_x_left) / dims.dx + 
					(v_bar_y_top^2 - v_bar_y_bot^2) / dims.dy
			end
		end
	end

	######
	###### GENERAL INDEXING - can be used for either x or y directions as long as you are not on a boundary
	######

	reg_bar_left(q::Vector{Float64}, indexer, bcs::BoundaryConditions, i::Int, j::Int)::Float64 =
		@inbounds (1/2) * (q[indexer[i-1,j]] + q[indexer[i,j]])

	reg_bar_right(q::Vector{Float64}, indexer, bcs::BoundaryConditions, i::Int, j::Int)::Float64 =
		@inbounds (1/2) * (q[indexer[i+1,j]] + q[indexer[i,j]])

	reg_bar_top(q::Vector{Float64}, indexer, bcs::BoundaryConditions, i::Int, j::Int)::Float64 =
		@inbounds (1/2) * (q[indexer[i,j+1]] + q[indexer[i,j]])
	reg_bar_bot(q::Vector{Float64}, indexer, bcs::BoundaryConditions, i::Int, j::Int)::Float64 =
		@inbounds (1/2) * (q[indexer[i,j-1]] + q[indexer[i,j]])

	######
	###### X DIRECTION INDEXING
	######

	#
	# inner domain indexing for v values above and below the current x value
	# 
	
	x_dir_v_bar_top(q::Vector{Float64}, indexer, bcs::BoundaryConditions, i::Int, j::Int)::Float64 =
		# dont need j+1 here since values at y=j are already above the current i value
		@inbounds (1/2) * (q[indexer[i+1,j]] + q[indexer[i,j]])

	x_dir_v_bar_bot(q::Vector{Float64}, indexer, bcs::BoundaryConditions, i::Int, j::Int)::Float64 =
		@inbounds (1/2) * (q[indexer[i+1,j-1]] + q[indexer[i,j-1]])

	#
	# boundary stuff top / bottom
	#

	# top of boundary indexing for y
	x_dir_v_bar_top_bc(q::Vector{Float64}, indexer, bcs::BoundaryConditions, i::Int, j::Int)::Float64 =
		@inbounds (1/2) * (bcs.v_t[i] + bcs.v_t[i+1])

	# bottom boundary indexing for y
	x_dir_v_bar_bot_bc(q::Vector{Float64}, indexer, bcs::BoundaryConditions, i::Int, j::Int)::Float64 =
		@inbounds (1/2) * (bcs.v_b[i] + bcs.v_b[i+1])

	
	# bot average with the ghost node is just the boundary condition
	x_dir_u_bar_bot_bc(q::Vector{Float64}, indexer, bcs::BoundaryConditions, i::Int, j::Int)::Float64 =
		@inbounds bcs.u_b[i]

	# top average with the ghost node is just the boundary condition
	x_dir_u_bar_top_bc(q::Vector{Float64}, indexer, bcs::BoundaryConditions, i::Int, j::Int)::Float64 =
		@inbounds bcs.u_t[i]

	#
	# boundary stuff left / right
	#

	x_dir_u_bar_left_bc(q::Vector{Float64}, indexer, bcs::BoundaryConditions, i::Int, j::Int)::Float64 =
		@inbounds (1/2) * (q[indexer[i,j]] + bcs.u_l[j])

	x_dir_u_bar_right_bc(q::Vector{Float64}, indexer, bcs::BoundaryConditions, i::Int, j::Int)::Float64 =
		@inbounds (1/2) * (q[indexer[i,j]] + bcs.u_r[j])
	
	######
	###### Y DIRECTION INDEXING
	######
	
	#
	# indexing x velocity points on the left and right sides of the current y velocity point
	# and averaging them to the level at which the y velocity is acting
	#

	y_dir_u_bar_left(q::Vector{Float64}, indexer, bcs::BoundaryConditions, i::Int, j::Int)::Float64 =
		@inbounds (1/2) * (q[indexer[i-1,j]] + q[indexer[i-1,j+1]])

	y_dir_u_bar_right(q::Vector{Float64}, indexer, bcs::BoundaryConditions, i::Int, j::Int)::Float64 =
		# values with x=i are already on the right side of our node
		@inbounds (1/2) * (q[indexer[i,j]] + q[indexer[i,j+1]])

	#
	# top and bottom boundary conditions for y velocity
	#

	y_dir_v_bar_bot_bc(q::Vector{Float64}, indexer, bcs::BoundaryConditions, i::Int, j::Int)::Float64 =
		@inbounds (1/2) * (bcs.v_b[i] + q[indexer[i,j]])

	y_dir_v_bar_top_bc(q::Vector{Float64}, indexer, bcs::BoundaryConditions, i::Int, j::Int)::Float64 =
		@inbounds (1/2) * (bcs.v_t[i] + q[indexer[i,j]])

	#
	# left and right conditions for y velocity
	#
	
	y_dir_v_bar_left_bc(q::Vector{Float64}, indexer, bcs::BoundaryConditions, i::Int, j::Int)::Float64 =
		@inbounds bcs.v_l[j]

	y_dir_v_bar_right_bc(q::Vector{Float64}, indexer, bcs::BoundaryConditions, i::Int, j::Int)::Float64 =
		@inbounds bcs.v_r[j]

	y_dir_u_bar_left_bc(q::Vector{Float64}, indexer, bcs::BoundaryConditions, i::Int, j::Int)::Float64 =
		@inbounds (1/2) * (bcs.u_l[j] + bcs.u_l[j+1])

	y_dir_u_bar_right_bc(q::Vector{Float64}, indexer, bcs::BoundaryConditions, i::Int, j::Int)::Float64 =
		@inbounds (1/2) * (bcs.u_r[j] + bcs.u_r[j+1])
end
