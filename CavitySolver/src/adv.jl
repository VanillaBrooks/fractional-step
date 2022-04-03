module advection
	using Printf

	using .Main.Structs: BoundaryConditions, Dims

	export adv

	struct IndexerX
		u_bar_x_left::Function
		u_bar_x_right::Function
		u_bar_y_top::Function
		u_bar_y_bot::Function
		v_bar_x_top::Function
		v_bar_x_bot::Function
	end

	struct IndexerY
		v_bar_x_left::Function
		v_bar_x_right::Function
		v_bar_y_top::Function
		v_bar_y_bot::Function
		u_bar_y_left::Function
		u_bar_y_right::Function
	end

	function adv(dims::Dims, bcs::BoundaryConditions, iu::Matrix{Int64}, iv::Matrix{Int64}, ip::Matrix{Int64}, q::Matrix{Float64})::Matrix{Float64}
		nq = zeros(dims.nu, 1)
		nx = dims.nx
		ny = dims.ny

		##
		## Inner Domain
		##

		x_dir_idx = IndexerX(reg_bar_left, reg_bar_right, reg_bar_top, reg_bar_bot, x_dir_v_bar_top, x_dir_v_bar_bot)
		calculate_nq(dims, bcs, iu, iv, q, nq, 2:nx-2, 2:ny-1, x_dir_idx)

		y_dir_idx = IndexerY(reg_bar_left, reg_bar_right, reg_bar_top, reg_bar_bot, y_dir_u_bar_left, y_dir_u_bar_right)
		calculate_nq(dims, bcs, iu, iv, q, nq, 2:nx-1, 2:ny-2, y_dir_idx)

		# ##
		# ## Bottom Inner
		# ##
		# 
		# # x dir
		# x_indexers = Indexer(reg_plus_i, reg_minus_i, reg_plus_j, u_bot_minus_j);
		# y_indexers = Indexer(reg_plus_i, reg_minus_i, reg_plus_j, v_bot_minus_j);
		# calculate_nq_x_dir(dims, bcs, iu, iv, q, nq, 2:nx-2, 1:1, x_indexers, y_indexers)
		# # y dir
		# calculate_nq_y_dir(dims, bcs, iu, iv, q, nq, 2:nx-1, 1:1, x_indexers, y_indexers)

		# ##
		# ## Bottom Left
		# ##
		# 
		# x_indexers = Indexer(reg_plus_i, u_left_minus_i, reg_plus_j, u_bot_minus_j);
		# y_indexers = Indexer(reg_plus_i, v_left_minus_i, reg_plus_j, v_bot_minus_j);
		# calculate_nq_x_dir(dims, bcs, iu, iv, q, nq, 1:1, 1:1, x_indexers, y_indexers)
		# # y dir
		# calculate_nq_y_dir(dims, bcs, iu, iv, q, nq, 1:1, 1:1, x_indexers, y_indexers)

		# ##
		# ## Bottom Right
		# ##
		# 
		# x_indexers = Indexer(u_right_plus_i, reg_minus_i, reg_plus_j, u_bot_minus_j);
		# y_indexers = Indexer(v_right_plus_i, reg_minus_i, reg_plus_j, v_bot_minus_j);
		# calculate_nq_x_dir(dims, bcs, iu, iv, q, nq, nx-1:nx-1, 1:1, x_indexers, y_indexers)
		# # y dir
		# # TODO: outside the boundary here if we index nx:nx
		# calculate_nq_y_dir(dims, bcs, iu, iv, q, nq, nx:nx, 1:1, x_indexers, y_indexers) 


		# ##
		# ## Left Inner
		# ##
		# 
		# # x dir
		# x_indexers = Indexer(reg_plus_i, u_left_minus_i, reg_plus_j, reg_minus_j);
		# y_indexers = Indexer(reg_plus_i, v_left_minus_i, reg_plus_j, reg_minus_j);
		# calculate_nq_x_dir(dims, bcs, iu, iv, q, nq, 1:1, 2:ny-1, x_indexers, y_indexers)
		# # y dir
		# calculate_nq_y_dir(dims, bcs, iu, iv, q, nq, 1:1, 2:ny-2, x_indexers, y_indexers)

		# ##
		# ## Right Inner
		# ##
		# 
		# # x dir
		# x_indexers = Indexer(u_right_plus_i, reg_minus_i, reg_plus_j, reg_minus_j);
		# y_indexers = Indexer(v_right_plus_i, reg_minus_i, reg_plus_j, reg_minus_j);
		# calculate_nq_x_dir(dims, bcs, iu, iv, q, nq, nx-1:nx-1, 2:ny-1, x_indexers, y_indexers)
		# # y dir
		# # TODO: outside the boundary here if we index nx:nx which is required for the v values
		# calculate_nq_y_dir(dims, bcs, iu, iv, q, nq, nx-1:nx-1, 2:ny-2, x_indexers, y_indexers)

		# ##
		# ## TOP Inner
		# ##
		# # x dir
		# x_indexers = Indexer(reg_plus_i, reg_minus_i, u_top_plus_j, reg_minus_j);
		# y_indexers = Indexer(reg_plus_i, reg_minus_i, v_top_plus_j, reg_minus_j);
		# # TODO: this bound needs to be 
		# calculate_nq_x_dir(dims, bcs, iu, iv, q, nq, 2:nx-1, ny:ny, x_indexers, y_indexers)
		# # y dir
		# # TODO: outside the boundary here if we index nx:nx which is required for the v values
		# calculate_nq_y_dir(dims, bcs, iu, iv, q, nq, 2:nx-1, ny-1:ny-1, x_indexers, y_indexers)

		## RETURN

		return nq
	end

	# calculate entire nq for loop
	function calculate_nq(
		dims::Dims,
		bcs::BoundaryConditions,
		iu::Matrix{Int64}, 
		iv::Matrix{Int64}, 
		q::Matrix{Float64},
		nq::Matrix{Float64},
		x_range::UnitRange{Int},
		y_range::UnitRange{Int},
		indexer::IndexerX
	)

		for i = x_range
			for j = y_range
				u_bar_x_left = indexer.u_bar_x_left(q, iu, bcs, i,j)
				u_bar_x_right = indexer.u_bar_x_right(q, iu, bcs, i,j)

				u_bar_y_top = indexer.u_bar_y_top(q, iu, bcs, i,j)
				u_bar_y_bot = indexer.u_bar_y_bot(q, iu, bcs, i,j)

				v_bar_x_top = indexer.v_bar_x_top(q, iv, bcs, i, j)
				v_bar_x_bot = indexer.v_bar_x_bot(q, iv, bcs, i, j)

				# calculate the final value
				nq[iu[i,j]] = (u_bar_x_right^2 - u_bar_x_left^2)/dims.dx + 
					(u_bar_y_top * v_bar_x_top - u_bar_y_bot * v_bar_x_bot) / dims.dy
			end
		end
		
	end

	function calculate_nq(
		dims::Dims,
		bcs::BoundaryConditions,
		iu::Matrix{Int64}, 
		iv::Matrix{Int64}, 
		q::Matrix{Float64},
		nq::Matrix{Float64},
		x_range::UnitRange{Int},
		y_range::UnitRange{Int},
		indexer::IndexerY
	)

		for i = x_range
			for j = y_range
				v_bar_x_left = indexer.v_bar_x_left(q, iv, bcs, i,j)
				v_bar_x_right = indexer.v_bar_x_right(q, iv, bcs, i,j)

				v_bar_y_top = indexer.v_bar_y_top(q, iv, bcs, i,j)
				v_bar_y_bot = indexer.v_bar_y_bot(q, iv, bcs, i,j)

				u_bar_y_left = indexer.u_bar_y_left(q, iu, bcs, i, j)
				u_bar_y_right = indexer.u_bar_y_right(q, iu, bcs, i, j)

				nq[iv[i,j]] = (u_bar_y_right * v_bar_x_right - u_bar_y_left * v_bar_x_left) / dims.dx + 
					(v_bar_y_top^2 - v_bar_y_bot^2) / dims.dy
			end
		end
	end

	######
	###### GENERAL INDEXING - can be used for either x or y directions
	######

	reg_bar_left(q::Matrix{Float64}, indexer::Matrix{Int64}, bcs::BoundaryConditions, i::Int, j::Int) =
		(1/2) * (q[indexer[i-1,j]] + q[indexer[i,j]])
	reg_bar_right(q::Matrix{Float64}, indexer::Matrix{Int64}, bcs::BoundaryConditions, i::Int, j::Int) =
		(1/2) * (q[indexer[i+1,j]] + q[indexer[i,j]])

	reg_bar_top(q::Matrix{Float64}, indexer::Matrix{Int64}, bcs::BoundaryConditions, i::Int, j::Int) =
		(1/2) * (q[indexer[i,j+1]] + q[indexer[i,j]])
	reg_bar_bot(q::Matrix{Float64}, indexer::Matrix{Int64}, bcs::BoundaryConditions, i::Int, j::Int) =
		(1/2) * (q[indexer[i,j-1]] + q[indexer[i,j]])

	######
	###### X DIRECTION INDEXING
	######

	#
	# inner domain indexing for v values above and below the current x value
	# 
	
	x_dir_v_bar_top(q::Matrix{Float64}, indexer::Matrix{Int64}, bcs::BoundaryConditions, i::Int, j::Int) =
		# dont need j+1 here since values at y=j are already above the current i value
		(1/2) * (q[indexer[i+1,j]] + q[indexer[i,j]])

	x_dir_v_bar_bot(q::Matrix{Float64}, indexer::Matrix{Int64}, bcs::BoundaryConditions, i::Int, j::Int) =
		(1/2) * (q[indexer[i+1,j-1]] + q[indexer[i,j-1]])

	#
	# boundary stuff
	#

	# top of boundary indexing for y
	x_dir_v_plus_j_bc(q::Matrix{Float64}, indexer::Matrix{Int64}, bcs::BoundaryConditions, i::Int, j::Int) =
		(1/2) * (bcs.v_t[i] + bcs.v_t[i+1])

	# bottom boundary indexing for y
	x_dir_v_minus_j_bc(q::Matrix{Float64}, indexer::Matrix{Int64}, bcs::BoundaryConditions, i::Int, j::Int) =
		(1/2) * (bcs.v_b[i] + bcs.v_b[i+1])

	#
	# left and right boundary conditions for x velocity
	#
	
	# right and left conditions for u
	u_left_minus_i(q::Matrix{Float64}, indexer::Matrix{Int64}, bcs::BoundaryConditions, i::Int, j::Int) =
		(1/2) * (bcs.u_l[j] + q[indexer[i,j]])

	u_right_plus_i(q::Matrix{Float64}, indexer::Matrix{Int64}, bcs::BoundaryConditions, i::Int, j::Int) =
		(1/2) * (bcs.u_r[j] + q[indexer[i,j]])
	
	######
	###### Y DIRECTION INDEXING
	######
	
	#
	# indexing x velocity points on the left and right sides of the current y velocity point
	# and averaging them to the level at which the y velocity is acting
	#

	y_dir_u_bar_left(q::Matrix{Float64}, indexer::Matrix{Int64}, bcs::BoundaryConditions, i::Int, j::Int) =
		(1/2) * (q[indexer[i-1,j]] + q[indexer[i-1,j+1]])

	y_dir_u_bar_right(q::Matrix{Float64}, indexer::Matrix{Int64}, bcs::BoundaryConditions, i::Int, j::Int) =
		# values with x=i are already on the right side of our node
		(1/2) * (q[indexer[i,j]] + q[indexer[i,j+1]])

	#
	# top and bottom boundary conditions for y velocity
	#

	v_bot_minus_j(q::Matrix{Float64}, indexer::Matrix{Int64}, bcs::BoundaryConditions, i::Int, j::Int) =
		(1/2) * (bcs.v_b[i] + q[indexer[i,j]])
	v_top_plus_j(q::Matrix{Float64}, indexer::Matrix{Int64}, bcs::BoundaryConditions, i::Int, j::Int) =
		(1/2) * (bcs.v_t[i] + q[indexer[i,j]])

	##
	## TOP AND BOTTOM SPECIAL BOUDNARIES 
	##
	
	# # top and bottom BCs for u
	# # the boundary value is the average with the ghost node
	# u_bot_minus_j(q::Matrix{Float64}, indexer::Matrix{Int64}, bcs::BoundaryConditions, i::Int, j::Int) =
	# 	bcs.u_b[i]

	# u_top_plus_j(q::Matrix{Float64}, indexer::Matrix{Int64}, bcs::BoundaryConditions, i::Int, j::Int) =
	#	bcs.u_t[i]
		

	##
	## LEFT AND RIGHT SPECIAL BOUDNARIES 
	##

	# left and right conditions for v
	v_left_minus_i(q::Matrix{Float64}, indexer::Matrix{Int64}, bcs::BoundaryConditions, i::Int, j::Int) =
		bcs.v_l[j]
	v_right_plus_i(q::Matrix{Float64}, indexer::Matrix{Int64}, bcs::BoundaryConditions, i::Int, j::Int) =
		bcs.v_r[j]


end
