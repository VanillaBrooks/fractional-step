module advection
	using Printf

	using .Main.Structs: BoundaryConditions, Dims

	export adv

	struct Indexer
		pi::Function
		mi::Function
		pj::Function
		mj::Function
	end

	function adv(dims::Dims, bcs::BoundaryConditions, iu::Matrix{Int64}, iv::Matrix{Int64}, ip::Matrix{Int64}, q::Matrix{Float64})::Matrix{Float64}
		nq = zeros(dims.nu, 1)
		nx = dims.nx
		ny = dims.ny

		##
		## Inner Domain
		##

		x_indexers = Indexer(reg_plus_i, reg_minus_i, reg_plus_j, reg_minus_j);
		y_indexers = Indexer(reg_plus_i, reg_minus_i, reg_plus_j, reg_minus_j);
		calculate_nq_x_dir(dims, bcs, iu, iv, q, nq, 2:nx-2, 2:ny-1,x_indexers, y_indexers)

		calc = final_calc_v
		calculate_nq_y_dir(dims, bcs, iu, iv, q, nq, 2:nx-1, 2:ny-2,x_indexers, y_indexers)

		##
		## Bottom Inner
		##
		
		# x dir
		x_indexers = Indexer(reg_plus_i, reg_minus_i, reg_plus_j, u_bot_minus_j);
		y_indexers = Indexer(reg_plus_i, reg_minus_i, reg_plus_j, v_bot_minus_j);
		calculate_nq_x_dir(dims, bcs, iu, iv, q, nq, 2:nx-2, 1:1, x_indexers, y_indexers)
		# y dir
		calculate_nq_y_dir(dims, bcs, iu, iv, q, nq, 2:nx-1, 1:1, x_indexers, y_indexers)

		##
		## Bottom Left
		##
		
		x_indexers = Indexer(reg_plus_i, reg_minus_i, reg_plus_j, u_bot_minus_j);
		y_indexers = Indexer(reg_plus_i, reg_minus_i, reg_plus_j, v_bot_minus_j);
		calculate_nq_x_dir(dims, bcs, iu, iv, q, nq, 2:nx-2, 1:1, x_indexers, y_indexers)
		# y dir
		calculate_nq_y_dir(dims, bcs, iu, iv, q, nq, 2:nx-1, 1:1, x_indexers, y_indexers)

		## RETURN

		return nq
	end

	# calculate entire nq for loop
	function calculate_nq_x_dir(
		dims::Dims,
		bcs::BoundaryConditions,
		iu::Matrix{Int64}, 
		iv::Matrix{Int64}, 
		#ip::Matrix{Float64}, 
		q::Matrix{Float64},
		nq::Matrix{Float64},
		x_range::UnitRange{Int},
		y_range::UnitRange{Int},
		u_idx::Indexer,
		v_idx::Indexer,
	)

		for i = x_range
			for j = y_range
				u_pi = u_idx.pi(q, iu, bcs, i,j)
				u_mi = u_idx.mi(q, iu, bcs, i,j)
				u_pj = u_idx.pj(q, iu, bcs, i,j)
				u_mj = u_idx.mj(q, iu, bcs, i,j)

				v_pi = v_idx.pi(q, iv, bcs, i,j)
				v_mi = v_idx.mi(q, iv, bcs, i,j)
				#v_pj = v_idx.pj(q, iv, bcs, i,j)
				#v_mj = v_idx.mj(q, iv, bcs, i,j)

				nq[iu[i,j]] = (u_pi^2 - u_mi^2 )/dims.dx + (u_pj * v_pi - u_mj * v_mi) / dims.dy
			end
		end
		
	end

	function calculate_nq_y_dir(
		dims::Dims,
		bcs::BoundaryConditions,
		iu::Matrix{Int64}, 
		iv::Matrix{Int64}, 
		#ip::Matrix{Float64}, 
		q::Matrix{Float64},
		nq::Matrix{Float64},
		x_range::UnitRange{Int},
		y_range::UnitRange{Int},
		u_idx::Indexer,
		v_idx::Indexer,
	)

		for i = x_range
			for j = y_range
				#u_pi = u_idx.pi(q, iu, bcs, i,j)
				#u_mi = u_idx.mi(q, iu, bcs, i,j)
				u_pj = u_idx.pj(q, iu, bcs, i,j)
				u_mj = u_idx.mj(q, iu, bcs, i,j)

				v_pi = v_idx.pi(q, iv, bcs, i,j)
				v_mi = v_idx.mi(q, iv, bcs, i,j)
				v_pj = v_idx.pj(q, iv, bcs, i,j)
				v_mj = v_idx.mj(q, iv, bcs, i,j)

				nq[iv[i,j]] = (v_pj^2 - v_mj^2 )/dims.dy + (u_pj * v_pi - u_mj * v_mi) / dims.dx
			end
		end
		
	end

	reg_plus_i(q::Matrix{Float64}, indexer::Matrix{Int64}, bcs::BoundaryConditions, i::Int, j::Int) =
		(1/2) * (q[indexer[i+1,j]] + q[indexer[i,j]])
	reg_minus_i(q::Matrix{Float64}, indexer::Matrix{Int64}, bcs::BoundaryConditions, i::Int, j::Int) =
		(1/2) * (q[indexer[i-1,j]] + q[indexer[i,j]])
	reg_plus_j(q::Matrix{Float64}, indexer::Matrix{Int64}, bcs::BoundaryConditions, i::Int, j::Int) =
		(1/2) * (q[indexer[i,j+1]] + q[indexer[i,j]])
	reg_minus_j(q::Matrix{Float64}, indexer::Matrix{Int64}, bcs::BoundaryConditions, i::Int, j::Int) =
		(1/2) * (q[indexer[i,j-1]] + q[indexer[i,j]])
	
	##
	## TOP AND BOTTOM SPECIAL BOUDNARIES 
	##
	
	# top and bottom BCs for u
	# the boundary value is the average with the ghost node
	u_bot_minus_j(q::Matrix{Float64}, indexer::Matrix{Int64}, bcs::BoundaryConditions, i::Int, j::Int) =
		bcs.u_b[i]
	u_top_plus_j(q::Matrix{Float64}, indexer::Matrix{Int64}, bcs::BoundaryConditions, i::Int, j::Int) =
		bcs.u_t[i]
		
	# top and bottom BCS for v
	# they involve averaging with the boundary value
	v_bot_minus_j(q::Matrix{Float64}, indexer::Matrix{Int64}, bcs::BoundaryConditions, i::Int, j::Int) =
		(1/2) * (bcs.v_b[i] + q[indexer[i,j]])
	v_top_plus_j(q::Matrix{Float64}, indexer::Matrix{Int64}, bcs::BoundaryConditions, i::Int, j::Int) =
		(1/2) * (bcs.v_t[i] + q[indexer[i,j]])

	##
	## LEFT AND RIGHT SPECIAL BOUDNARIES 
	##

	# left and right conditions for v
	v_left_minus_i(q::Matrix{Float64}, indexer::Matrix{Int64}, bcs::BoundaryConditions, i::Int, j::Int) =
		bcs.v_l[j]
	v_right_plus_i(q::Matrix{Float64}, indexer::Matrix{Int64}, bcs::BoundaryConditions, i::Int, j::Int) =
		bcs.v_r[j]

	# right and left conditions for u
	u_left_minus_i(q::Matrix{Float64}, indexer::Matrix{Int64}, bcs::BoundaryConditions, i::Int, j::Int) =
		(1/2) * (bcs.u_l[j] + q[indexer[i,j]])
	u_right_plus_i(q::Matrix{Float64}, indexer::Matrix{Int64}, bcs::BoundaryConditions, i::Int, j::Int) =
		(1/2) * (bcs.v_r[j] + q[indexer[i,j]])

end
