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
		calc = final_calc_u
		calculate_nq(dims, bcs, iu, iv, q, nq, 2:nx-2, 2:ny-1,x_indexers, y_indexers, calc)

		# calc = final_calc_v
		# calculate_nq(dims, bcs, iu, iv, q, nq, 2:nx-1, 2:ny-2,x_indexers, y_indexers, calc)

		# ##
		# ## Bottom Inner
		# ##
		# 
		# x_indexers = Indexer(reg_plus_i, reg_minus_i, reg_plus_j, u_bot_bc_i);
		# y_indexers = Indexer(reg_plus_i, reg_minus_i, reg_plus_j, reg_minus_j);
		# calculate_nq(dims, bcs, iu, iv, q, nq, 2:nx-2, 1:1,x_indexers, y_indexers, final_calc_u)
		# calculate_nq(dims, bcs, iu, iv, q, nq, 2:nx-2, 1:1,x_indexers, y_indexers, final_calc_v)


		## RETURN

		return nq
	end

	# calculate entire nq for loop
	function calculate_nq(
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
		final_calc::Function
	)

		for i = x_range
			for j = y_range
				u_pi = u_idx.pi(q, iu, bcs, i,j)
				u_mi = u_idx.mi(q, iu, bcs, i,j)
				u_pj = u_idx.pj(q, iu, bcs, i,j)
				u_mj = u_idx.mj(q, iu, bcs, i,j)

				v_pi = v_idx.pi(q, iv, bcs, i,j)
				v_mi = v_idx.mi(q, iv, bcs, i,j)
				v_pj = v_idx.pj(q, iv, bcs, i,j)
				v_mj = v_idx.mj(q, iv, bcs, i,j)

				# update the nq value at this point based on the generic function
				final_calc(nq, iu, iv, dims, i, j, u_pi, u_mi, u_pj, u_mj, v_pi, v_mi, v_pj, v_mj)
			end
		end
		
	end


	function final_calc_u(
		nq::Matrix{Float64},
		iu::Matrix{Int64},
		iv::Matrix{Int64},
		dims::Dims,
		i::Int,
		j::Int,
		u_pi::Float64,
		u_mi::Float64,
		u_pj::Float64,
		u_mj::Float64,
		v_pi::Float64,
		v_mi::Float64,
		v_pj::Float64,
		v_mj::Float64)

		nq[iu[i,j]] = (u_pi^2 - u_mi^2 )/dims.dx + (u_pj * v_pi - u_mj * v_mi) / dims.dy
	end

	function final_calc_v(
		nq::Matrix{Float64},
		iu::Matrix{Int64},
		iv::Matrix{Int64},
		dims::Dims,
		i::Int,
		j::Int,
		u_pi::Float64,
		u_mi::Float64,
		u_pj::Float64,
		u_mj::Float64,
		v_pi::Float64,
		v_mi::Float64,
		v_pj::Float64,
		v_mj::Float64)

		nq[iv[i,j]] = (v_pj^2 - v_mj^2 )/dims.dy + (u_pj * v_pi - u_mj * v_mi) / dims.dx
	end

	reg_plus_i(q::Matrix{Float64}, indexer::Matrix{Int64}, bcs::BoundaryConditions, i::Int, j::Int) =
		(1/2) * (q[indexer[i+1,j]] + q[indexer[i,j]])
	reg_minus_i(q::Matrix{Float64}, indexer::Matrix{Int64}, bcs::BoundaryConditions, i::Int, j::Int) =
		(1/2) * (q[indexer[i-1,j]] + q[indexer[i,j]])
	reg_plus_j(q::Matrix{Float64}, indexer::Matrix{Int64}, bcs::BoundaryConditions, i::Int, j::Int) =
		(1/2) * (q[indexer[i,j+1]] + q[indexer[i,j]])
	reg_minus_j(q::Matrix{Float64}, indexer::Matrix{Int64}, bcs::BoundaryConditions, i::Int, j::Int) =
		(1/2) * (q[indexer[i,j-1]] + q[indexer[i,j]])
	
	v_top_bc_i(q::Matrix{Float64}, indexer::Matrix{Int64}, bcs::BoundaryConditions, i::Int, j::Int) =
		bcs.v_t[i]
	v_top_bc_j(q::Matrix{Float64}, indexer::Matrix{Int64}, bcs::BoundaryConditions, i::Int, j::Int) =
		bcs.v_t[j]

	v_top_bc_minus_i(q::Matrix{Float64}, indexer::Matrix{Int64}, bcs::BoundaryConditions, i::Int, j::Int) =
		bcs.v_t[i-1]
	v_top_bc_plus_i(q::Matrix{Float64}, indexer::Matrix{Int64}, bcs::BoundaryConditions, i::Int, j::Int) =
		bcs.v_t[i+1]

	u_bot_bc_i(q::Matrix{Float64}, indexer::Matrix{Int64}, bcs::BoundaryConditions, i::Int, j::Int) =
		bcs.u_b[i]
		

end
