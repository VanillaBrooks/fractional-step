module divergence
	using .Main.Structs: Dims, BoundaryConditions
	using .Main.indexing: Indexable

	export div_

	struct Indexer
		bottom_y::Function
		center_y::Function
		left_x::Function
		center_x::Function
	end

	function div_(
		dims::Dims, 
		bcs::BoundaryConditions, 
		p::Vector{Float64},
		iu::IU,
		iv::IV,
		ip::IP,
	)::Vector{Float64} where IU <: Indexable where IV <: Indexable where IP <: Indexable
		div = zeros(dims.np)

		# indexing in this routine is based on the locations of the pressure
		# since div is indexed exclusively on ip[ , ]

		nx = dims.nx
		ny = dims.ny

		###
		### Inner domain
		###
		
		indexer = Indexer(reg_bottom, reg_center, reg_left, reg_center)
		calculate_divergence(dims, bcs, p, div, 2:nx-1, 2:ny-1, iu, iv, ip, indexer)

		###
		### Bottom Inner
		###

		indexer = Indexer(y_dir_bottom_bc, reg_center, reg_left, reg_center)
		calculate_divergence(dims, bcs, p, div, 2:nx-1, 1:1, iu, iv, ip, indexer)

		###
		### Bottom Left (Pinned)
		###

		div[ip[1,1]] = 0.0

		###
		### Bottom Right
		###

		indexer = Indexer(y_dir_bottom_bc, reg_center, reg_left, x_dir_right_bc)
		calculate_divergence(dims, bcs, p, div, nx:nx, 1:1, iu, iv, ip, indexer)

		###
		### Left Inner
		###
		
		indexer = Indexer(reg_bottom, reg_center, x_dir_left_bc, reg_center)
		calculate_divergence(dims, bcs, p, div, 1:1, 2:ny-1, iu, iv, ip, indexer)

		###
		### Right Inner
		###
		
		indexer = Indexer(reg_bottom, reg_center, reg_left, x_dir_right_bc)
		calculate_divergence(dims, bcs, p, div, nx:nx, 2:ny-1, iu, iv, ip, indexer)

		###
		### Top Inner
		###
		
		indexer = Indexer(reg_bottom, y_dir_top_bc, reg_left, reg_center)
		calculate_divergence(dims, bcs, p, div, 2:nx-1, ny:ny, iu, iv, ip, indexer)

		###
		### Top Left
		###
		
		indexer = Indexer(reg_bottom, y_dir_top_bc, x_dir_left_bc, reg_center)
		calculate_divergence(dims, bcs, p, div, 1:1, ny:ny, iu, iv, ip, indexer)

		###
		### Top Right
		###
		
		indexer = Indexer(reg_bottom, y_dir_top_bc, reg_left, x_dir_right_bc)
		calculate_divergence(dims, bcs, p, div, nx:nx, ny:ny, iu, iv, ip, indexer)

		###
		### RETURN
		###

		return div
	end

	@inline function calculate_divergence(
		dims::Dims,
		bcs::BoundaryConditions,
		p::Vector{Float64},
		div::Vector{Float64},
		x_range::UnitRange{Int64},
		y_range::UnitRange{Int64},
		iu::IU,
		iv::IV,
		ip::IP,
		indexer::Indexer,
	) where IU <: Indexable where IV <: Indexable where IP <: Indexable
		for i = x_range
			for j = y_range 
				center_x = indexer.center_x(p, bcs, iu, i, j)
				center_y = indexer.center_y(p, bcs, iv, i, j)

				left_x = indexer.left_x(p, bcs, iu, i,j)
				bottom_y = indexer.bottom_y(p, bcs, iv, i,j)

				# backward difference in each direction
				diff_x = (center_x - left_x) / dims.dx
				diff_y = (center_y - bottom_y) / dims.dy

				div[ip[i,j]] = diff_x + diff_y
			end
		end
	end

	####
	#### GENERAL INDEXERS - can be used wherever we are not hitting a boundary
	####

	reg_center(p::Vector{Float64}, bcs::BoundaryConditions, indexer, i::Int64, j::Int64)::Float64 = 
		p[indexer[i,j]]
	reg_left(p::Vector{Float64}, bcs::BoundaryConditions, indexer, i::Int64, j::Int64)::Float64 = 
		p[indexer[i-1,j]]
	reg_bottom(p::Vector{Float64}, bcs::BoundaryConditions, indexer, i::Int64, j::Int64)::Float64 = 
		p[indexer[i,j-1]]

	####
	#### Y Velocity boundary conditions for top / bottom
	####

	y_dir_bottom_bc(p::Vector{Float64}, bcs::BoundaryConditions, indexer, i::Int64, j::Int64)::Float64 = 
		bcs.v_b[i]

	y_dir_top_bc(p::Vector{Float64}, bcs::BoundaryConditions, indexer, i::Int64, j::Int64)::Float64 = 
		bcs.v_t[i]

	####
	#### X Velocity boundary conditions for top / bottom
	####
	
	x_dir_left_bc(p::Vector{Float64}, bcs::BoundaryConditions, indexer, i::Int64, j::Int64)::Float64 = 
		bcs.u_l[j]

	x_dir_right_bc(p::Vector{Float64}, bcs::BoundaryConditions, indexer, i::Int64, j::Int64)::Float64 = 
		bcs.u_r[j]
		
end
