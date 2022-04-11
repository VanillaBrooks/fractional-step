module gradient
	using .Main.Structs: Dims, BoundaryConditions
	using .Main.indexing: Indexable
	using Printf
	import Base.Threads.@threads

	export grad
	
	function grad!(
		dims::Dims, 
		p::Vector{Float64}, 
		iu::IU,
		iv::IV,
		ip::IP,
		grad::Vector{Float64}
		) where IU <: Indexable where IV <: Indexable where IP <: Indexable

		#grad = zeros(dims.nu)

		nx = dims.nx
		ny = dims.ny

		## dp / dx

		@threads for i = 2:nx
			for j = 1:ny
				sum = p[ip[i,j]] - p[ip[i-1,j]]

				grad[iu[i-1,j]] = sum / dims.dx
			end
		end

		# dp / dy
	
		@threads for i = 1:nx
			for j = 2:ny
				sum = p[ip[i,j]] - p[ip[i,j-1]]

				grad[iv[i,j-1]] = sum / dims.dy
			end
		end

		#return grad
	end
end
