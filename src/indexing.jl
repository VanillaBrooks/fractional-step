module indexing
	import Base

	using .Main.Structs: Dims, create_dims

	export Indexable, check_indexing, Iu, Iv, Ip, IuT, IvT, IpT

	abstract type Indexable end

	struct StoredIndexing{T} <: Indexable
		inner::Matrix{T}
	end

	function Base.getindex(store::StoredIndexing, i::Int64, j::Int64)::Int64
		return store.inner[i,j]
	end

	struct Iu <: Indexable
		dims::Dims
	end

	function Base.getindex(rt::Iu, i::Int64, j::Int64)::Int64
		elem_per_col = rt.dims.ny
		column_num = i - 1
		return (elem_per_col * column_num) + j
	end

	struct Iv <: Indexable
		dims::Dims
	end

	function Base.getindex(rt::Iv, i::Int64, j::Int64)::Int64
		offset = (rt.dims.nx - 1) * rt.dims.ny
		elem_per_col = rt.dims.ny-1
		column_num = i - 1
		return (elem_per_col * column_num) + j + offset
	end

	struct Ip <: Indexable
		dims::Dims
	end

	function Base.getindex(rt::Ip, i::Int64, j::Int64)::Int64
		offset = 0
		elem_per_col = rt.dims.ny
		column_num = i - 1
		return (elem_per_col * column_num) + j + offset
	end

	struct IuT <: Indexable
		dims::Dims
	end

	function Base.getindex(rt::IuT, i::Int64, j::Int64)::Int64
		offset = 0
		elem_per_row = rt.dims.nx-1
		row_num = j - 1
		return (elem_per_row * row_num) + i + offset
	end

	struct IvT <: Indexable
		dims::Dims
	end

	function Base.getindex(rt::IvT, i::Int64, j::Int64)::Int64
		offset = (rt.dims.nx - 1) * rt.dims.ny
		elem_per_row = rt.dims.nx
		row_num = j - 1
		return (elem_per_row * row_num) + i + offset
	end

	struct IpT <: Indexable
		dims::Dims
	end

	function Base.getindex(rt::IpT, i::Int64, j::Int64)::Int64
		offset = 0
		elem_per_row = rt.dims.nx
		row_num = j - 1
		return (elem_per_row * row_num) + i + offset
	end

	function check_indexing()
		nx = 10
		ny = 12
		
		dims = create_dims(1.0, 1.0, nx, ny)
		
		iu, iv = create_iu_iv(dims)
		ip = create_ip(dims)

		iu_linear = Iu(dims)
		iv_linear = Iv(dims)
		ip_linear = Ip(dims)

		check_helper(nx, ny, iu,iv,ip, iu_linear, iv_linear, ip_linear)

		iu, iv = create_iu_iv_aditya(dims)
		ip = create_ip_aditya(dims)

		iu_linear = IuT(dims)
		iv_linear = IvT(dims)
		ip_linear = IpT(dims)

		check_helper(nx, ny, iu,iv,ip, iu_linear, iv_linear, ip_linear)
	end

	function check_helper(nx::Int, ny::Int, 
		iu::StoredIndexing{Int},
		iv::StoredIndexing{Int},
		ip::StoredIndexing{Int},
		iu_linear,
		iv_linear,
		ip_linear,
	)
		#
		# check IP
		#

		for i = 1:nx-1
			for j = 1:ny
				matrix_val = iu[i,j]
				linear_val = iu_linear[i,j]

				println("matrix: ", matrix_val, " linear ", linear_val)

				@assert matrix_val == linear_val
			end
		end

		#
		# Check IV
		#

		for i = 1:nx
			for j = 1:ny-1
				matrix_val = iv[i,j]
				linear_val = iv_linear[i,j]

				println("matrix: ", matrix_val, " linear ", linear_val)

				@assert matrix_val == linear_val
			end
		end

		#
		# Check IP
		#

		for i = 1:nx
			for j = 1:ny
				matrix_val = ip[i,j]
				linear_val = ip_linear[i,j]

				println("matrix: ", matrix_val, " linear ", linear_val)

				@assert matrix_val == linear_val
			end
		end
	end

	function create_iu_iv_aditya(dims::Dims)::Tuple{StoredIndexing{Int}, StoredIndexing{Int}}
		iu = zeros(Int, dims.nx-1, dims.ny)
		iv = zeros(Int, dims.nx, dims.ny-1)

		k = 1
		# init iu values
		for j = 1:dims.ny
			for i = 1:dims.nx -1
				iu[i,j] = k
				k += 1
			end
		end

		# init iv values
		for j = 1:dims.ny -1
			for i = 1:dims.nx
				iv[i,j] = k
				k += 1
			end
		end

		return StoredIndexing(iu),StoredIndexing(iv)
	end

	function create_iu_iv(dims::Dims)::Tuple{StoredIndexing{Int}, StoredIndexing{Int}}
		iu = zeros(Int, dims.nx-1, dims.ny)
		iv = zeros(Int, dims.nx, dims.ny-1)

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

		return StoredIndexing(iu),StoredIndexing(iv)
	end

	function create_ip_aditya(dims::Dims)::StoredIndexing{Int} 
		ip = zeros(Int, dims.nx, dims.ny)

		k = 1
		#for i = 1:dims.nx
		#	for j = 1:dims.ny
		for j = 1:dims.ny
			for i = 1:dims.nx
				ip[i,j] = k
				k += 1
			end
		end
		
		return StoredIndexing(ip)
	end

	function create_ip(dims::Dims)::StoredIndexing{Int} 
		ip = zeros(Int, dims.nx, dims.ny)

		k = 1
		for i = 1:dims.nx
			for j = 1:dims.ny
				ip[i,j] = k
				k += 1
			end
		end
		
		return StoredIndexing(ip)
	end

end
