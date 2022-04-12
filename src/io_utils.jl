module IoUtils
	using .Main.Structs: Dims, BoundaryConditions
	export init_flowfield, FlowfieldWriter, close_file
	using HDF5

	struct FlowfieldWriter{T}
		dims::Dims
		bcs::BoundaryConditions
		h5file::HDF5.File
		velocity_dset::HDF5.Dataset
		pressure_dset::HDF5.Dataset
		time_dset::HDF5.Dataset
		iu
		iv
		ip
		u::Matrix{T}
		v::Matrix{T}
		p::Matrix{T}
		u_tmp::Matrix{T}
		v_tmp::Matrix{T}
	end

	function init_flowfield(dims::Dims, bcs::BoundaryConditions, iu, iv, ip, num_flowfields::Int)
		println("creating a writer efor ", num_flowfields, " flowfields")
		#
		# matricies
		#

		u_temp_store = zeros(dims.nx-1, dims.ny)
		v_temp_store = zeros(dims.nx, dims.ny-1)
		u_matrix = zeros(dims.nx, dims.ny)
		v_matrix = zeros(dims.nx, dims.ny)
		p_matrix = zeros(dims.nx, dims.ny)

		#
		# hdf5 init stuff
		#
		filename = "./results/flowfield.h5"
		Base.Filesystem.rm(filename, force=true)
		flowfield_file_id = h5open(filename, "w")
		velocity_dataset = create_dataset(flowfield_file_id, "velocity", datatype(Float64), dataspace(num_flowfields, 2, dims.nx, dims.ny))
		pressure_dataset = create_dataset(flowfield_file_id, "pressure", datatype(Float64), dataspace(num_flowfields, dims.nx, dims.ny))
		time_dataset = create_dataset(flowfield_file_id, "times", datatype(Float64), dataspace(num_flowfields, 1))

		return FlowfieldWriter(
			dims, bcs, 
			flowfield_file_id, velocity_dataset, pressure_dataset, time_dataset,
			iu, iv, ip, 
			u_matrix, v_matrix, p_matrix,
			u_temp_store, v_temp_store
		)
	end

	function close_file(writer::FlowfieldWriter{T}) where T <: AbstractFloat
		close(writer.h5file)
	end

	function write_flowfield(writer::FlowfieldWriter{T}, q::Vector{T}, p::Vector{T}, step::Int, t::T) where T <: AbstractFloat
		println("writing io field number ", step)

		nx = writer.dims.nx
		ny = writer.dims.ny


		for i = 1:nx-1
			for j = 1:ny
				writer.u_tmp[i,j] = q[writer.iu[i,j]]
			end
		end

		for i = 1:nx
			for j = 1:ny-1
				writer.v_tmp[i,j] = q[writer.iv[i,j]]
			end
		end

		for i = 1:nx
			for j = 1:ny
				writer.p[i,j] = p[writer.ip[i,j]]
			end
		end

		writer.u[2:nx-1, 1:ny] = (writer.u_tmp[2:nx-1, 1:ny] .+ writer.u_tmp[1:nx-2, 1:ny]) ./ 2
		writer.u[1, 1:ny] = (writer.u_tmp[1, 1:ny] + writer.bcs.u_l[1:ny]) ./ 2
		writer.u[nx, 1:ny] = (writer.u_tmp[nx-1, 1:ny] + writer.bcs.u_r[1:ny]) ./ 2

		writer.v[1:nx, 2:ny-1] = (writer.v_tmp[1:nx, 2:ny-1] .+ writer.v_tmp[1:nx, 1:ny-2]) ./ 2
		writer.v[1:nx, 1] = (writer.v_tmp[1:nx, 1] + writer.bcs.v_b[1:nx]) ./ 2
		writer.v[1:nx, ny] = (writer.v_tmp[1:nx, ny-1] + writer.bcs.v_t[1:nx]) ./ 2

		# write the results to the output arrays for post processessing
		writer.velocity_dset[step, 1, :, :] = writer.u
		writer.velocity_dset[step, 2, :, :] = writer.v
		writer.pressure_dset[step, :, :] = writer.p
		writer.time_dset[step, 1] = t
	end
end
