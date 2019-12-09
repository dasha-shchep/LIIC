# Functions for LIIC
module Interpol

using DelimitedFiles

export cartesian, distance, internal, xyz2matrix, writeFile

function cartesian(inMat1,inMat2,steps)
    # For an N atom system, this function returns an N x 3 x nsteps array
    diffMat = (inMat2 - inMat1) / (steps - 1)
    natoms = size(diffMat)[1]
    interpollatedArray = Array{Float32,3}(undef,natoms,3,steps)
    for i in 0:(steps - 1)
        interpollatedArray[:,:,i+1] = inMat1 + i * diffMat
    end
    return interpollatedArray
end

function distance()
    println("execute dii")
end

function internal(int_arr_1,int_arr_2,steps,header)
	diffVec = (int_arr_2[:,2]-int_arr_1[:,2])/steps
	intlArr = Vector(undef,steps)
	for j in 1:steps
		intlArr[j] = header * "Variables:\n" * intlVec(int_arr_1,diffVec,j)
	end
	return intlArr
	# open("temp","w") do io
	# for i in 0:steps
	# 	internals = string(int_arr_1[:,1]," ",(int_arr_1[:,2]+i*diffVec),"\n"))
	# 	print(internals)
		# write(io,header)
		# writedlm(io,internals)
		# babel_xyz = read(`obabel -igzmat "temp" -oxyz`,String)
	# end
	# end
end

function intlVec(int_arr_1,diffVec,step)
	intlCrd = Vector(undef,length(diffVec))
	for i in 1:length(diffVec)
		intlCrd[i] = string(int_arr_1[i,1]," ",(int_arr_1[i,2]+step*diffVec[i]),"\n")
	end
	return join(intlCrd)
end

function writeFile(output_array,atom_names,file)
    scanLength = size(output_array)[3]
    atom_number = size(output_array)[1]
    open(file,"w") do io
        for i in 1:scanLength
	    geometry = hcat(atom_names,output_array[:,:,i])
            write(io, string(atom_number))
            writedlm(io, "\n")
            writedlm(io, geometry,'\t')
        end
    end
end

function xyz2matrix(input_xyz)
    # Imports .xyz file format as data frame, converts it to an N x 3 Matrix
    raw_xyz =   readdlm(input_xyz)
    just_coords = Array{Float32}(raw_xyz[2:end,2:4])
    atom_names = Array{String}(raw_xyz[2:end,1])
    return just_coords, atom_names
end

function xyz2internal(input_xyz)
    # Imports .xyz file format to zmat, returns internal coord vector
    babel_zmat = read(`obabel -ixyz $(input_xyz) -ogzmat`,String)
    regexmatch = match(r"Variables:",babel_zmat)
	coord_segment_zmat = babel_zmat[regexmatch.offset:end]
	header_segment = babel_zmat[1:(regexmatch.offset-1)]
	internals_array = readdlm(IOBuffer(coord_segment_zmat),skipstart=1)
    return internals_array, header_segment
end

end
