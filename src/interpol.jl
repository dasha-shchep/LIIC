# Functions for LIIC
module Interpol

using DelimitedFiles
export cartesian, distance, internal, xyz2matrix, writeFile

function cartesian(inMat1,inMat2,steps)
	"""
    For an N atom system, this function returns a 3D, N x 3 x nsteps, array
	"""
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

function internal(zmat_f,zmat_l,steps)
	int_f = zmat_f.IntVars
	int_l = zmat_l.IntVars
	natoms = zmat_l.Number
	nvars = length(zmat_l.IntVars)
	difference = Vector(undef,nvars)
	for i in 1:nvars
		toAngle = (zmat_l.IntVars[i]+180.)
		fromAngle = (zmat_f.IntVars[i]+180.)
		diff = abs(toAngle - fromAngle)
		if diff < 180.
			difference[i] = toAngle - fromAngle
		else
			if fromAngle > toAngle
				fromAngle = fromAngle - 360.
				difference[i] = toAngle - fromAngle
			elseif toAngle > fromAngle
				toAngle = toAngle - 360.
				difference[i] = toAngle - fromAngle
			end
		end
		# diff = zmat_l.IntVars[i]-zmat_f.IntVars[i]
		# if -180. < diff < 180.
		# 	difference[i] = diff
		# elseif diff < 180.
		# 	difference[i] = diff + 360.
		# elseif diff > 180.
		# 	difference[i] = diff - 360.
		# else
		# 	println("ERROR")
		# end
	end
	difference = difference / (steps-1)
	liic = Array{Float64,2}(undef,nvars,steps)
	arrayOfZMats = Array{ZMatrix, 1}(undef, steps)
	for step in 1:steps
		liic[:,step] = int_f + (step-1)*difference
		arrayOfZMats[step] = ZMatrix(zmat_f.Atoms,liic[:,step],zmat_f.VarNames,natoms)
	end
	@assert isapprox(int_l,liic[:,steps],atol=1e-5)
	return arrayOfZMats
end

# OLD IMPLEMENTATION THAT USED OPEN BABEL
function internal_babel(int_arr_1,int_arr_2,steps,header)
	"""
	Takes array of Z-matrix internal coordinates and interpolates between them
	to return an array of xyzs.
	"""
	numInternalCrd = size(int_arr_1)[1]
	natoms = Int((numInternalCrd + 6)/3)
	diffVec = Vector(undef,numInternalCrd)
	for i in 1:numInternalCrd
		if ( -180. <= (int_arr_2[i,2]-int_arr_1[i,2]) <= 180. )
			diffVec[i] = (int_arr_2[i,2]-int_arr_1[i,2])/steps
		elseif (int_arr_2[i,2] - int_arr_1[i,2]) < -180.
			diffVec[i] = (int_arr_2[i,2] - int_arr_1[i,2] + 360.)/steps
		elseif (int_arr_2[i,2]-int_arr_1[i,2]) > 180.
			diffVec[i] = (int_arr_2[i,2] - int_arr_1[i,2] - 360.)/steps
		else
			println("Something has gone terribly wrong")
		end
	end
	zmatArray = Vector(undef,steps)
	interXYZ = Array{Float32,3}(undef,natoms,3,steps)
	atom_names = Array{String,1}(undef,natoms)
	for j in 1:steps
		zmatArray[j] = header * "Variables:\n" * intlVec(int_arr_1,diffVec,j)
		open("tmp","w") do io
			write(io,zmatArray[j])
		end
		convXyz = read(`obabel -igzmat tmp -oxyz`,String)
		run(`rm tmp`)
		just_coords = Array{Float32}(readdlm(IOBuffer(convXyz),skipstart=2)[:,2:4])
		atom_names = Array{String}(readdlm(IOBuffer(convXyz),skipstart=2)[:,1])
		interXYZ[:,:,j] = just_coords
	end
	return interXYZ, atom_names
end
# END OF OLD IMPLEMENTATION

function intlVec(int_arr_1,diffVec,step)
	"""
	Function generates set of internal coordinates as a string at a particular
	step during the interpollation to concatenate with a header string to output
	a Z-matrix format string.
	"""
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
	"""
    Imports .xyz file format as data frame, converts it to an N x 3 Matrix
	"""
    raw_xyz =   readdlm(input_xyz)
    just_coords = Array{Float32}(raw_xyz[2:end,2:4])
    atom_names = Array{String}(raw_xyz[2:end,1])
    return just_coords, atom_names
end

function xyz2internal(input_xyz)
	"""
	Imports .xyz file format to zmat, returns internal coord vector
	"""
    babel_zmat = read(`obabel -ixyz $(input_xyz) -ogzmat`,String)
    regexmatch = match(r"Variables:",babel_zmat)
	coord_segment_zmat = babel_zmat[regexmatch.offset:end]
	header_segment = babel_zmat[1:(regexmatch.offset-1)]
	internals_array = readdlm(IOBuffer(coord_segment_zmat),skipstart=1)
    return internals_array, header_segment
end

end
