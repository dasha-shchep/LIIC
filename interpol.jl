# Functions for LIIC
module Interpol

using DelimitedFiles

export cartesian, distance, internal, xyz2matrix, writeFile

function cartesian(inMat1,inMat2,steps)
    # For an N atom system, this function returns an N x 3 x nsteps array
    difMat = (inMat2 - inMat1) / (steps - 1)
    natoms = size(difMat)[1]
    interpollatedArray = Array{Float32,3}(undef,natoms,3,steps)
    for i in 0:(steps - 1)
        interpollatedArray[:,:,i+1] = inMat1 + i * difMat
    end
    return interpollatedArray
end

function distance()
    println("execute dii")
end

function internal()
    println("execute iii")

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
    natoms  =   size(raw_xyz)[1]
    just_coords = Array{Float32}(raw_xyz[2:end,2:4])
    atom_names = Array{String}(raw_xyz[2:end,1])
    return just_coords, atom_names
end

function xyz2internal(input_xyz)
    # Imports .xyz file format to zmat, returns internal coord vector
    babel_zmat = read(`obabel -ixyz $(input_xyz) -ogzmat`,String)
    regexmatch = match(r"Variables:",babel_zmat)
	segment_zmat = babel_zmat[regexmatch.offset:end]
	internals_array = readdlm(IOBuffer(segment_zmat),skipstart=1)
    return internals_array
end

end
