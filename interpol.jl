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

function writeFile(output_array,file)
    scanLength = size(output_array)[3]
    atom_number = size(output_array)[1]
    println(typeof(atom_number))
    open(file,"w") do io
        for i in 1:scanLength
            write(io, string(atom_number))
            writedlm(io, "\n")
            writedlm(io, output_array[:,:,i],'\t')
        end
    end
end

function xyz2matrix(input_xyz)
    # Imports .xyz file format as data frame, converts it to an N x 3 Matrix
    raw_xyz =   readdlm(input_xyz)
    natoms  =   size(raw_xyz)[1]
    just_coords = Array{Float32}(raw_xyz[2:end,2:4])
    return just_coords
end

end
