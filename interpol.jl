# Functions for LIIC
module Interpol

using DelimitedFiles

export cartesian
export distance
export internal
export xyz2matrix

function cartesian(inMat1,inMat2,steps)
    # For an N atom system, this function returns an N x nsteps x 3 array
    difMat = (inMat2 - inMat1) / (steps - 1)
    interpollatedArray = Array{Float64,3}
    return interpollatedArray
end

function distance()
    println("execute dii")
end

function internal()
    println("execute iii")
end

function xyz2matrix(input_xyz)
    # Imports .xyz file format as data frame, converts it to an N x 3 Matrix
    raw_xyz =   readdlm(input_xyz)
    natoms  =   size(raw_xyz)[1]
    just_coords = Array{Float64}(raw_xyz[2:natoms,2:4])
    return just_coords
end

end
