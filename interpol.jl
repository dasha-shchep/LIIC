# Functions for LIIC
module Interpol

export cartesian_interpolation

function cartesian_interpolation(inMat1,inMat2,steps)
    # For an N atom system, this function returns an N x nsteps x 3 array
    difMat = (inMat2 - inMat1) / (steps - 1)
    return difMat
end

end
