# Module to perform optional kabsch rotation to minimise distance between two structures before interpolating between them
module Kabsch

using DelimitedFiles

using LinearAlgebra
using Statistics

export kabsch_rotate

function translate_to_centroid(coord_matrix)
    # Normalises the molecular coordinates by centering them.
    center = [mean(coord_matrix[:,1]);mean(coord_matrix[:,2]);mean(coord_matrix[:,3])]
    centroid = transpose(center)
    translated_geom = broadcast(-,coord_matrix,centroid)
    return translated_geom
end

function cross_covariance_matrix(Pmatrix,Qmatrix)
    # Cross covariance matrix gives measure of variability between two matrices
    CCmatrix = transpose(Pmatrix) * Qmatrix
    return CCmatrix
end

function optimal_rotation_matrix(CCmatrix)
    # Returns 3x3 matrix that can be applied to P to get Q
    # Previous implementation, not general (for cases of non-invertible matrices)
    # ORmatrix = sqrt(transpose(CCmatrix)*CCmatrix)*inv(CCmatrix)
    # Using SVD instead
    u,sing,vt = svd(CCmatrix)
    if sign(det(transpose(vt)*transpose(u))) == 1.0
        ORmatrix = transpose(vt)*transpose(u)
    elseif sign(det(transpose(vt)*transpose(u))) == -1.0
        mat = [1. 0. 0.;0. 1. 0.; 0. 0. -1.]
        ORmatrix = transpose(vt)*mat*transpose(u)
    else
        println("Error: Issue with SVD routine in finding optimal rotation matrix")
    end
    return ORmatrix
end

function xyz2matrix(input_xyz)
    # Imports .xyz file format as data frame, and converts it to an N x 3 Matrix
    raw_xyz=readdlm(input_xyz)
    natoms=size(raw_xyz)[1]
    just_coords=Array{Float64}(raw_xyz[2:natoms,2:4])
    return just_coords
end

function kabsch_rotate(Pgeom,Qgeom)
    # Returns the two geometries in xyz format, one of which has been translated, the
    # other translated and rotated 
    # Pgeom = xyz2matrix(p_xyz)
    # Qgeom = xyz2matrix(q_xyz)

    normalisedP = (translate_to_centroid(Pgeom))
    normalisedQ = (translate_to_centroid(Qgeom))
    
    xcov = cross_covariance_matrix(normalisedP,normalisedQ)
    orot = optimal_rotation_matrix(xcov)
    
    num_atoms = size(normalisedP)[1]
    rotated = zeros(Float64,num_atoms,3)
    
    for i=1:num_atoms
        rotated[i,:] = orot*normalisedP[i,:]
    end

    RMSD_value = norm(rotated-normalisedQ)

    println("The RMSD between these two structures is ",RMSD_value)

    return rotated, normalisedQ
end  

end



