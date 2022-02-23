# Module to perform optional kabsch rotation to minimise distance between two structures before interpolating between them
module Kabsch

using DelimitedFiles

using LinearAlgebra
using Statistics

export run_kabsch

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
    ORmatrix = sqrt(transpose(CCmatrix)*CCmatrix)*inv(CCmatrix)
    return ORmatrix
end

function run_kabsch(Pgeom,Qgeom)
    # Returns the two geometries in xyz format, one of which has been translated, the
    # other translated and rotated
    

end  

end

