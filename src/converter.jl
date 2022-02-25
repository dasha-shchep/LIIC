#Utility to convert between xyz and internal (Z-matrix) coordinates
module Converter
"""
Contains functions that convert between xyz and z-matrix formats to serve as
input in the LIIC routines.

"""

export xyz_to_zmat, zmat_to_xyz

struct Molecule{I,F} 
    Atoms::Vector
    Coord::Array
    Number::Int64
end

function bond(a1,a2)
    return 0
end

function angle(a1,a2,a3)
    return 0
end

function dihedral(a1,a2,a3,a4)
    return 0
end

function xyz_to_zmat(xyz_filename)
    zmat_array = []
    return zmat_array
end

function zmat_to_xyz()
    xyz_array = 3
    return xyz_array
end

end