#Utility to convert between xyz and internal (Z-matrix) coordinates
module Converter
"""
Contains functions that convert between xyz and z-matrix formats to serve as
input in the LIIC routines.

"""

export xyz_to_zmat, zmat_to_xyz

struct Molecule
    Atoms::Vector{String}
    Coord::Array{Float64,2}
    Number::Int64
end

function import_molecule(xyz_file)
    io = open(xyz_file)
    num_atom = parse(Int,readline(io))
    readline(io)
    # string_xyz = read(io,String)
    # Insert regex to get data
    close(io)
    return Molecule(atom_names,coordinates,num_atom)
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

function xyz_to_zmat(molecule)
    zmat_array = []
    return zmat_array
end

function zmat_to_xyz()
    xyz_array = 3
    return xyz_array
end

end