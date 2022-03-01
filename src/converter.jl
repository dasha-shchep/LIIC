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
    readline(io);
    reg=r"([a-zA-Z]+)\s+(-?\d+.?\d*)\s+(-?\d+.?\d*)\s+(-?\d+.?\d*)"
    atom_names=Vector{String}()
    atom_coords=Array{Float64,2}(undef,num_atom,3)
    for num in 1:num_atom
        m=match(reg,readline(io))
        push!(atom_names,string(m.captures[1]::SubString))
        atom_coords[num,:]=[parse(Float64,m.captures[2]),parse(Float64,m.captures[3]),parse(Float64,m.captures[4])]
    end
    close(io)
    return Molecule(atom_names,atom_coords,num_atom)
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