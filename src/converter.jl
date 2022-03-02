#Utility to convert between xyz and internal (Z-matrix) coordinates
module Converter
"""
Contains functions that convert between xyz and z-matrix formats to serve as
input in the LIIC routines.
"""

using LinearAlgebra

export xyz_to_zmat, zmat_to_xyz

struct Molecule
    Atoms::Vector{String}
    Coord::Array{Float64,2}
    Number::Int64
end

struct ZMatrix
    Atoms::Vector{String}
    Bonds::Vector{Float64}
    Angles::Vector{Float64}
    Dihedrals::Vector{Float64}
    VarNames::Vector{String}
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

function bond(a1::Vector{Float64},a2::Vector{Float64})
    return norm(a1-a2)
end

function angle(a1::Vector{Float64},a2::Vector{Float64},a3::Vector{Float64})
    bond1 = a1-a2
    bond2 = a3-a2
    angle = acosd(dot(bond1,bond2)/(norm(bond1)*norm(bond2)))
    return angle
end

function dihedral(a1::Vector{Float64},a2::Vector{Float64},a3::Vector{Float64},a4::Vector{Float64})
    bond1 = a2-a1
    bond2 = a3-a2
    bond3 = a4-a3
    plane1 = cross(bond1,bond2)
    plane2 = cross(bond2,bond3)
    normalize = norm(plane1)*norm(plane2)
    dihedral = acosd(dot(plane1,plane2)/normalize)
    return dihedral
end

function write_zmat(zmatrix,withvars)
    if withvars==false
        println("bleh")
    else
        println("blah")
    end
end    

function xyz_to_zmat(molecule)
    natoms     = molecule.Number
    atoms      = molecule.Atoms
    bonds      = Vector{Float64}()
    angles     = Vector{Float64}()
    dihedrals  = Vector{Float64}()
    varnames   = Vector{String}()
    for atom in 1:(natoms-1)
        bl     = bond(molecule.Coord[atom,:],molecule.Coord[atom+1,:])
        push!(bonds,bl)
        # if natoms > 2
        #     agl  = angle(molecule.Coord[atom,:],molecule.Coord[atom+1,:],molecule.Coord[atom+2,:])
        #     push!(angles,agl)
        #     if natoms > 3
        #         dhl  = dihedral(molecule.Coord[atom,:],molecule.Coord[atom+1,:],molecule.Coord[atom+2,:])
        #         push!(dihedrals,dhl)
        #     end
        # end
    end
    zmatrix    = ZMatrix(atoms,bonds,angles,dihedrals,varnames)
    return zmatrix
end

function zmat_to_xyz()
    xyz_array = 3
    return xyz_array
end

end