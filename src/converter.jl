#Utility to convert between xyz and internal (Z-matrix) coordinates
module Converter
"""
Contains functions that convert between xyz and z-matrix formats to serve as
input in the LIIC routines.
"""

using LinearAlgebra
using Formatting

export xyz_to_zmat, zmat_to_xyz

struct Molecule
    Atoms::Vector{String}
    Coord::Array{Float64,2}
    Number::Int64
end

struct ZMatrix
    Atoms::Vector{String}
    IntVars::Vector{Float64}
    VarNames::Vector{String}
    Number::Int64
end

function import_molecule(xyz_file::String)
    io = open(xyz_file)
    natoms = parse(Int,readline(io))
    readline(io);
    reg=r"([a-zA-Z]+)\s+(-?\d+.?\d*)\s+(-?\d+.?\d*)\s+(-?\d+.?\d*)"
    atoms=Vector{String}()
    atom_coords=Array{Float64,2}(undef,natoms,3)
    for num in 1:natoms
        m=match(reg,readline(io))
        push!(atoms,string(m.captures[1]::SubString))
        atom_coords[num,:]=[parse(Float64,m.captures[2]),parse(Float64,m.captures[3]),parse(Float64,m.captures[4])]
    end
    close(io)
    return Molecule(atoms,atom_coords,natoms)
end

function bond(a1::Vector{Float64},a2::Vector{Float64})
    return norm(a1-a2)::Float64
end

function angle(a1::Vector{Float64},a2::Vector{Float64},a3::Vector{Float64})
    bond1 = a1-a2
    bond2 = a3-a2
    angle = acosd(dot(bond1,bond2)/(norm(bond1)*norm(bond2)))
    return angle::Float64
end

function dihedral(a1::Vector{Float64},a2::Vector{Float64},a3::Vector{Float64},a4::Vector{Float64})
    bond1 = a2-a1
    bond2 = a3-a2
    bond3 = a4-a3
    plane1 = cross(bond1,bond2)
    plane2 = cross(bond2,bond3)
    normalize = norm(plane1)*norm(plane2)
    dihedral = acosd(dot(plane1,plane2)/normalize)
    return dihedral::Float64
end

function write_zmat(zmat::ZMatrix,withvars::Bool=true)
    if withvars==true
        natoms = length(zmat.Atoms)
        for atom in 1:natoms
            if atom < 2
                printfmtln("{:3s}",zmat.Atoms[atom])
            elseif atom < 3
                printfmtln("{:3s} {:>4d}  {:11s}",zmat.Atoms[atom],(atom-1),zmat.VarNames[atom-1])
            elseif atom < 4
                printfmtln("{:3s} {:>4d}  {:11s} {:>4d}  {:11s}",zmat.Atoms[atom],(atom-1),
                zmat.VarNames[atom-1],(atom-2),zmat.VarNames[atom])
            else
                printfmtln("{:3s} {:>4d}  {:11s} {:>4d}  {:11s} {:>4d}  {:11s}",zmat.Atoms[atom],
                (atom-1),zmat.VarNames[3*atom-8],(atom-2),zmat.VarNames[3*atom-7],(atom-3),zmat.VarNames[3*atom-6])
            end
        end
        println("  Variables:")
        for var in 1:length(zmat.IntVars)
            printfmtln("{:3s} {:>11.8f}",zmat.VarNames[var],zmat.IntVars[var])
        end
    else
        natoms = length(zmat.Atoms)
        for atom in 1:natoms
            if atom < 2
                printfmtln("{:3s}",zmat.Atoms[atom])
            elseif atom < 3
                printfmtln("{:3s} {:>4d}  {:>11.8f}",zmat.Atoms[atom],(atom-1),zmat.IntVars[atom-1])
            elseif atom < 4
                printfmtln("{:3s} {:>4d}  {:>11.8f} {:>4d}  {:>11.8f}",zmat.Atoms[atom],(atom-1),
                zmat.IntVars[atom-1],(atom-2),zmat.IntVars[atom])
            else
                printfmtln("{:3s} {:>4d}  {:>11.8f} {:>4d}  {:>11.8f} {:>4d}  {:>11.8f}",zmat.Atoms[atom],
                (atom-1),zmat.IntVars[3*atom-8],(atom-2),zmat.IntVars[3*atom-7],(atom-3),zmat.IntVars[3*atom-6])
            end
        end
    end
end    

function xyz_to_zmat(molecule)
    natoms     = molecule.Number
    atoms      = molecule.Atoms
    intvars      = Vector{Float64}()
    varnames   = Vector{String}()
    for atom in 1:natoms
        if atom > 1
            bl     = bond(molecule.Coord[atom-1,:],molecule.Coord[atom,:])
            push!(intvars,bl)
            push!(varnames,"R$atom")
        end
        if atom > 2
            agl  = angle(molecule.Coord[atom-2,:],molecule.Coord[atom-1,:],molecule.Coord[atom,:])
            push!(intvars,agl)
            push!(varnames,"A$atom")
        end
        if atom > 3
            dhl  = dihedral(molecule.Coord[atom-3,:],molecule.Coord[atom-2,:],molecule.Coord[atom-1,:],molecule.Coord[atom,:])
            push!(intvars,dhl)
            push!(varnames,"D$atom")
        end
    end
    zmatrix    = ZMatrix(atoms,intvars,varnames,natoms)
    return zmatrix
end

function zmat_to_xyz(zmat::ZMatrix)
    xyz=Array{Float64,2}
    return mol::Molecule
end

end