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
    bond1 = a3 - a2
    bond2 = a2 - a1
    sintheta = norm(cross(bond1,bond2))
    costheta = dot(bond1,bond2)
    theta = atan(sintheta,costheta)
    angle = 180.0 * theta/pi
    # angle = acosd(dot(bond1,bond2)/(norm(bond1)*norm(bond2)))
    return angle::Float64
end

function dihedral(a1::Vector{Float64},a2::Vector{Float64},a3::Vector{Float64},a4::Vector{Float64})
    bond1 = a4-a3
    bond2 = a3-a2
    bond3 = a2-a1
    plane1 = cross(bond1,bond2)
    plane2 = cross(bond2,bond3) 
    dihedral = acosd(dot(plane1,plane2)/(norm(plane1)*norm(plane2)))
    if dot(cross(plane1,plane2),bond2) < 0
        dihedral = -dihedral
    end
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

function rotation_matrix(axis,angle)
    # Rotation about the origin in 3 dimensions using the Euler-Rodrigues formula
    axis = -axis / norm(axis)
    a = cosd(angle*0.5)
    (b,c,d) = axis * sind(angle*0.5)
    rotM = [a*a+b*b-c*c-d*d 2*(b*c-a*d) 2*(b*d+a*c);
    2*(b*c+a*d) a*a+c*c-b*b-d*d 2*(c*d-a*b); 
    2*(b*d-a*c) 2*(c*d+a*b) a*a+d*d-b*b-c*c]
    return rotM::Matrix{Float64}
end

function zmat_to_xyz(zmat::ZMatrix)
    natoms=zmat.Number
    atom_coords = Array{Float64,2}(undef,natoms,3)

    if natoms > 0
        atom_coords[1,:]=[0.0,0.0,0.0]

        if natoms > 1
            r2 = zmat.IntVars[1]
            atom_coords[2,:]=[r2,0.0,0.0]

            if natoms > 2
                vec = atom_coords[2,:] - atom_coords[1,:]
                r3 = zmat.IntVars[2]
                a3 = zmat.IntVars[3]
                vec = r3 * vec / norm(vec)
                atom_coords[3,:] = atom_coords[2,:] + (rotation_matrix([0, 0, 1], a3) * vec)

                if natoms > 3
                    for atom in 4:natoms
                        rn = zmat.IntVars[3*atom-8]
                        an = zmat.IntVars[3*atom-7]
                        dn = zmat.IntVars[3*atom-6]

                        # ox = rn * cosd(an)
                        # oy = rn * cosd(dn) * sind(an)
                        # oz = rn * sind(an) * sind(dn)

                        # Coordinates of the last three atoms
                        am1 = atom_coords[atom-1,:]
                        am2 = atom_coords[atom-2,:]
                        am3 = atom_coords[atom-3,:]

                        # Vectors between them
                        vec1 = am1 - am2
                        vec2 = am2 - am3
                        # Plane normal defined by these vectors
                        nm = cross(vec1,vec2)
                        
                        # Vector of length rn pointing along vec1
                        vecrn = rn * vec1 / norm(vec1)

                        # Rotate vecrn by angle an around the normal
                        vecan = rotation_matrix(nm,an)*vecrn
                        
                        # Rotate vecan around dihedral angle
                        vecdn = rotation_matrix(vec1,dn)*vecan
                        
                        # Add this vector to previous coordinate
                        atom_coords[atom,:] = am1 + vecdn
                    end
                end 
            end
        end
    end

    mol=Molecule(zmat.Atoms,atom_coords,natoms)

    return mol::Molecule
end

end