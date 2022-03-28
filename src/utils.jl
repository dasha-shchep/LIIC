module Utils
# Module containing the necessary structs and functions for LIIC execution

using LinearAlgebra
using Statistics
using Formatting

export ZMatrix, Molecule, import_molecule, write_file

struct ZMatrix
    Atoms::Vector{String}
    IntVars::Vector{Float64}
    VarNames::Vector{String}
    Number::Int64
end

struct Molecule
    Atoms::Vector{String}
    Coord::Array{Float64,2}
    Number::Int64
end

function translate_to_centroid(coord_matrix)
    # Normalises the molecular coordinates by centering them.
    center = mean(coord_matrix,dims=1)
    translated_geom = coord_matrix .- center
    return translated_geom
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

function kabsch_rotate(Pxyz,Qxyz)
    # Returns the two geometries in xyz format, one of which has been translated, the
    # other translated and rotated 
    Pgeom = import_molecule(Pxyz)
    Qgeom = import_molecule(Qxyz)

    normalisedP = (translate_to_centroid(Pgeom.Coord))
    normalisedQ = (translate_to_centroid(Qgeom.Coord))
    
    # Calculate cross covariance matrix
    xcov = transpose(normalisedP)*normalisedQ
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
        natoms = zmat.Number
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

function write_xyz(mol::Molecule,io=nothing)
    natoms = mol.Number
    if io !== nothing
        printfmtln(io,"{:>4d}",natoms)
        println(io)
        for atom in 1:natoms
            printfmtln(io,"{:3s}   {:>11.8f}   {:>11.8f}    {:>11.8f}",mol.Atoms[atom],mol.Coord[atom,1],mol.Coord[atom,2],mol.Coord[atom,3])
        end
    elseif io === nothing
        printfmtln("{:>4d}",natoms)
        println()
        for atom in 1:natoms
            printfmtln("{:3s}   {:>11.8f}   {:>11.8f}    {:>11.8f}",mol.Atoms[atom],mol.Coord[atom,1],mol.Coord[atom,2],mol.Coord[atom,3])
        end
    else
        println("Error in write_xyz")
    end
end

function write_liic(outArray,filename::String)
    io = open(filename,"w")
    for i in 1:length(outArray)
        write_xyz(outArray[i],io)
    end
    close(io)
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

function cartesian(mol1::Molecule,mol2::Molecule,steps::Int64)
	"""
    This function returns an "steps" length array of Molecule objects
	"""
    # difference = (mol2.Coord - mol1.Coord) / (steps - 1)
    natoms = mol1.Number
	liic = Array{Float64,3}(undef,natoms,3,steps)
    arrayOfMolecules = Array{Molecule, 1}(undef, steps)
    for step in 1:steps
		fracA = (steps-step)/(steps-1)
		fracB = (step-1)/(steps-1)
		liic[:,:,step] = mol1.Coord*fracA + mol2.Coord*fracB
        arrayOfMolecules[step] = Molecule(mol1.Atoms,liic[:,:,step],natoms) 
    end
    return arrayOfMolecules
end

function internal_babel(int_arr_1,int_arr_2,steps,header)
	"""
	Takes array of Z-matrix internal coordinates and interpolates between them
	to return an array of xyzs.
	"""
	numInternalCrd = size(int_arr_1)[1]
	natoms = Int((numInternalCrd + 6)/3)
	diffVec = Vector(undef,numInternalCrd)
	for i in 1:numInternalCrd
		if ( -180. <= (int_arr_2[i,2]-int_arr_1[i,2]) <= 180. )
			diffVec[i] = (int_arr_2[i,2]-int_arr_1[i,2])/steps
		elseif (int_arr_2[i,2] - int_arr_1[i,2]) < -180.
			diffVec[i] = (int_arr_2[i,2] - int_arr_1[i,2] + 360.)/steps
		elseif (int_arr_2[i,2]-int_arr_1[i,2]) > 180.
			diffVec[i] = (int_arr_2[i,2] - int_arr_1[i,2] - 360.)/steps
		else
			println("Something has gone terribly wrong")
		end
	end
	zmatArray = Vector(undef,steps)
	interXYZ = Array{Float32,3}(undef,natoms,3,steps)
	atom_names = Array{String,1}(undef,natoms)
	for j in 1:steps
		zmatArray[j] = header * "Variables:\n" * intlVec(int_arr_1,diffVec,j)
		open("tmp","w") do io
			write(io,zmatArray[j])
		end
		convXyz = read(`obabel -igzmat tmp -oxyz`,String)
		run(`rm tmp`)
		just_coords = Array{Float32}(readdlm(IOBuffer(convXyz),skipstart=2)[:,2:4])
		atom_names = Array{String}(readdlm(IOBuffer(convXyz),skipstart=2)[:,1])
		interXYZ[:,:,j] = just_coords
	end
	return interXYZ, atom_names
end

end