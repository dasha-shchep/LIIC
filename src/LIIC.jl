# Julia LIIC & other useful scripts
# Written on 26 November 2019
# Update on 28 March 2022
# Execute by running `julia liic.jl --flags path_to_geom1 path_to_geom2`

using ArgParse

include("utils.jl")
import .Utils

# include("interpol.jl")
# include("kabsch.jl")
# include("converter.jl")
# include("structs.jl")
# using .Interpol, .Kabsch, .Converter, ..Structs

parse_settings = ArgParseSettings()

@add_arg_table parse_settings begin
    "--cartesian", "-c"
        help = "Specifies cartesian coordinate interpolation"
        action = :store_true
    "--distance", "-d"
        help = "Specifies distance coordinate interpolation"
        action = :store_true
    "--internal", "-i"
        help = "Specifies interatomic distance interpolation"
        action = :store_true
    "--steps", "-s"
        help = "Number of steps between structures in interpolation"
        arg_type = Int
        default = 20
    "--output", "-o"
        help = "Path to output .xyz file of the LIIC"
        arg_type = String
        default = "liic_output.xyz"
    "--kabsch", "-k"
        help = "Whether kabsch should run before interpolation"
        action = :store_true
    "geom1"
        help = "path to .xyz file containing first geometry"
        required = true
    "geom2"
        help = "path to .xyz file containing first geometry"
        required = true
end

function main()
    parsed_args = parse_args(ARGS, parse_settings)
    stp = parsed_args["steps"]
    file = parsed_args["output"]
    in_file_1 = parsed_args["geom1"]
    in_file_2 = parsed_args["geom2"]
    MoleculeA = Utils.import_molecule(in_file_1)
    MoleculeB = Utils.import_molecule(in_file_2)
    @assert MoleculeA.Atoms == MoleculeB.Atoms
    if parsed_args["cartesian"]
        println("Performing cartesian interpolation...")
        # in_mat_1, atomNames1 = Interpol.xyz2matrix(in_file_1)
        # in_mat_2, atomNames2 = Interpol.xyz2matrix(in_file_2)
        # Perform kabsch algorithm on geom2 if kabsch==True
        if parsed_args["kabsch"]
            println("Running Kabsch algorithm on second geometry prior to LIIC")
            MoleculeA,MoleculeB = Utils.kabsch_rotate(MoleculeA,MoleculeB)
        end
        outputArray = Utils.cartesian(MoleculeA,MoleculeB,stp)
    elseif parsed_args["distance"]
        println("Performing interpolation in internal distance matrix...")
    elseif parsed_args["internal"]
        println("Performing interpolation in internal coordinates...")
        ZMatA = Utils.xyz_to_zmat(MoleculeA)
        ZMatB = Utils.xyz_to_zmat(MoleculeB)
        @assert ZMatA.VarNames == ZMatB.VarNames
        outputArray = Utils.internal(ZMatA,ZMatB,stp)
		# internalCoords1,header = Utils.xyz2internal(in_file_1)
		# internalCoords2,header = Utils.xyz2internal(in_file_2)
		# @assert internalCoords1[:,1] == internalCoords2[:,1]
		# outputArray, atomNames1 = Interpol.internal_babel(internalCoords1,internalCoords2,stp,header)
    else
        println("ERROR")
    end

    Utils.write_liic(outputArray,file)

end

# Off we ... off we ... off we fucking GO!
main()
