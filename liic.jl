# Julia LIIC Script
# Written on 26 November 2019
# Execute by running `julia liic.jl --flags path_to_geom1 path_to_geom2`

using LinearAlgebra
using ArgParse

include("interpol.jl")
import .Interpol

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
        default = 10
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
    in_file_1 = parsed_args["geom1"]
    in_file_2 = parsed_args["geom2"]
    if parsed_args["cartesian"]
        println("Performing cartesian interpolation...")
        in_mat_1 = Interpol.xyz2matrix(in_file_1)
        in_mat_2 = Interpol.xyz2matrix(in_file_2)
        foo = Interpol.cartesian(in_mat_1,in_mat_2,stp)
        show(foo)
    elseif parsed_args["distance"]
        println("Performing interpolation in internal distance matrix...")
    elseif parsed_args["internal"]
        println("Performing interpolation in internal coordinates...")
    else
        println("what")
    end

    # println(parsed_args["cartesian"])
    # println("Parsed args:")
    # for (arg,val) in parsed_args
    #     println("  $arg  =>  $val")
    # end
end

# Off we ... off we ... off we fucking GO!
main()
