# Julia LIIC Script
# Written on 26 November 2019
# Execute by running `julia liic.jl`

using DelimitedFiles
using LinearAlgebra
using ArgParse

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

function xyz2matrix(input_xyz)
    # Imports .xyz file format as data frame, converts it to an N x 3 Matrix
    raw_xyz=readdlm(input_xyz)
    natoms=size(raw_xyz)[1]
    just_coords=Array{Float64}(raw_xyz[2:natoms,2:4])
    return just_coords
end

function cartesian_interpolation()
    println("Performing cartesian interpolation...")
end

function distance_interpolation()
    println("execute dii")
end

function internal_interpolation()
    println("execute iii")
end

function main()
    parsed_args = parse_args(ARGS, parse_settings)
    geometry = parsed_args["geom1"]
    println(geometry)
    if parsed_args["cartesian"]
        cartesian_interpolation()
        input_matrix_a = xyz2matrix(parsed_args["geom1"])
        # show(input_matrix_a)
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
