using Pkg; Pkg.activate(@__DIR__)
using Glob
using ArgParse

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--filenames_dir"
            help = "Directory containing the input .txt filenames"
            arg_type = String
            default = "filenames/"

        "--results_dir"
            help = "Directory where output .csv files will be saved"
            arg_type = String
            default = "results/"

        "--stride"
            help = "Stride parameter for vida_pol.jl"
            arg_type = Int
            default = 16

        "--blur"
            help = "Blur parameter (in uas) for vida_pol.jl"
            arg_type = Float64
            default = 15.0
    end

    return parse_args(s)
end

args = parse_commandline()

function main()
    filenames_dir = args["filenames_dir"]*"/"
    results_dir = args["results_dir"]*"/"
    stride = args["stride"]
    blur = args["blur"]

    # Ensure the results directory exists
    mkpath(results_dir)

    # Get sorted list of .txt files in the filenames directory
    files = sort(glob("*.txt", filenames_dir))

    # Loop over each file and process it
    for file in files
        results_filename = replace(basename(file), ".txt" => ".csv")
        results_filepath = joinpath(results_dir, results_filename)

        cmd = `julia vida_pol.jl --imfiles $file --outname $results_filepath --stride $stride --blur $blur`
        run(cmd)
    end
end

main()