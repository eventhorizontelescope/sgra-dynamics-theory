using Pkg; Pkg.activate(@__DIR__)
using Glob
using ArgParse

function parse_command_line()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--base_dir"
            help = "Base directory containing the data files"
            arg_type = String
            default = "/bd6/eht/Illinois_SgrA_v3check/230GHz/"

        "--filenames_dir"
            help = "Directory where output .txt files will be saved"
            arg_type = String
            default = "./filenames/"
    end

    return parse_args(ARGS, s)
end

args = parse_command_line()

function main()
    base_dir = args["base_dir"]*"/"
    filenamedir = args["filenames_dir"]*"/"

    # Ensure the filenames directory exists
    mkpath(filenamedir)

    # Define parameter values
    windows = ["w3", "w4", "w5"]
    fluxes = ["Ma", "Sa"]
    spins = ["-0.94", "-0.5", "0", "+0.5", "+0.94"]
    Rh_values = ["Rh1", "Rh10", "Rh40", "Rh160"]
    inclinations = ["i10", "i30", "i50", "i70", "i90", "i110", "i130", "i150", "i170"]

    # Loop over all parameter combinations
    for window in windows, flux in fluxes, spin in spins, Rhigh in Rh_values, inclination in inclinations
        # Construct the model file pattern
        model = "$(flux)$(spin)_$(window)/img_*_$(Rhigh)_$(inclination).h5"
        
        # Generate the output filename
        model_txt = "$(flux)$(spin)_$(window)_$(Rhigh)_$(inclination).txt"
        output_filename = joinpath(filenamedir, model_txt)

        # Get sorted file paths
        files = sort(glob(model, base_dir))

        # Write to the output file if there are matching files
        if !isempty(files)
            open(output_filename, "w") do f
                for file in files
                    println(f, file)
                end
            end
            println("Saved: $output_filename")
        else
            println("No files found for: $model")
        end
    end
end

main()