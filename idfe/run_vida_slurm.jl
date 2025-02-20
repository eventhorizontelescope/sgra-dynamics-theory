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

        "--partition"
            help = "SLURM partition to submit jobs to"
            arg_type = String
            default = "short"

        "--time"
            help = "Maximum runtime for each job (e.g., 00:30:00 for 30 minutes)"
            arg_type = String
            default = "00:30:00"

        "--cpus"
            help = "Number of CPUs per job"
            arg_type = Int
            default = 1
    end

    return parse_args(s)
end

args = parse_commandline()

function submit_job(file, results_filepath, stride, blur, partition, time, cpus)
    slurm_script = """
    #!/bin/bash
    #SBATCH --job-name=vida_pol
    #SBATCH --output=${results_filepath}.log
    #SBATCH --error=${results_filepath}.err
    #SBATCH --partition=$partition
    #SBATCH --time=$time
    #SBATCH --cpus-per-task=$cpus

    module load julia
    julia vida_pol.jl --imfiles $file --outname $results_filepath --stride $stride --blur $blur
    """

    slurm_filename = "submit_$(basename(file)).slurm"
    open(slurm_filename, "w") do f
        write(f, slurm_script)
    end

    run(`sbatch $slurm_filename`)
end

function main()
    filenames_dir = args["filenames_dir"]
    results_dir = args["results_dir"]
    stride = args["stride"]
    blur = args["blur"]
    partition = args["partition"]
    time = args["time"]
    cpus = args["cpus"]

    # Ensure results directory exists
    mkpath(results_dir)

    # Get sorted list of .txt files in the filenames directory
    files = sort(glob("*.txt", filenames_dir))

    # Submit separate SLURM jobs for each file
    for file in files
        results_filename = replace(basename(file), ".txt" => ".csv")
        results_filepath = joinpath(results_dir, results_filename)
        submit_job(file, results_filepath, stride, blur, partition, time, cpus)
    end
end

main()