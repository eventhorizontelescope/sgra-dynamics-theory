stride=64  # Save VIDA output every 16th .h5 file
blur=15.0  # Blur the images by 15.0 uas
base_dir=/bd6/eht/Illinois_SgrA_v3check/230GHz/ # Directory of the v3 library
filenames_dir=./cache/GRMHD_filepaths/ # Directory to save the filenames.txt file
results_dir=./cache/IDFE_results/   # Directory to save the VIDA output for filenames.txt file

# Create the filenames.txt files for the v3 library
julia make_filenames_txt.jl --base_dir $base_dir --filenames_dir $filenames_dir

# Run VIDA on the filenames.txt files for the v3 library
julia run_vida.jl --filenames_dir $filenames_dir --results_dir $results_dir --stride $stride --blur $blur

# If you want to use SLURM for each filename in filenames.txt, you can use the following command

#julia run_vida_slurm.jl --filenames_dir $filenames_dir --results_dir $results_dir --partition long --time 01:00:00 --cpus 4
