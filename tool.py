from General.args import *
from General.read_sequences import *
from PD_mul_var.run_relaxed import *
from Non_relaxed.run_non_relaxed import *
from PD_mul_nh.run_mul_nh import *
from PD_single.run_pd_single import *


def main():

    # parse user arguments
    args = get_args()

    # read protein coding-sequences from file
    mutreg_regions, full_sequences, protein_names = read_sequences(args.file_path)

    if args.version =="Non_relaxed":
        print("Running Non-relaxed version")
        run_non_relaxed(mutreg_regions,full_sequences,protein_names,args)
    elif args.version == "PD-mul-nh":
        print("Running PD-mul-nh version")
        run_extension(full_sequences, mutreg_regions, protein_names,args)
    elif args.version=="PD-single":
        print("Running PD-single version")
        run_single(full_sequences[0], mutreg_regions[0], protein_names[0],args)  # only runs on first sequence
    else: # PD-mul-var is default
        print("Running PD-mul-var version")
        run_relaxed_ilp(full_sequences[0],mutreg_regions[0],protein_names[0],args) # only runs on first sequence

if __name__ == '__main__':
    main()


