from General.args import *
from General.read_sequences import *
from Relaxed.run_relaxed import *
from Non_relaxed.run_non_relaxed import *
from Extension.run_extension import *


def main():

    # parse user arguments
    args = get_args()

    # read protein coding-sequences from file
    mutreg_regions, full_sequences, protein_names = read_sequences(args.file_path)

    if args.version =="Non_relaxed":
        run_non_relaxed(mutreg_regions,full_sequences,protein_names,args)
    elif args.version == "Extension":
        run_extension(full_sequences, mutreg_regions, protein_names,args)
    else: # relaxed version is default
        run_relaxed_ilp(full_sequences[0],mutreg_regions[0],protein_names[0],args) # only runs on first sequence

if __name__ == '__main__':
    main()


