from General.args import *
from General.utils import *
from PD_var_ILP.run_var_ilp import *
from PD_mul_ILP.run_mul_ilp import *
from PD_mul_greedy.run_mul_greedy import *
from PD_single_LPath.run_pd_single_LPath import *



def main():

    # parse user arguments
    args = get_args()

    # read protein coding-sequences from file
    mutreg_regions, full_sequences, protein_names = read_sequences(args.file_path)
    
    if args.version =="PD-mul-ILP":
        print("Running PD-mul-IL version")
        run_mul_ilp(mutreg_regions,full_sequences,protein_names,args)
    elif args.version == "PD-mul-Greedy":
        print("Running PD-mul-Greedy version")
        run_mul_greedy(full_sequences, mutreg_regions, protein_names,args)
    elif args.version=="PD-var-ILP":
        print("Running PD-var-ILP version")
        run_var_ilp(full_sequences[0], mutreg_regions[0], protein_names[0],args)
    else: 
        print("Running PD-single-LPath version")
        run_longest_path(full_sequences[0], mutreg_regions[0], protein_names[0],args)  # only runs on first sequence
        


if __name__ == '__main__':
    main()


