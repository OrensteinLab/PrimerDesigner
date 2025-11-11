from General.args import *
from General.utils import *
from PD_var_ILP.run_var_ilp import *
from PD_mul_ILP.run_mul_ilp import *
from PD_mul_greedy.run_mul_greedy import *
from PrimerDesigner.PD_single_LPath.run_pd_single_LPath import *



def main():

    # parse user arguments
    args = get_args()

    # read protein coding-sequences from file
    mutreg_regions, full_sequences, protein_names = read_sequences(args.file_path)
    
    if args.version =="Non_relaxed":
        print("Running Non-relaxed version")
        run_mul_ilp(mutreg_regions,full_sequences,protein_names,args)
    elif args.version == "PD-mul-nh":
        print("Running PD-mul-nh version")
        run_mul_greedy(full_sequences, mutreg_regions, protein_names,args)
    elif args.version=="PD-single":
        print("Running PD-single version")
        
    else: # PD-mul-var is default
        print("Running PD-single-LPath version")
        run_shortest_path(full_sequences[0], mutreg_regions[0], protein_names[0],args)  # only runs on first sequence
        


if __name__ == '__main__':
    main()


