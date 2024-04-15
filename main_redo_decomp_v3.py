import argparse

from postprocess_v3 import redo_decomp

if __name__ in '__main__':
    parser = argparse.ArgumentParser()
    
    parser.add_argument(
        '--nstars', 
        default=16, 
        help='Number of stars in the runs', 
        type=int,
        )
    
    parser.add_argument(
        '--start_id', 
        default=0, 
        help='Which run to start with (optional)-  alllows you to skip pre-redone rons',
        type=int,
        )
    
    args = parser.parse_args()
    
    redo_decomp(n_stars=args.nstars, start_num=args.start_id)

