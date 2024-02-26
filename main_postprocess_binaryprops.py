import argparse

from postprocess_fb_properties import get_final_binary_properties_for_n_stars

if __name__ in '__main__':
    parser = argparse.ArgumentParser()
    
    parser.add_argument(
        '--nstars', 
        default=16, 
        help='Number of stars in the runs', 
        type=int,
        )
    
    args = parser.parse_args()
    
    get_final_binary_properties_for_n_stars(n_stars=args.nstars)

