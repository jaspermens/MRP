import argparse

from run_decomp_alice import make_buncha_data


if __name__ in '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--nstars', 
        default=16, 
        help='Number of stars in the runs', 
        )
    parser.add_argument(
        '--nchkpt', 
        default=128, 
        help='checkpoint frequency',
        )
    parser.add_argument(
        '--nruns', 
        default=10, 
        help='Number of runs',
        )
    parser.add_argument(
        '--min_hardness_kt', 
        default=1, 
        help='Minimum hardness for binary membership',
        )
    args = parser.parse_args()
    
    make_buncha_data(n_runs=args.nruns, n_stars=args.nstars, snapshot_frequency=args.nchkpt, min_hardness_kt=args.min_hardness_kt)

