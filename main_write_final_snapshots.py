import argparse
import os

from helpers import get_final_snapshot, get_final_snapshot_filename, get_run_ids_for_n_stars
from amuse.lab import write_set_to_file


def write_final_snapshot(n_stars: int, run_id: int):
    final_snapshot_filename = get_final_snapshot_filename(n_stars=n_stars, run_id=run_id)
    if os.path.exists(final_snapshot_filename):
        print(f'file {final_snapshot_filename} already exists!')
        return
    
    final_snapshot = get_final_snapshot(n_stars=n_stars, run_id=run_id)
    write_set_to_file(final_snapshot, final_snapshot_filename, "amuse")


def write_heads_for_nstars(n_stars: int):
    run_ids = get_run_ids_for_n_stars(n_stars=n_stars)
    for rnid in run_ids:
        write_final_snapshot(n_stars=n_stars, run_id=rnid)


if __name__ in '__main__':
    parser = argparse.ArgumentParser()
    
    parser.add_argument(
        '--nstars', 
        default=16, 
        help='Number of stars in the runs', 
        type=int,
        )
    
    args = parser.parse_args()
    
    write_heads_for_nstars(n_stars=args.nstars)

