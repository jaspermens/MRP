import numpy as np
from amuse.lab import Particles, nbody_system
from helpers import custom_tqdm, read_snapshot_file, get_run_directory
from decompose_multiples import find_composite_multiples
from run_decomp_alice import get_max_hardness_in_decomp_list
import os
from amuse.io.base import IoException

def get_quantities_for_snapshot(snapshot: Particles, initial_ke, header=False) -> str | tuple[str]:
    if header:
        return 'time - T/T_0 - hardest_hardness - decomposition'
    
    time = snapshot.get_timestamp()
    ids = find_composite_multiples(plummer=snapshot.copy(), initial_ke=initial_ke, min_hardness_kt=1)
    decomp = [id for id in ids if len(id) > 5]
    relative_temperature = snapshot.kinetic_energy()/initial_ke
    hardest_hardness = get_max_hardness_in_decomp_list(decomp_list=decomp)

    return (
            f'{time.value_in(nbody_system.time):.3f}', 
            f'{relative_temperature:.3f}', 
            f'{hardest_hardness:.3f}', 
            str(decomp),
            )


def redo_decomp_for_run(run_id: int, n_stars: int) -> None:
    run_directory = get_run_directory(n_stars=n_stars, run_id=run_id)
    decomp_filename = f'{run_directory}/decompositions.txt'

    snapshots = read_snapshot_file(run_directory=run_directory)
    initial_ke = snapshots[0].kinetic_energy()

    header = get_quantities_for_snapshot(snapshot=None, initial_ke=None, header=True)
    csvlines = np.zeros(shape=(len(snapshots), len(header.split(' - '))), dtype='U2000')

    for i, snapshot in enumerate(custom_tqdm(snapshots, total=len(snapshots))):
        csvlines[i,:] = get_quantities_for_snapshot(snapshot=snapshot, initial_ke=initial_ke)

    np.savetxt(decomp_filename, csvlines, fmt='%5s - %5s - %5s - %s', header=header)


def claim_run(run_id: int, n_stars: int):
    run_directory = get_run_directory(n_stars=n_stars, run_id=run_id)
    claimed_flag = f'{run_directory}/claimed.txt'

    os.mknod(claimed_flag) # will raise a FileExistsError if the file already exists
    if os.path.exists(f'{run_directory}/decompositions.txt'):
        raise FileExistsError


def redo_decomp(n_stars: int, start_num: int = 0):
    run_id = start_num
    while True:
        try:
            # try to claim the run
            claim_run(run_id=run_id, n_stars=n_stars)
            # if successful, redo the decomposition (and put the result in the claimed file)
            redo_decomp_for_run(run_id=run_id, n_stars=n_stars)
            run_id += 1

        except FileExistsError:
            # if the file already exists, then a run is either done or in progress
            # move on to the next one
            run_id += 1

        except FileNotFoundError:
            # then the run directory doesn't exist, so we're done (or something is borked)
            print('end of the line! No more runs to do!')
            break


if __name__ == '__main__':
    redo_decomp(n_stars=16, start_num=599)
