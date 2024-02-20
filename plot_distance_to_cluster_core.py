import numpy as np
import matplotlib.pyplot as plt

from amuse.lab import Particles
from helpers import get_run_directory
from heritage_from_history import get_hardest_binary_from_decomp
from amuse.lab import read_set_from_file
from amuse.lab import nbody_system

def get_final_binary_distance_in_snapshot(snapshot: Particles, binary_ids: list[int] | list[str]) -> float:
    primary_id, secondary_id = binary_ids
    final_binary = Particles()
    final_binary.add_particle(snapshot[snapshot.id == primary_id])
    final_binary.add_particle(snapshot[snapshot.id == secondary_id])
    distance_to_cc = final_binary.center_of_mass().length().value_in(nbody_system.length)

    return distance_to_cc


def get_final_binary_for_run(n_stars: int, run_id: int):
    directory = get_run_directory(n_stars=n_stars, run_id=run_id)
    history_fn = f'{directory}/history.txt'
    snapshot_fn = f'{directory}/snapshots.log'

    with open(history_fn, 'r') as file:
        final_decomp = file.read().split('\n')[-3].split(' - ')[1].split(':')[1]
        _, final_binary_ids = get_hardest_binary_from_decomp(decomp=final_decomp)
    
    final_snapshot = read_set_from_file(snapshot_fn, copy_history=False)

    return final_snapshot, final_binary_ids


def plot_coredistances():
    n_bodies = 64
    run_ids = range(100)

    bindists = []
    for run_id in run_ids:
        print(run_id)
        final_snapshot, final_binary = get_final_binary_for_run(n_stars=n_bodies, run_id=run_id)
        bindist = get_final_binary_distance_in_snapshot(snapshot=final_snapshot, binary_ids=final_binary)
        bindists.append(bindist)

    plt.hist(bindists)
    plt.show()


if __name__ == '__main__':
    plot_coredistances()
    # get_final_binary_distance_in_snapshot()
    # get_final_binary_for_run(n_stars=16, run_id=0)
