import numpy as np
import matplotlib.pyplot as plt
from amuse.lab import nbody_system

from helpers import read_history_csv, custom_tqdm, get_run_ids_for_n_stars, snapshot_at_time, read_snapshot_file, get_run_directory, get_output_path
from plotconfig import *

def get_core_radius_history_for_run(n_stars: int, run_id: int):
    output_path = get_output_path()
    filename = f'{output_path}/n{n_stars}/run_{run_id:>04}/core_radius_n{n_stars}_run{run_id}.npy'
    try:
        time, core_radius = np.load(file=filename)
        return time, core_radius
    
    except FileNotFoundError:
        pass
    snapshots = read_snapshot_file(run_directory=get_run_directory(n_stars=n_stars, run_id=run_id))
    rcs = []
    times = []
    for snapshot in custom_tqdm(snapshots[::10], total=len(snapshots[::10])):
        core_radius = snapshot.densitycentre_coreradius_coredens(number_of_neighbours=4)[1].value_in(nbody_system.length)
        time = snapshot.get_timestamp().value_in(nbody_system.time)
        rcs.append(core_radius)
        times.append(time)

    np.save(file=filename, arr=np.array([times, rcs]))

    return np.array([times, rcs])


def plot_core_radius_history():
    n_stars = 16    
    run_ids = range(3)

    fig, ax = plt.subplots(1,1, figsize=[8,6])

    for run_id in run_ids:
        times, rcs = get_core_radius_history_for_run(n_stars=n_stars, run_id=run_id)
        ax.plot(times - times[-1], rcs, lw=1, alpha=.5)

    ax.set_xlabel('time until core collapse [nbody]')

    plt.show()


if __name__ == '__main__':
    plot_core_radius_history()