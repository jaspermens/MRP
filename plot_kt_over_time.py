import numpy as np
import matplotlib.pyplot as plt

from amuse.lab import read_set_from_file
from amuse.lab import Particles
from amuse.lab import nbody_system

from helpers import get_run_directory

from tqdm import tqdm


def get_ke_history_for_run(n_bodies: int, run_number: int) -> tuple[np.array, np.array]:
    run_directory = get_run_directory(n_stars=n_bodies, run_id=run_number)
    filename = f'{run_directory}/decompositions.txt'
    times, ke_history, *_ = np.genfromtxt(filename, delimiter=' - ', dtype=str, skip_header=1, usecols=(0, 1), unpack=True)
    return times.astype(float), ke_history.astype(float)

def get_ke_history_depr(n_bodies, run_number):
    filename = f'/home/jasper/studie/MRP/output/n{n_bodies}/run_{run_number:>04}/snapshots.log'
    all_snapshots = read_set_from_file(filename=filename, copy_history=False, close_file=False)
    print('snapshots read')
    times = []
    ke_history = []
    c = 16
    for snapshot in tqdm(all_snapshots.history):
        if c < 16:
            c += 1
            continue
        ke_history.append(get_ke_from_snapshot(snapshot))
        times.append(snapshot.get_timestamp().value_in(nbody_system.time))
        c = 0

    return np.array(times), np.array(ke_history) / ke_history[0]


def get_ke_from_snapshot(bodies: Particles):
    return bodies.kinetic_energy().value_in(nbody_system.energy)


def plot_ke_history(n_bodies: int = 16, run_number: int = 2):
    fig, ax = plt.subplots(1, 1, figsize=[6,4])
    times, ke_hist = get_ke_history_for_run(n_bodies=n_bodies, run_number=run_number)
    ax.plot(times, ke_hist, c='black')
    ax.set_xlabel('Time (nbody)')
    ax.set_ylabel(r'$T/T_0$')
    plt.show()

def plot_buncha_ke_history():
    n_runs = 20
    fig, (ax16, ax64) = plt.subplots(2, 1, sharex=True)
    for run_id in range(n_runs):
        ax16.plot(*get_ke_history_for_run(n_bodies=64, run_number=run_id), c='black', alpha=.3)
        ax64.plot(*get_ke_history_for_run(n_bodies=16, run_number=run_id), c='black', alpha=.3)

    ax16.set_title('N=16')
    # ax64.set_title('N=64')
    ax64.set_xlabel('Time [nbody]')

    plt.show()

def plot_fourier_ke_history(n_bodies: int = 16, run_number: int = 2):
    fig, ax = plt.subplots(1, 1, figsize=[6,4])
    times, ke_hist = get_ke_history_for_run(n_bodies=n_bodies, run_number=run_number)
    freqs, pgram = get_lombscargle(t=times, m=ke_hist/ke_hist[0])
    ax.semilogx(freqs, pgram, c='black')
    # ax.set_xlabel('Time (nbody)')
    # ax.set_ylabel(r'$T/T_0$')
    plt.show()


def get_lombscargle(t, m):
    from scipy.signal import lombscargle

    freqs = np.logspace(np.log10(0.01), np.log10(5), num=10_000)

    pgram = lombscargle(t, m, freqs*2*np.pi, precenter=True, normalize=False)

    return freqs, pgram


if __name__ == '__main__':
    # plot_log_ke_history(run_number=3)
    plot_ke_history(run_number=0, n_bodies=16)

