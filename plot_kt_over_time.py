import numpy as np
import matplotlib.pyplot as plt

from amuse.lab import read_set_from_file
from amuse.lab import Particles
from amuse.lab import nbody_system

from tqdm import tqdm


def get_ke_history(n_bodies, run_number):
    filename = f'/home/jasper/studie/MRP/output/n{n_bodies}/run_{run_number:>04}/snapshots.log'
    all_snapshots = read_set_from_file(filename=filename, copy_history=False, close_file=False)
    print('snapshots read')
    times = []
    ke_history = []
    for snapshot in tqdm(all_snapshots.history):
        ke_history.append(get_ke_from_snapshot(snapshot))
        times.append(snapshot.get_timestamp().value_in(nbody_system.time))
    return np.array(times), np.array(ke_history)


def get_ke_from_snapshot(bodies: Particles):
    return bodies.kinetic_energy().value_in(nbody_system.energy)


def plot_ke_history(n_bodies: int = 16, run_number: int = 2):
    fig, ax = plt.subplots(1, 1, figsize=[6,4])
    times, ke_hist = get_ke_history(n_bodies=n_bodies, run_number=run_number)
    ax.plot(times, ke_hist / ke_hist[0], c='black')
    ax.set_xlabel('Time (nbody)')
    ax.set_ylabel(r'$T/T_0$')
    plt.show()


def plot_log_ke_history(n_bodies: int = 16, run_number: int = 2):
    fig, ax = plt.subplots(1, 1, figsize=[6,4])
    times, ke_hist = get_ke_history(n_bodies=n_bodies, run_number=run_number)
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
    plot_ke_history(run_number=3, n_bodies=64)

