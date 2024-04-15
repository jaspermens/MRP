from plotconfig import *
import numpy as np
import matplotlib.pyplot as plt
from helpers import get_run_directory, dists_for_run


def get_dist_chunk(dists, t_start, chunk_size):
    start_id = np.minimum(t_start * 128, len(dists))
    end_id = np.minimum((t_start + chunk_size) * 128, len(dists))
    chunk = dists[start_id:end_id]
    pad_width = ((0, chunk_size*128 - (end_id-start_id)), (0, 0))
    padded = np.pad(chunk, pad_width=pad_width, mode='constant', constant_values=0, )
    return padded

def get_dists_at_time(n_stars: int, t_start, chunk_size):
    n_runs = 150

    chunks = np.zeros((n_runs, chunk_size * 128, n_stars))
    for run_id in range(n_runs):
        run_directory = get_run_directory(run_id=run_id, n_stars=n_stars)
        dists = np.load(f"{run_directory}/distances.npy")
        chunk = get_dist_chunk(dists=dists, t_start=t_start, chunk_size=chunk_size)   
        chunks[run_id] = chunk

    return chunks


def hist_dists_over_time() -> None:
    chunk_size = 10
    window_step = 1
    n_stars = 16

    for t in np.arange(0, 200, window_step):
        dists = get_dists_at_time(n_stars=n_stars, t_start=t, chunk_size=chunk_size)
        cleaned = dists[dists > 0].flatten()
        plt.hist(cleaned, bins=np.arange(0, 20, .1))
        plt.savefig(f'plots/distance_hists/n{n_stars}/time{t}.png')
        plt.close()


def step_dists() -> None:
    num_runs = 100

    distances16 = np.zeros((num_runs, 16))
    distances64 = np.zeros((num_runs, 64))

    fig, ax = plt.subplots()
    for run_id in range(num_runs):
        distances16[run_id] = dists_for_run(run_id, 16)[-1]
        distances64[run_id] = dists_for_run(run_id, 64)[-1]

        # ax.step(x = np.linspace(0,16,16), y=np.sort(distances16[-10]), c='black', alpha=.1)
        # ax.step(x = np.linspace(0, 16, 64), y=np.sort(distances64[-10]), c='red', alpha=.1)
    
    ax.step(np.sort(distances16, axis=1).T, np.linspace(0, 16, 16), c='black', alpha=.1)
    ax.step(np.sort(distances64, axis=1).T, np.linspace(0, 16, 64), c='red', alpha=.1)

    ax.set_yscale('symlog')
    ax.set_ylim(0)
    plt.show()
    # plt.savefig(f'plots/distance_steps/n{n_stars}/cc_dists_{run_id}.png')
    # plt.close()

def hist_nth_star_dists() -> None:
    num_runs = 1000

    distances16 = np.zeros((num_runs, 16))
    distances64 = np.zeros((num_runs, 64))

    for run_id in range(num_runs):
        distances16[run_id] = dists_for_run(run_id, 16)[0]
        distances64[run_id] = dists_for_run(run_id, 64)[0]

    fig, axes = plt.subplots(5,2, sharex=True, layout='constrained')
    histbins = np.logspace(-2, 1.5, 35, base=10)

    for n, (ax16, ax64) in enumerate(axes):
        ax16.set_ylabel(f"star {n}")
        ax16.hist(distances16[:, n], bins=histbins)
        ax64.hist(distances64[:, n], bins=histbins)

    ax16.set_xscale('log')
    axes[0][0].set_title('N=16')
    axes[0][1].set_title('N=64')
    ax16.set_xlabel("Distance to binary COM")
    ax64.set_xlabel("Distance to binary COM")
    fig.suptitle('Radial distributions at start')
    plt.show()

if __name__ == "__main__":
    hist_nth_star_dists()