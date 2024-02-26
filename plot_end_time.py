import numpy as np
import matplotlib.pyplot as plt

from helpers import read_history_csv, custom_tqdm, get_run_ids_for_n_stars, get_output_path
from plotconfig import *


def get_end_time_for_run(n_stars: int, run_id: int):
    times, *_ = read_history_csv(n_stars=n_stars, run_id=run_id)
    return times[-1]


def get_end_times_for_n_stars(n_stars: int):
    npz_filename = f'{get_output_path()}/end_times_n{n_stars}.npy'
    try:
        end_times = np.load(file=npz_filename)
        print(f'File found! re-using the {len(end_times)} end times from {npz_filename}...')
        return end_times
    
    except FileNotFoundError:
        print(f'File {npz_filename} not found- redoing')
        pass

    run_ids = get_run_ids_for_n_stars(n_stars=n_stars)
    
    end_times = np.zeros_like(run_ids, dtype=float)
    for i,run_id in custom_tqdm(enumerate(run_ids), total=len(run_ids)):
        end_times[i] = get_end_time_for_run(n_stars=n_stars, run_id=run_id)

    np.save(file=npz_filename, arr=end_times)
    return end_times


def plot_end_time_cdfs(n_stars: int | list[int]):
    def color_for_n(n_stars: int):
        match n_stars:
            case 16: return 'red'
            case 14: return 'firebrick'
            case 12: return 'maroon'
            case 64: return 'blue'
            case 72: return 'darkblue'
            case 80: return 'navy'
    if isinstance(n_stars, int):
        n_stars = [n_stars]

    fig, ax = plt.subplots(1,1, figsize=[5,5])

    for n in n_stars:
        t_end = get_end_times_for_n_stars(n_stars=n)
        ax.ecdf(t_end, label=f'N={n}', c=color_for_n(n_stars=n))

    ax.legend()
    ax.set_title('CDF of core collapse time')
    ax.set_ylim(0, 1.05)
    ax.set_xlim(0, 201)
    ax.set_xlabel('Time of binary formation')

    plt.show()

if __name__ == '__main__':
    plot_end_time_cdfs(n_stars=np.sort([16, 64, 12, 14, 72, 80]))