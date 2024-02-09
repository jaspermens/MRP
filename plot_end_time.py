import numpy as np
import matplotlib.pyplot as plt

from helpers import get_run_directory


def get_end_time_for_run(n_stars: int, run_id: int):
    directory = get_run_directory(n_stars=n_stars, run_id=run_id)
    history_fn = f'{directory}/history.txt'
    with open(history_fn, 'r') as file:
        end_time = file.read().split('\n')[-1]
        print(end_time)
        return float(end_time)


def plot_end_histogram():
    run_ids = range(120)
    n_stars = 64
    end_times = []
    for run_id in run_ids:
        end_time = get_end_time_for_run(n_stars=n_stars, run_id=run_id)
        end_times.append(end_time)

    plt.hist(end_times, bins=20)
    plt.show()

if __name__ == '__main__':
    # get_end_time_for_run(n_stars=16, run_id=0)
    plot_end_histogram()