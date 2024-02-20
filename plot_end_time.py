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
    run_ids = range(600)
    fig, (ax16, ax64) = plt.subplots(2, 1, sharex=True)
    end_times_16 = []
    end_times_64 = []
    for run_id in run_ids:
        end_times_16.append(get_end_time_for_run(n_stars=16, run_id=run_id))
        end_times_64.append(get_end_time_for_run(n_stars=64, run_id=run_id))

    ax16.hist(end_times_16, bins=np.arange(0, 200, 10))
    ax64.hist(end_times_64, bins=np.arange(0, 200, 10))
    ax64.set_xlabel(r'$T_{cc}$')

    ax16.set_title('N=16')
    ax64.set_title('N=64')
    fig.suptitle('Histogram of simulation end times')
    # plt.hist(end_times, bins=20)
    plt.show()

if __name__ == '__main__':
    # get_end_time_for_run(n_stars=16, run_id=0)
    plot_end_histogram()