import matplotlib.pyplot as plt
import numpy as np

from helpers import custom_tqdm, get_run_directory, read_history_csv, get_run_ids_for_n_stars
from heritage_from_history import compile_heritage
from itertools import groupby
import scipy

def get_hbh_for_run(run_id: int, n_stars: int) -> tuple[np.array, np.array]:
    # times, _, hardnesses, *_ = read_history_csv(run_id=run_id, n_stars=n_stars)
    times, hardnesses, _ = compile_heritage(n_stars=n_stars, run_id=run_id)
    return times, hardnesses


def split_hardness_curve(run_id, n_stars, boundary):
    times, hardnesses = get_hbh_for_run(run_id=run_id, n_stars=n_stars)
    if len(times) == 0:
        return np.array([[0, 0, 0], [0, 0, 0], [0, 0, 0]])
    # buttord = 1
    ord, wn = scipy.signal.buttord(wp=10, ws=.5, gpass=1, gstop=18, fs=128)
    print(wn)
    b, a = scipy.signal.butter(ord, wn, 'highpass', fs=128)
    fast_signal = scipy.signal.filtfilt(b, a, np.gradient(hardnesses))

    b, a = scipy.signal.butter(ord, wn, 'lowpass', fs=128)
    slow_signal = scipy.signal.filtfilt(b, a, np.gradient(hardnesses))

    return times, fast_signal, slow_signal
    

def get_interaction_widths(n_stars, run_id, boundary):
    times, fast_signal, slow_signal = split_hardness_curve(run_id=run_id, n_stars=n_stars, boundary=boundary)
    
    abs_gradient_fast = np.gradient(fast_signal)
    abs_gradient_slow = np.gradient(slow_signal)
    
    fig, (
        (ax0, ax0a), 
        (ax3, ax4),
        (ax1, ax2), 
          ) = plt.subplots(3,2, sharex=False, figsize=[12,12], layout='constrained')
    # ax0.plot(times, fast_signal + slow_signal, c='grey')
    ax0.plot(*get_hbh_for_run(run_id=run_id, n_stars=n_stars), c='grey')
    ax0a.plot(*get_hbh_for_run(run_id=run_id+1, n_stars=n_stars), c='grey')
    ax1.plot(*get_hbh_for_run(run_id=run_id+2, n_stars=n_stars), c='grey')
    ax2.plot(*get_hbh_for_run(run_id=run_id+3, n_stars=n_stars), c='grey')
    ax3.plot(*get_hbh_for_run(run_id=run_id+4, n_stars=n_stars), c='grey')
    ax4.plot(*get_hbh_for_run(run_id=run_id+5, n_stars=n_stars), c='grey')
    # ax1.plot(times, abs_gradient_fast, c='black')
    # ax2.plot(times, abs_gradient_slow, c='red')    
    # ax3.plot(times, fast_signal, c='black')
    # ax4.plot(times, slow_signal, c='red')
    plt.show()


def plot_split_curve(run_id, n_stars, boundary):
    times, slow_signal, fast_signal = split_hardness_curve(run_id=run_id, n_stars=n_stars, boundary=boundary)
    fig, (ax1, ax2) = plt.subplots(1,2, sharey=True)
    
    ax1.plot(times, slow_signal)
    ax2.plot(times, fast_signal)

    plt.show()



def plot_buncha_split_curves(n_stars: int, boundary: float):
    run_ids = range(100)

    fig, (ax1, ax2) = plt.subplots(1,2, sharey=False)
    
    for rnid in run_ids:
        times16, fast_signal16, slow_signal16 = split_hardness_curve(n_stars=16, run_id=rnid, boundary=boundary)
        times64, fast_signal64, slow_signal64 = split_hardness_curve(n_stars=64, run_id=rnid, boundary=boundary)

        relax_time = 0.1*64/np.log(64)
        ax1.plot(times16 - times16[0], np.abs(np.gradient(slow_signal16)), alpha=.1, c='black')
        ax2.plot(times16 - times16[0], np.abs(np.gradient(fast_signal16)), alpha=.1, c='black')
        ax1.plot((times64 - times64[0]) / relax_time, np.abs(np.gradient(slow_signal64)), alpha=.1, c='red')
        ax2.plot((times64 - times64[0]) / relax_time, np.abs(np.gradient(fast_signal64)), alpha=.1, c='red')

    plt.show()



if __name__ == "__main__":
    # split_hardness_curve(run_id=3, n_stars=16, boundary=1)  
    get_interaction_widths(n_stars=64, run_id=0 , boundary=48)
    # plot_buncha_split_curves(n_stars=64, boundary=1)