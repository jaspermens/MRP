import matplotlib.pyplot as plt
import numpy as np

from helpers import custom_tqdm, get_run_directory
from heritage_from_history import compile_heritage
from plotconfig import *


def count_stars_in_heritage(heritage: np.array):
    """counts the number of unique star ids involved in any previous binary from this heritage"""
    # heritage = np.array(heritage).flatten()
    # unique_stars = list(set(heritage))
    # return len(unique_stars)
    unique = np.unique(heritage)
    return len(unique)

def count_switches_in_heritage(heritage: np.array):
    """counts how many times the binary composition has changed along the way"""
    swaps = 0
    current_binary = heritage[0]
    for binary in heritage:
        if not any(current_binary) in binary:
            swaps += 1

    return swaps


def hist_family_size_for_n_stars(n_stars: int = 16, n_runs: int = 600):        
    run_ids = range(n_runs)
    family_sizes = []
    for rnid in custom_tqdm(run_ids, total=len(run_ids)):
        _, _, heritage = compile_heritage(run_id=rnid, n_stars=16)
        family_sizes.append(count_stars_in_heritage(heritage=np.array(heritage)))
    plt.hist(family_sizes)
    plt.show()

def hist_family_size(n_runs: int = 10):
    fig, ax = plt.subplots(1, 1)
    
    family_size_16 = []
    family_size_64 = []
    for rnid in custom_tqdm(range(n_runs), total=n_runs):
        _, _, heritage16 = compile_heritage(run_id=rnid, n_stars=16)
        _, _, heritage64 = compile_heritage(run_id=rnid, n_stars=64)

        family_size_16.append(count_stars_in_heritage(heritage=np.array(heritage16)))
        family_size_64.append(count_stars_in_heritage(heritage=np.array(heritage64)))

    print(np.min(family_size_16))
    print(np.min(family_size_64))

    ax.hist(family_size_16, bins=np.arange(0, 15, 1), label='N=16', alpha=0.8, align='left')
    ax.hist(family_size_64, bins=np.arange(0, 15, 1), label='N=64', alpha=0.8, align='left')
    ax.legend()
    ax.set_xlabel(r'family size')
    ax.set_title('Number of stars involved in direct predecessor binaries')
    # fig.suptitle('Histogram of family size')
    # plt.hist(end_times, bins=20)
    plt.show()


if __name__ == "__main__":
    hist_family_size(n_runs = 1000)