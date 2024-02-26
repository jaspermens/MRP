import numpy as np
import matplotlib.pyplot as plt
from plotconfig import * 
from postprocess_fb_properties import get_final_binary_properties_for_n_stars


def plot_coredistances():
    n_stars = 16

    coredists, _ = get_final_binary_properties_for_n_stars(n_stars=n_stars)

    plt.hist(coredists)
    plt.show()


if __name__ == '__main__':
    plot_coredistances()
    # get_final_binary_distance_in_snapshot()
    # get_final_binary_for_run(n_stars=16, run_id=0)
