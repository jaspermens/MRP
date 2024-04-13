import numpy as np
import matplotlib.pyplot as plt
from plotconfig import * 
from postprocess_fb_properties import get_final_binary_properties_for_n_stars
from plot_end_time import get_end_times_for_n_stars

    
def scatter_coredist_vs_endtime(n_stars: int | list[int] | None = None):
    if isinstance(n_stars, int):
        n_stars = [n_stars]
    
    if n_stars is None:
        n_stars = [12, 14, 16, 64, 72, 80]

    fig, ax = plt.subplots(1,1, figsize=[5,5])

    for n in n_stars:
        coredists, _ = get_final_binary_properties_for_n_stars(n_stars=n)
        endtimes = get_end_times_for_n_stars(n_stars=n)
        print(len(np.argwhere((coredists < 0) & (endtimes > 199))))
        ax.scatter(endtimes, coredists, label=f'N={n}', c=color_for_n(n_stars=n), alpha=.5, s=2)
    
    plt.show()


def plot_coredist_cdfs(n_stars: int | list[int] | None = None):
    if isinstance(n_stars, int):
        n_stars = [n_stars]
    
    if n_stars is None:
        n_stars = [12, 14, 16, 64, 72, 80]
        
    fig, ax = plt.subplots(1,1, figsize=[5,5])

    for n in n_stars:
        coredists, _ = get_final_binary_properties_for_n_stars(n_stars=n)
        ax.ecdf(coredists[coredists > -1], label=f'N={n}', c=color_for_n(n_stars=n))

    ax.legend()
    ax.set_title('CDF of distances to cluster core @ T_cc')
    ax.set_ylim(0.01, 1.01)
    ax.set_xlim(left=0.01)
    ax.set_xlabel('Core distance')

    plt.show()

if __name__ == '__main__':
    # plot_coredist_cdfs()

    scatter_coredist_vs_endtime(n_stars=None)
    # get_final_binary_distance_in_snapshot()
    # get_final_binary_for_run(n_stars=16, run_id=0)
