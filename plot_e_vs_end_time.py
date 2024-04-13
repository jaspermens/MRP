import numpy as np
import matplotlib.pyplot as plt
from plotconfig import * 
from postprocess_fb_properties import get_final_binary_properties_for_n_stars
from plot_end_time import get_end_times_for_n_stars

from scipy.stats import kstest


def scatter_r_vs_e(n_stars, ax):
    _, e = get_final_binary_properties_for_n_stars(n_stars=n_stars)
    t_end = get_end_times_for_n_stars(n_stars=n_stars)
    mask = (e > -1) 
    # [slice(None,None,10)]
    ax.scatter(e[mask], t_end[mask], c=color_for_n(n_stars=n_stars), alpha=0.3, s=4)

def plot_coredist_cdfs(n_stars: int | list[int] | None = None):
    if isinstance(n_stars, int):
        n_stars = [n_stars]
    
    if n_stars is None:
        n_stars = [12, 14, 16, 64, 72, 80]
        
    fig, axes = plt.subplots(1,len(n_stars), figsize=[4*len(n_stars),4], sharey=True, layout='constrained')

    for n, ax in zip(n_stars, axes):
        # put_cdf_with_errorbars_on_ax(n_stars=n, ax=ax)
        scatter_r_vs_e(n_stars=n, ax=ax)


    plt.show()

if __name__ == '__main__':
    # compute_kstest(n1=16, n2=64)
    plot_coredist_cdfs(n_stars=[16, 64])
    # get_final_binary_distance_in_snapshot()
    # get_final_binary_for_run(n_stars=16, run_id=0)
