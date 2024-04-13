import numpy as np
import matplotlib.pyplot as plt
from plotconfig import * 
from postprocess_fb_properties import get_final_binary_properties_for_n_stars

from scipy.stats import kstest

def compute_kstest(n1: int, n2: int) -> None:
    _, samples1 = get_final_binary_properties_for_n_stars(n_stars=n1)
    _, samples2 = get_final_binary_properties_for_n_stars(n_stars=n2)
    print(kstest(samples1, samples2))


def put_cdf_with_errorbars_on_ax(n_stars, ax):
    from sksurv.nonparametric import kaplan_meier_estimator

    _, eccentricities = get_final_binary_properties_for_n_stars(n_stars=n_stars)
    ecc_clean = eccentricities[eccentricities > 0] 
    time, survival_prob, conf_int = kaplan_meier_estimator(
        # np.ones_like(eccentricities, dtype=bool), eccentricities, conf_type="log-log"
        np.ones_like(ecc_clean, dtype=bool), ecc_clean, conf_type="log-log"
        )

    plotcolor = color_for_n(n_stars)
    ax.step(time, 1 - survival_prob, where="post", color=plotcolor, label=f"N={n_stars}")
    ax.fill_between(time, 1 - conf_int[1], 1 - conf_int[0], alpha=0.1, step="post", color=plotcolor)


def put_cdf_on_ax(n_stars, ax):
    from sksurv.nonparametric import kaplan_meier_estimator

    _, eccentricities = get_final_binary_properties_for_n_stars(n_stars=n_stars)
    ecc_clean = eccentricities[eccentricities > 0] 

    plotcolor = color_for_n(n_stars)
    ax.ecdf(ecc_clean, color=plotcolor, label=f"N={n_stars}")


def plot_coredist_cdfs(n_stars: int | list[int] | None = None):
    if isinstance(n_stars, int):
        n_stars = [n_stars]
    
    if n_stars is None:
        n_stars = [12, 14, 16, 64, 72, 80]
        
    fig, ax = plt.subplots(1,1, figsize=[5,5])

    for n in n_stars:
        # put_cdf_with_errorbars_on_ax(n_stars=n, ax=ax)
        put_cdf_on_ax(n_stars=n, ax=ax)

    ax.plot(x:= np.arange(0.01,1, .0001), x**2, c='black')
    ax.legend()
    ax.set_title('CDF of hard binary eccentricity @ T_cc')
    ax.set_ylim(0.01, 1.01)
    ax.set_xlim(left=0.01)
    ax.set_xlabel('e')

    plt.show()

if __name__ == '__main__':
    # compute_kstest(n1=16, n2=64)
    plot_coredist_cdfs()
    # get_final_binary_distance_in_snapshot()
    # get_final_binary_for_run(n_stars=16, run_id=0)
