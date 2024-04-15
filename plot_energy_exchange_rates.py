from plotconfig import *
import numpy as np
import matplotlib.pyplot as plt
from helpers import get_run_directory, dists_for_run, edots_for_run, hardnesses_for_run
import scipy



def hist_energy_exchange_fractions():
    x_axis = np.linspace(1, 5, 100)

    def get_contribution_fractions(n_stars: int):
        n_runs = 1
        all_contribution_fractions = np.zeros((1, 100))
        for run_id in range(n_runs):
            edots = edots_for_run(run_id=10, n_stars=n_stars)
            cdf = scipy.stats.ecdf(np.abs(edots).flatten()).cdf.evaluate(x_axis)#*np.linspace(.01, 20.01, n_stars)
            all_contribution_fractions[run_id] = cdf

        return all_contribution_fractions
    
    ns_stars = [16, 64]
    fig, axes = plt.subplots(1, len(ns_stars), layout="constrained", sharey=True)
    for n_stars, ax in zip(ns_stars, axes):
        cfs = get_contribution_fractions(n_stars=n_stars)
        ax.step(x_axis, cfs.T, alpha=1, c='black')
        ax.set_title(f'N={n_stars}')
        ax.set_xlabel(r"$\dot{E}$")
        # ax.set_ylabel("")
    
    # axes[1].set_xlim(0, 16)
    plt.show()

def compare_hardness_hist_and_edots():
    n_stars = 16
    run_id = 10

    edots = edots_for_run(run_id=run_id, n_stars=n_stars)
    hardnesses = hardnesses_for_run(run_id=run_id, n_stars=n_stars)

    # plt.plot(np.diff(hardnesses))
    plt.plot(edots, label=[f"{c}" for c in range(len(edots[0]))])
    # plt.plot(hardnesses)
    plt.legend()
    plt.show()

if __name__ == "__main__":
    # compare_hardness_hist_and_edots()
    hist_energy_exchange_fractions()