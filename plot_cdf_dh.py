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

def plot_cdf_dh(n_stars = None):
    def get_all_interactions_for_runs(n_stars):
        run_ids = get_run_ids_for_n_stars(n_stars=n_stars)
        # run_ids = range(100)

        filename = f'output/all_dh_n{n_stars}.npy'
        try:
            all_dh = np.load(file=filename)
            return all_dh
        except FileNotFoundError:
            pass

        all_dh = np.zeros((len(run_ids)), dtype=object)
        for i, run_id in enumerate(run_ids):
            _, dh = get_interactions_for_run(n_stars=n_stars, run_id=run_id)
            if len(dh) == 0:
                continue
            all_dh[i] = dh

        all_dh = np.concatenate(all_dh).ravel()

        np.save(file=filename, arr=all_dh)
        return all_dh
    
    def get_interactions_for_run(run_id, n_stars):
        nonlocal DNFS
        
        times, hardnesses = get_hbh_for_run(run_id=run_id, n_stars=n_stars)
        if len(hardnesses) == 0 or hardnesses[0] < 10:
            DNFS.append(f'N{n_stars} run {run_id}')
            return [0], [0]
        
        times = times - times[0]
        dh = np.abs(np.diff(hardnesses))
        times = times[:-1]

        return times, dh
    
    DNFS = []


    if isinstance(n_stars, int):
        n_stars = [n_stars]
    
    if n_stars is None:
        n_stars = [12, 14, 16, 64, 72, 80]
        
    fig, ax = plt.subplots(1,1, figsize=[5,5])

    for n in n_stars:
        dh = get_all_interactions_for_runs(n_stars=n)
        # dh = np.log10(np.ravel(dh[dh>1]))
        dh = np.ravel(dh)
        # dh = np.ravel(dh[(dh < 1) & (dh > 0)])

        ax.ecdf(dh, label=f'N={n}', c=color_for_n(n_stars=n))
    
    ax.legend()
    ax.set_title('CDF dh until T_cc')
    # ax.set_ylim(0.0001, 1.01)
    # ax.set_xlim(1.00,)
    ax.set_xscale('log')
    ax.set_xlabel('dH')

    plt.show()

if __name__ == "__main__":
    # hist_encounter_duration()
    # plot_hardness_over_time()
    plot_cdf_dh()

    # plot_stacked_hardnesses()
