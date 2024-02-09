import matplotlib.pyplot as plt
import numpy as np

from helpers import custom_tqdm, get_run_directory
from heritage_from_history import compile_heritage

def plot_hardness_over_time():
    n_stars = 16
    run_ids = range(20)
    for rnid in custom_tqdm(run_ids, total=len(run_ids)):    
        times, hardnesses, _ = compile_heritage(run_id=rnid)
        fig, ax = plt.subplots(1,1, figsize=[8, 5])
        ax.plot(times, hardnesses)
        ax.set_title('Hardness over time for the final hard binary\n(and its direct predecessor binaries)')
        ax.set_ylabel('Hardness $(kT_0)$')
        ax.set_xlabel('time [nbody]')
        # ax.set_xlim(left=0)
        directory = get_run_directory(n_stars=n_stars, run_id=rnid)
        plt.savefig(f'{directory}/hardness_history_{rnid:>04}.png')
        plt.close()

    
def plot_buncha_hardnesses_over_time():
    run_ids = range(80)
    fig, (ax16, ax64) = plt.subplots(2,1, figsize=[8, 5], sharex=True, sharey=True, layout='constrained')

    for rnid in custom_tqdm(run_ids, total=len(run_ids)):    
        times16, hardnesses16, _ = compile_heritage(run_id=rnid, n_stars=16)
        times16 = np.array(times16) - times16[0]
        ax16.plot(times16, hardnesses16, lw=1, c='black', alpha=.2)        
        
        times64, hardnesses64, _ = compile_heritage(run_id=rnid, n_stars=64)
        times64 = np.array(times64) - times64[0]
        ax64.plot(times64, hardnesses64, lw=1, c='black', alpha=.2)

    fig.suptitle('Hardness over time for the final hard binary\n(and its direct predecessor binaries)')
    ax64.set_ylabel('Hardness $(kT_0)$')
    ax16.set_ylabel('Hardness $(kT_0)$')
    ax16.set_title('N=16')
    ax64.set_title('N=64')
    ax64.set_xlabel('time [nbody]')
    ax64.set_xscale('symlog')
    ax64.set_xlim(right=.05)
    ax64.set_ylim(top=20)
    plt.show()


if __name__ == "__main__":
    plot_buncha_hardnesses_over_time()