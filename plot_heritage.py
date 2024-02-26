import matplotlib.pyplot as plt
import numpy as np

from helpers import custom_tqdm, get_run_directory, read_history_csv, get_run_ids_for_n_stars
from heritage_from_history import compile_heritage
from itertools import groupby
import scipy

def plot_hardness_over_time():
    n_stars = 16
    run_ids = range(100)
    for rnid in custom_tqdm(run_ids, total=len(run_ids)):    
        times, hardnesses, _ = compile_heritage(run_id=rnid, n_stars=n_stars)
        fig, ax = plt.subplots(1,1, figsize=[8, 5])
        ax.plot(times, hardnesses)
        ax.set_title('Hardness over time for the final hard binary\n(and its direct predecessor binaries)')
        ax.set_ylabel('Hardness $(kT_0)$')
        ax.set_xlabel('time [nbody]')
        # ax.set_xlim(left=0)
        directory = get_run_directory(n_stars=n_stars, run_id=rnid)
        plt.savefig(f'{directory}/hardness_history_{rnid:>04}.png')
        plt.close()


def hist_encounter_energy():
    run_ids = range(600)
    times16 = []
    interactions16 = []
    interactions64 = []

    for rnid in custom_tqdm(run_ids, total=len(run_ids)):
        _, h16, _ = compile_heritage(run_id=rnid, n_stars=16)
        _, h64, _ = compile_heritage(run_id=rnid, n_stars=64)
        
        dh16 = np.diff(h16)
        interactions16 = np.append(interactions16, dh16)

        dh64 = np.diff(h64)
        interactions64 = np.append(interactions64, dh64)

    fig, ax = plt.subplots(1, 1)
    ax.set_title("Interaction energy distribution (nonzero only)")
    ax.set_xlabel("Interaction energy [$kT_0$]")
    ax.set_ylabel("Frequency")

    ax.hist(interactions16[np.abs(interactions16) > 0.01], bins=np.arange(-5, 5, .5), density=True, alpha=.7, label='N=16')
    ax.hist(interactions64[np.abs(interactions64) > 0.01], bins=np.arange(-5, 5, .5), density=True, alpha=.7, label='N=64')
    ax.legend()

    plt.show()

def hist_encounter_duration():

    # def plot_history_on_ax(ax, run_id, n_stars):
    #     times, hardnesses, _ = compile_heritage(run_id=run_id, n_stars=n_stars)
    #     times = np.array(times) - times[0]
    #     dh = np.diff(hardnesses)
    #     times = times[1:]
        
    #     ax.plot(times, np.sign(dh),
    #                  lw=1, 
    #                  c='black', 
    #                  alpha=.1,
    #             )
        
    # run_ids = range(60)
    # fig, (ax16, ax64) = plt.subplots(2,1, figsize=[8, 5], sharex=False, sharey=True, layout='constrained')
    # for rnid in custom_tqdm(run_ids, total=len(run_ids)):    
    #     plot_history_on_ax(ax=ax16, n_stars=16, run_id=rnid)
    #     plot_history_on_ax(ax=ax64, n_stars=64, run_id=rnid)


    # plt.show()
    run_ids = range(600)
    times16 = []
    interactions16 = []
    interactions64 = []

    for rnid in custom_tqdm(run_ids, total=len(run_ids)):
        _, h16, _ = compile_heritage(run_id=rnid, n_stars=16)
        _, h64, _ = compile_heritage(run_id=rnid, n_stars=64)
        
        dh16 = np.abs(np.sign(np.diff(h16)))
        peaks16, _ = scipy.signal.find_peaks(dh16)
        widths16 = scipy.signal.peak_widths(dh16, peaks16)
        interactions16 = np.append(widths16, dh16)

        dh64 = np.abs(np.sign(np.diff(h64)))
        peaks64, _ = scipy.signal.find_peaks(dh64)
        widths64 = scipy.signal.peak_widths(dh64, peaks64)
        interactions64 = np.append(widths64, dh64)

    fig, ax = plt.subplots(1, 1)
    ax.set_title("Interaction energy distribution (nonzero only)")
    ax.set_xlabel("Interaction energy [$kT_0$]")
    ax.set_ylabel("Frequency")

    ax.hist(interactions16, density=True, alpha=.7, label='N=16')
    ax.hist(interactions64, density=True, alpha=.7, label='N=64')
    ax.legend()

    plt.show()


def get_hbh_for_run(run_id: int, n_stars: int) -> tuple[np.array, np.array]:
    times, _, hardnesses, *_ = read_history_csv(run_id=run_id, n_stars=n_stars)
    return times, hardnesses


def plot_buncha_hardnesses_over_time():
    def plot_history_on_ax(ax, run_id, n_stars):
        nonlocal DNFS
        # times, hardnesses, _ = compile_heritage(run_id=run_id, n_stars=n_stars)
        times, hardnesses = get_hbh_for_run(run_id=run_id, n_stars=n_stars)

        if len(hardnesses) == 0 or hardnesses[-1] < 10:
            DNFS.append(f'N{n_stars} run {run_id}')
            return
        
        # times = times - times[-1]
        # hardnesses_smoothed = scipy.ndimage.maximum_filter1d(hardnesses, size=50)
        # mask = ((hardnesses < 10) & (hardnesses > 0))
        mask = hardnesses > 0
        dh = np.abs(np.diff(hardnesses))
        # ax.plot(times[mask], hardnesses[mask],
        ax.scatter(times[:-1][dh>0], dh[dh>0],
        # ax.scatter(times[1:][dh>0.1], hardnesses[1:][dh>0.1],
                     lw=1, 
                     c='black', 
                     alpha=.1,
                     s=.1,
                )
        
    run_ids = range(100)
    DNFS = []
    fig, (ax16, ax64) = plt.subplots(2,1, figsize=[8, 15], sharex=True, sharey=True, layout='constrained')
    for rnid in custom_tqdm(run_ids, total=len(run_ids)):    
        plot_history_on_ax(ax=ax16, n_stars=16, run_id=rnid)
        plot_history_on_ax(ax=ax64, n_stars=64, run_id=rnid)

    print(DNFS)
    fig.suptitle('Hardness over time for the final hard binary\n(and its direct predecessor binaries)')
    ax64.set_ylabel('Hardness [$kT_0$]')
    ax16.set_ylabel('Hardness [$kT_0$]')
    ax16.set_title('N=16')
    ax64.set_title('N=64')
    ax64.set_xlabel('time [nbody]')
    ax64.set_yscale('log')
    # ax64.set_xscale('symlog')
    # ax64.set_xlim(right=.05)
    # ax64.set_ylim(-1, 25)
    plt.show()



def plot_stacked_hardnesses():
    def get_all_interactions_for_runs(n_stars):
        # run_ids = get_run_ids_for_n_stars(n_stars=n_stars)
        run_ids = range(100)

        filename = f'output/all_dh_n{n_stars}.npy'
        try:
            all_dh = np.load(file=filename)
            return all_dh
        except FileNotFoundError:
            pass

        all_dh = np.zeros((len(run_ids), 128*200 - 1))
        for i, run_id in enumerate(run_ids):
            _, dh = get_interactions_for_run(n_stars=n_stars, run_id=run_id)
            if len(dh) == 0:
                continue
            all_dh[i] = np.pad(dh, [0, all_dh.shape[1]-len(dh)])

        np.save(file=filename, arr=all_dh)
        return all_dh
    
    def get_interactions_for_run(run_id, n_stars):
        nonlocal DNFS
        
        times, hardnesses = get_hbh_for_run(run_id=run_id, n_stars=n_stars)

        if len(hardnesses) == 0 or hardnesses[-1] < 10:
            DNFS.append(f'N{n_stars} run {run_id}')
            return [], []
        
        times = times - times[-1]
        dh = np.abs(np.diff(hardnesses))
        times = times[:-1]

        return times, dh
    
    def plot_history_on_ax(ax, n_stars: int):
        all_dh = get_all_interactions_for_runs(n_stars=n_stars)
        all_times = np.ones_like(all_dh)
        all_times[:] = np.arange(0, 200 - 1/128, 1/128)

        hist2d, *_ = np.histogram2d(all_times.ravel()[all_dh.ravel() > 0], np.log10(all_dh.ravel()[all_dh.ravel() > 0]), bins=[100,10], density=True)
        # print(hist2d)
        # ax.imshow(hist2d.T)
        # all_dh = np.digitize(all_dh, bins=[0, 1, 2])
        # print(all_dh)
        # print(np.unique(all_dh, return_counts=True, axis=1))
        ax.stackplot(np.linspace(0, 200, 100), hist2d.T[::-1])
        # ax.stackplot(all_times[0].astype(int), all_dh)
        
    run_ids = range(10)
    DNFS = []
    fig, (ax16, ax64) = plt.subplots(2,1, figsize=[8, 15], sharex=True, sharey=True, layout='constrained')
    plot_history_on_ax(ax=ax16, n_stars=16)
    plot_history_on_ax(ax=ax64, n_stars=64)

    print(DNFS)
    fig.suptitle('Hardness over time for the final hard binary\n(and its direct predecessor binaries)')
    ax64.set_ylabel('Hardness [$kT_0$]')
    ax16.set_ylabel('Hardness [$kT_0$]')
    ax16.set_title('N=16')
    ax64.set_title('N=64')
    ax64.set_xlabel('time [nbody]')
    # ax64.set_xscale('symlog')
    # ax64.set_xlim(right=.05)
    # ax64.set_ylim(-1, 25)
    plt.show()

if __name__ == "__main__":
    # hist_encounter_duration()
    # plot_hardness_over_time()
    # plot_buncha_hardnesses_over_time()
    plot_stacked_hardnesses()
