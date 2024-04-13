import numpy as np
import matplotlib.pyplot as plt

from helpers import read_history_csv, custom_tqdm, get_run_ids_for_n_stars, get_output_path
from plotconfig import *

from sksurv.nonparametric import kaplan_meier_estimator

import scipy


def fit_lognormal(end_times):
    x_axis = np.unique(end_times, return_counts=False)

    fit_params = scipy.stats.lognorm.fit(end_times)
    y_axis = scipy.stats.lognorm.cdf(x_axis, *fit_params)

    fig, ax = plt.subplots(1,1)
    ax.ecdf(end_times)
    ax.plot(x_axis, y_axis)
    plt.show()


def plot_lognorm_fits():
    def plot_fit_on_ax(n_stars, ax):
        end_times = get_end_times_for_n_stars(n_stars=n_stars)
        end_times_data = scipy.stats.CensoredData.right_censored(end_times, (end_times > 199.5))
        x_axis = np.unique(end_times, return_counts=False)

        fit_params = scipy.stats.lognorm.fit(end_times_data)
        y_axis = scipy.stats.lognorm.cdf(x_axis, *fit_params) #* len(end_times[end_times < 199])/len(end_times)

        ax.plot(x_axis, y_axis, c=color_for_n(n_stars), zorder=100)#, label=f'N={n_stars}')
        print(n_stars)
        print(fit_params)


    lowns = [12, 14, 16, 32, 128]
    highns = [64, 72, 80, 128, 32]

    fig, (axlow, axhigh) = plt.subplots(1,2, layout='constrained')

    for lown, highn in zip(lowns, highns):
        plot_fit_on_ax(n_stars=lown, ax=axlow)
        plot_fit_on_ax(n_stars=highn, ax=axhigh)

        put_cdf_with_errorbars_on_ax(n_stars=lown, ax=axlow)
        put_cdf_with_errorbars_on_ax(n_stars=highn, ax=axhigh)

    axlow.legend()
    axhigh.legend()
    axlow.set_title('low N')
    axhigh.set_title('high N')
    axlow.set_xlabel('CC time')
    axhigh.set_xlabel('CC time')
    plt.show()        

def compute_kstest(n1: int, n2: int) -> None:
    samples1 = get_end_times_for_n_stars(n_stars=n1)
    samples2 = get_end_times_for_n_stars(n_stars=n2)
    print(scipy.stats.kstest(samples1, samples2))




def get_end_time_for_run(n_stars: int, run_id: int):
    times, *_ = read_history_csv(n_stars=n_stars, run_id=run_id)
    return times[-1]

from helpers import get_final_snapshot_filename
from amuse.lab import read_set_from_file, nbody_system
def get_end_time_from_final_snapshot(n_stars: int, run_id: int):
    snapshot_fn = get_final_snapshot_filename(n_stars=n_stars, run_id=run_id)
    plummer = read_set_from_file(snapshot_fn, format='amuse', copy_history=False, close_file=True)

    end_time = plummer.get_timestamp().value_in(nbody_system.time)
    return end_time


def get_end_times_for_n_stars(n_stars: int):
    npz_filename = f'{get_output_path()}/end_times_n{n_stars}.npy'
    try:
        end_times = np.load(file=npz_filename)
        return end_times # - 9/128 # 10/128 patience, but some systems init w/ 10kT and we don't want div by 0
    
    except FileNotFoundError:
        print(f'File {npz_filename} not found- redoing')
        pass

    run_ids = get_run_ids_for_n_stars(n_stars=n_stars)
    
    end_times = np.zeros_like(run_ids, dtype=float)
    for i,run_id in custom_tqdm(enumerate(run_ids), total=len(run_ids)):
        # end_times[i] = get_end_time_for_run(n_stars=n_stars, run_id=run_id)
        end_times[i] = get_end_time_from_final_snapshot(n_stars=n_stars, run_id=run_id)

    np.save(file=npz_filename, arr=end_times)
    return end_times


def plot_end_time_without_errorbars(n_stars: int | list[int]):
    if isinstance(n_stars, int):
        n_stars = [n_stars]

    fig, ax = plt.subplots(1,1, figsize=[5,5])

    for n in n_stars:
        t_end = get_end_times_for_n_stars(n_stars=n)
        ax.ecdf(t_end, label=f'N={n}', c=color_for_n(n_stars=n))

    ax.legend()
    ax.set_title('CDF of core collapse time')
    ax.set_ylim(0.01, 1.05)
    ax.set_xlim(0.01, 201)
    ax.set_xlabel('Time of binary formation')

    plt.show()



def put_cdf_with_errorbars_on_ax(n_stars, ax):
    end_times = get_end_times_for_n_stars(n_stars=n_stars)
    time, survival_prob, conf_int = kaplan_meier_estimator(
        np.ones_like(end_times, dtype=bool), end_times, conf_type="log-log"
        )

    plotcolor = color_for_n(n_stars)
    ax.step(time, 1 - survival_prob, where="post", color=plotcolor, label=f"N={n_stars}")
    ax.fill_between(time, 1 - conf_int[1], 1 - conf_int[0], alpha=0.1, step="post", color=plotcolor)


def plot_end_time_with_errorbars():
    # ns_stars = [12, 14, 16, 32, 64, 72, 80]
    ns_stars = [12, 14, 16, 32, 64, 72, 80, 128]

    fig, ax = plt.subplots(1, 1)
    for n_stars in ns_stars:
        put_cdf_with_errorbars_on_ax(n_stars=n_stars, ax=ax)

    ax.legend()
    ax.set_title('CDF of core collapse time for different N')
    ax.set_ylim(0.01, 1.05)
    ax.set_xlim(0.01, 201)
    ax.set_xlabel('Time of hard binary formation [N-body units]')
    plt.show()


def surface_plot():
    from matplotlib import cm
    import colorcet 
    ns_stars = [12, 14, 16, 32, 64, 72, 80, 128]
    sample_points = np.arange(200, step=.1)
    # sample_points = np.logspace(start=-3, stop=np.log10(200), num=1000, base=10)

    def cdf_for_n_stars(n_stars):
        end_times = get_end_times_for_n_stars(n_stars=n_stars)
        end_times_data = scipy.stats.CensoredData.right_censored(end_times, (end_times >199.5))
        cdf = scipy.stats.ecdf(end_times_data).cdf.evaluate(sample_points)
        return cdf
    
    end_time_cdfs = np.array([cdf_for_n_stars(n) for n in ns_stars])

    interpolator = scipy.interpolate.interp1d(x=np.array(ns_stars), y=end_time_cdfs, axis=0)
    interp_points = np.logspace(np.log2(12), np.log2(128), num=100, base=2)
    interp_points = np.linspace(12, 128, num=1000)
    end_times_interpd = interpolator(interp_points)
    print(end_times_interpd.shape)
    fig, ax = plt.subplots(1,1, subplot_kw={"projection":'3d'})
    # X, Y = np.meshgrid(sample_points, range(len(ns_stars)))
    X, Y = np.meshgrid(np.log10(sample_points), interp_points)
    ax.plot_surface(X, Y, end_times_interpd, cmap='cet_bmy', antialiased=False,
                     rcount=100, ccount=150,
                    #  vmin=0.6, 
                    #  vmax=.9,
                       )
    # ax.set_yticks(ticks = np.log2(ns_stars), labels=ns_stars)
    # ax.set_yticks(ticks = range(len(ns_stars)), labels=ns_stars)

    plt.show()


def plot_pdfs_stacked():
    ns_stars = [12, 14, 16, 32, 64, 72, 80, 128]
    sample_points = np.arange(200, step=.1)

    def pdf_for_n_stars(n_stars):
        end_times = get_end_times_for_n_stars(n_stars=n_stars)
        print(f"{n_stars = }, {np.min(end_times) = }")
        end_times_data = scipy.stats.CensoredData.right_censored(end_times, (end_times >199.5))

        fit_params = scipy.stats.lognorm.fit(end_times_data)
        pdf = scipy.stats.lognorm.pdf(sample_points, *fit_params)
        return pdf
    
    # fig, axes = plt.subplots(len(ns_stars), 1, sharex=True)
    fig, ax = plt.subplots(1, 1, sharex=True)
    plt.subplots_adjust(hspace=-.7)

    # for n_stars, ax in zip(ns_stars, axes):
    for n_stars in ns_stars:
        ax.plot(sample_points, pdf_for_n_stars(n_stars=n_stars), color=color_for_n(n_stars=n_stars))
    ax.set_xscale('log')
    # ax.set_axis_off()
    plt.show()

if __name__ == '__main__':
    # fit_lognormal(end_times=get_end_times_for_n_stars(n_stars=16))
    # plot_lognorm_fits()
    # compute_kstest(n1=14, n2=16)
    # surface_plot()
    # plot_pdfs_stacked()
    plot_end_time_with_errorbars()
    # plot_end_time_cdfs(n_stars=np.sort([16, 64, 12, 14, 72, 80]))