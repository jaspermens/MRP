import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt
import os

from helpers import stitch_movie, check_clean_directory


RUN_DIRECTORY = 'output/n64'
N_STARS = 64
SNAPSHOT_FREQUENCY = 100

def histogram_suspects(history, from_time=0, to_time=10, save_fn: str = None, title='"Binary membership distribution \nFull run'):
    membership_history = parse_history(history, from_time=from_time, to_time=to_time)
    fig, ax = plt.subplots(1,1, figsize=[10,5], dpi=100)
    ax.hist(membership_history, 
            bins=np.arange(start=0, 
                           stop=N_STARS, 
                           step=1), 
            orientation='vertical')
    
    ax.set_title(title)

    if save_fn:
        plt.savefig(save_fn)
        plt.close()
    else:
        plt.show()


def histogram_suspects_full_history():
    history = read_history()
    t_end, _ = time_ids_from_line(history[-1])
    print(t_end)
    histogram_suspects(from_time=0, to_time=t_end)


def animate_histogram(window_size: float = 5, 
                      window_cadence: float = .1, 
                      t_start: float = 0,
                      t_end: float = None):
    history = read_history()
    if not t_end:
        t_end, _ = time_ids_from_line(history[-1])

    img_directory = f'{RUN_DIRECTORY}/hist_movie'
    check_clean_directory(path = img_directory)

    window_starts = np.arange(0, t_end - window_size, window_cadence)
    for i, t in tqdm(enumerate(window_starts),
                     ascii ="-+$0█", ncols=80, total=len(window_starts)):
        
        plot_title = f"Binary membership distribution \nTime {t}"
        img_filename = f'{img_directory}/movie-{i:>04}.png'
        
        histogram_suspects(history, from_time = t, 
                           to_time = t + window_size, 
                           save_fn = img_filename, 
                           title = plot_title)
    
    stitch_movie(img_directory, filename='test', out_directory=RUN_DIRECTORY, delete_images=True)


def read_history():
    history_fn = f'{RUN_DIRECTORY}/history.txt'
    with open(history_fn, 'r') as file:
        # split into lines and skip the header:
        history = file.read().split('\n')[2:-1] 

    return history


def parse_history(history, from_time=0, to_time=10):
    """
    Reads the lines of the history file 
    and returns which stars are binary members at each time
    
    Assumes the following format for history.txt:
        time: 0.00 - decomposition: []
    """
    # times = np.empty_like(np.linspace(from_time, to_time, round(SNAPSHOT_FREQUENCY*(to_time-from_time))))
    
    start_id = round(from_time * SNAPSHOT_FREQUENCY)
    end_id = round(to_time * SNAPSHOT_FREQUENCY)
    decomps = np.empty(end_id-start_id, dtype=object)
    
    for i, line in enumerate(history[start_id:end_id]):
        _, suspects = time_ids_from_line(line)
        decomps[i] = suspects

    full_array = np.concatenate(decomps)

    return full_array

def time_ids_from_line(line):
    time, decomp = line.split(' - ')
    decomp_raw = decomp.split(':')[1]
    time = float(time.split(':')[1])
    # print(float(time))
    # time = float(line[6:10])
    # decomp_raw = line[30:-3]
    
    decomp_str = decomp_raw
    annoying_chars = list(",[]'")
    for char in annoying_chars:
        decomp_str = decomp_str.replace(char, '')
    
    decomp_array = np.array(decomp_str.split(), dtype=str)
    ids = decomp_array[np.where(['.' not in d for d in decomp_array])].astype(int)

    return time, ids


if __name__ == '__main__':
    animate_histogram(t_start = 50, t_end = 80)