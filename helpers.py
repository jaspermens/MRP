import os
from amuse.lab import Particles, read_set_from_file
from tqdm import tqdm
import numpy as np


def running_on_alice() -> bool:
    return os.path.exists('/home/s2015242/data1/output')


def get_run_directory(n_stars: int, run_id: int) -> str:
    end = f'n{n_stars}/run_{run_id:>04}'
    return f'{get_output_path()}/{end}'


def get_output_path() -> str: 
    if running_on_alice():
        return f'/home/s2015242/data1/output'
    else:
        return f'/home/jasper/studie/MRP/output'
    

def get_final_snapshot_filename(n_stars: int, run_id: int):
    directory = get_run_directory(n_stars=n_stars, run_id=run_id)
    filename = f'{directory}/final_snapshot_n{n_stars}_run{run_id}.amuse'
    return filename


def stitch_movie(image_directory: str, filename: str, out_directory: str, delete_images=False) -> None:
    command = f"ffmpeg -start_number 0 -i {image_directory}/movie-%04d.png -c:v libx264 -vb 20M -r 20 -pix_fmt yuv420p -filter:v 'setpts=2*PTS' -y {out_directory}/{filename}.mp4"
    print(command)
    # command = "ffmpeg -framerate 5 -pattern_type glob -i '{current_path}/output-figures/*.png' -filter:v 'setpts=2*PTS'-vb 20M -c:v libx264 -pix_fmt yuv420p -y {outdir}/movie-{filename}.mp4"

    os.system(command)

    if delete_images:
        os.system(f"rm -rf {image_directory}")


def check_clean_directory(path: str):
    if os.path.exists(path):
        print(f'Directory {path} already exists!')
        print("Deleting contents and starting fresh...")
        os.system(f"rm -rf {path}")
        
    os.mkdir(path)


def snapshot_at_time_depr(snapshot_directory: str, time: float) -> Particles:
    filename = f'{snapshot_directory}/snapshot_time{float(time):.2f}'
    return read_set_from_file(filename=filename, format='csv')


def snapshot_at_time(snapshots, time: float, snapshot_frequency: int = 128) -> Particles:
    target = snapshots[round(time*snapshot_frequency)]
    return target


def read_snapshot_file(run_directory: str):
    all_snapshots = read_set_from_file(f"{run_directory}/snapshots.log", format='amuse', copy_history=False, close_file=False)
    snapshots = np.array([snapshot for snapshot in all_snapshots.history])
    return snapshots

def get_final_snapshot_depr(run_id: int, n_stars: int) -> Particles: 
    run_directory = get_run_directory(n_stars=n_stars, run_id=run_id)
    all_snapshots = read_set_from_file(f"{run_directory}/snapshots.log", format='amuse', copy_history=False, close_file=True)
    return all_snapshots

def get_final_snapshot(run_id: int, n_stars: int) -> Particles:
    final_snapshot_fn = get_final_snapshot_filename(n_stars=n_stars, run_id=run_id)
    final_snapshot = read_set_from_file(final_snapshot_fn, format='amuse', copy_history=False, close_file=True)
    return final_snapshot

def read_history_csv(run_id: int, n_stars: int) -> np.array:
    run_directory = get_run_directory(n_stars=n_stars, run_id=run_id)
    filename = f'{run_directory}/decompositions.txt'
    times, ke_history, hardest_hardness, decomp = np.genfromtxt(filename, dtype=str, delimiter=' - ', skip_header=1, unpack=True)
    return times.astype(float), ke_history.astype(float), hardest_hardness.astype(float), decomp


def get_run_ids_for_n_stars(n_stars: int):
    if n_stars in [16, 32, 64]:
        run_ids = range(1000)
    elif n_stars in [72, 80, 12, 14, 128]:
        run_ids = range(500)
    else:
        raise NotImplementedError
    
    return run_ids


if running_on_alice():
    custom_tqdm = lambda x, total: x
else:
    custom_tqdm = lambda x, total: tqdm(x,ascii =" ▏▎▍▌▋▊▉█", ncols=80, total=total)

# " ▖▘▝▗▚▞█"
# " ▏▎▍▌▋▊▉█"



def test_final_snapshot_speed():
    run_id = 1
    n_stars = 16
    snapshot = read_snapshot_file(run_directory=get_run_directory(n_stars=n_stars, run_id=run_id))[-1]
    final_time = snapshot.get_timestamp()

    for _ in custom_tqdm(range(1), total=10):
        snapshot = get_final_snapshot(run_id=run_id, n_stars=n_stars)
        assert snapshot.get_timestamp() == final_time

    for _ in custom_tqdm(range(1), total=10):
        snapshot = read_snapshot_file(run_directory=get_run_directory(n_stars=n_stars, run_id=run_id))[-1]
        assert snapshot.get_timestamp() == final_time


if __name__ == '__main__':
    test_final_snapshot_speed()
    # rundir = get_run_directory(n_stars=16, run_id=1)
    # snaps = read_snapshot_file(run_directory=rundir)