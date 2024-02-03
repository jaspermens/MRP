import os
from amuse.lab import Particles, read_set_from_file
from tqdm import tqdm

custom_tqdm = lambda x, total: tqdm(x,ascii =" ▏▎▍▌▋▊▉█", ncols=80, total=total)

# " ▖▘▝▗▚▞█"
# " ▏▎▍▌▋▊▉█"

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


def snapshot_at_time(snapshot_directory: str, time: float) -> Particles:
    filename = f'{snapshot_directory}/snapshot_time{float(time):.2f}'
    return read_set_from_file(filename=filename, format='csv')

