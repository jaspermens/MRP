from amuse.datamodel import Particles
from amuse.lab import nbody_system
from amuse.lab import read_set_from_file

import matplotlib.pyplot as plt
import numpy as np

import matplotlib.pyplot as plt
from decompose_multiples import replace_a_binary
from helpers import stitch_movie, check_clean_directory, custom_tqdm

RUN_DIRECTORY = 'output/n16'
SNAPSHOT_DIRECTORY = RUN_DIRECTORY + '/snapshots'


def snapshot_at_time(time: float) -> Particles:
    filename = f'{SNAPSHOT_DIRECTORY}/snapshot_time{time:.2f}'
    return read_set_from_file(filename=filename, format='csv')


def get_binaries_to_plot(plummer, initial_ke, min_hardness_kt=1):
    """
    Succesively combine binaries in hierarchical configurations
    return a nested list
    Is this easier from posvel?
    """
    binaries_to_plot = []
    binary_replaced = True
    while binary_replaced:
        binary_replaced = replace_a_binary(plummer, initial_ke=initial_ke, hardness_kt=min_hardness_kt, for_plot=True)
        # (id1, id2), hardness_kt = binary_replaced
        if binary_replaced is not None:
            binaries_to_plot.append(binary_replaced)

    return binaries_to_plot


def pos_all(plummer):
    return (plummer.x.value_in(nbody_system.length), 
            plummer.y.value_in(nbody_system.length),
            plummer.z.value_in(nbody_system.length))


def pos_members(plummer, initial_ke):
    multiples = get_binaries_to_plot(plummer=plummer.copy(), initial_ke=initial_ke, min_hardness_kt=1)
    members = []
    hardnesses = []
    ids = []
    nstars = len(plummer)
    for multiple, hardness in multiples:
        for iden in multiple.astype(int):
            if iden > nstars:
                continue
            ids.append(iden)
            members.append(plummer[iden])
            hardnesses.append(hardness)

    posx = [m.x.value_in(nbody_system.length) for m in members]
    posy = [m.y.value_in(nbody_system.length) for m in members]
    posz = [m.z.value_in(nbody_system.length) for m in members]

    return posx, posy, posz, hardnesses, ids


def plot_at_time(time: float, ax, initial_ke, focus_star_id=None, focus_size=None):
    plummer = snapshot_at_time(time)

    all_x, all_y, all_z = pos_all(plummer)
    ax.scatter(all_x, all_y, all_z, c='grey')
    # for i, (x, y, z) in enumerate(zip(all_x, all_y, all_z)):
        # ax.scatter(x, y, z, c='grey', marker=f'${i}$')

    ax.scatter(0, 0, 0, c='red', alpha=0.5, marker="*")
    
    member_x, member_y, member_z, hardness, _ = pos_members(plummer, initial_ke=initial_ke)
    ax.scatter(member_x, member_y, member_z, c=hardness, cmap='viridis', vmin=0, vmax=10)

    # for x, y, z, hardness, name in zip(*pos_members(plummer, initial_ke=initial_ke)):
        # ax.scatter(x, y, z, c=hardness, cmap='viridis', vmin=0, vmax=20, marker=f"${name}$")
        
    ax.set_aspect('equal')
    plot_title = f'Time {float(time):.2f}'
    if not focus_star_id:
        img_center_x, img_center_y, img_center_z = 0, 0, 0
        plot_title += f' focus: star {focus_star_id}'
    else:
        focus_star = plummer[plummer.id == focus_star_id]
        img_center_x = focus_star.x.value_in(nbody_system.length)
        img_center_y = focus_star.y.value_in(nbody_system.length)
        img_center_z = focus_star.z.value_in(nbody_system.length)

    if not focus_size:
        maxdist = np.percentile(np.abs(np.concatenate([all_x, all_y, all_z])), q=90)
    else:
        maxdist = focus_size

    ax.set_xlim(img_center_x-maxdist, img_center_x + maxdist)
    ax.set_ylim(img_center_y-maxdist, img_center_y + maxdist)
    ax.set_zlim(img_center_z-maxdist, img_center_z + maxdist)

    ax.set_title(plot_title)


def make_buncha_plots(img_directory: str, from_time: float, 
                      to_time: float, cadence: float, 
                      focus_star_id: int, focus_size: float) -> None:
    fig = plt.figure(figsize=[10,10], dpi=100)
    ax = fig.add_subplot(projection='3d')

    initial_plummer = snapshot_at_time(0.00)
    initial_ke = initial_plummer.kinetic_energy()

    for i, t in custom_tqdm(enumerate(np.arange(from_time, to_time, cadence)),
                            total=(to_time-from_time)/cadence):
        fig = plt.figure(figsize=[10,10], dpi=100)
        ax = fig.add_subplot(projection='3d')
        plot_at_time(t, ax=ax, initial_ke=initial_ke, focus_size=focus_size, focus_star_id=focus_star_id)
        plt.savefig(f'{img_directory}/movie-{i:>04}.png')
        plt.close()


def make_movie(from_time=0, 
               to_time=10, 
               cadence=0.01, 
               delete_images: bool = False, 
               movie_filename: str = 'movie-temp',
               focus_star_id: int = None,
               focus_size: float = None,) -> None:
    img_directory = f'{RUN_DIRECTORY}/movie_images'

    check_clean_directory(img_directory)

    make_buncha_plots(img_directory=img_directory,
                      from_time=from_time, to_time=to_time, cadence=cadence,
                      focus_star_id = focus_star_id,
                      focus_size = focus_size)
    
    stitch_movie(img_directory, filename=movie_filename, out_directory=RUN_DIRECTORY, delete_images=delete_images)


if __name__ == '__main__':
    focus_id = 12
    t_start = 98
    t_end = 102
    cadence = .01
    movie_fn = f'star{focus_id}_{t_start}to{t_end}'

    make_movie(from_time=t_start, 
               to_time=t_end, 
               cadence=cadence, 
               delete_images=True, 
               movie_filename=movie_fn, 
               focus_star_id = focus_id, 
               focus_size=.5,
               )

    