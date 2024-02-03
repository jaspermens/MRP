import numpy as np
import matplotlib.pyplot as plt

from amuse.lab import nbody_system
from amuse.lab import Particles

from decompose_multiples import replace_a_binary
from helpers import snapshot_at_time, check_clean_directory, stitch_movie, custom_tqdm

RUN_DIRECTORY = 'output/n64'
SNAPSHOT_DIRECTORY = f'{RUN_DIRECTORY}/snapshots'

INITIAL_KE = snapshot_at_time(SNAPSHOT_DIRECTORY, time=0).kinetic_energy()


def replace_hardest_binary(plummer: Particles) -> Particles:
    _ = replace_a_binary(plummer=plummer, initial_ke=INITIAL_KE, hardness_kt=1, for_plot=True)


def tree_lengths_at_time(time: float) -> list:
    plummer = snapshot_at_time(SNAPSHOT_DIRECTORY, time=time)
    replace_hardest_binary(plummer)
    full_mstlength = plummer.minimum_spanning_tree_length().value_in(nbody_system.length)
    comparticle = plummer[plummer.id == 10000]
    lengths = [0]
    build_up_set = Particles()
    build_up_set.add_particle(comparticle)
    tear_down_set = plummer.copy()
    tear_down_set.remove_particle(comparticle)

    for i in range(len(plummer)):
        next_nn = build_up_set.nearest_neighbour(neighbours=tear_down_set)
        build_up_set.add_particle(next_nn)
        tear_down_set.remove_particle(next_nn)
        mstlength = build_up_set.minimum_spanning_tree_length()
        lengths.append(mstlength.value_in(nbody_system.length))
        if lengths[-2] - lengths[-1] > full_mstlength/20:
            break
        print(f'{i} - {lengths}')

    return lengths        
    # return hb.minimum_spanning_tree_length()

def distances_at_time(time: float) -> list:
    plummer = snapshot_at_time(SNAPSHOT_DIRECTORY, time=time)
    replace_hardest_binary(plummer)
    
    comparticle = plummer[plummer.id == 10000]
    distances = plummer.distances_squared(other_particles=comparticle).flatten().sqrt().value_in(nbody_system.length)

    return np.sort(distances)

def plot_tree_length_at_time(time: float):
    tree_lengths = distances_at_time(time)
    plt.hist(tree_lengths, cumulative=True, bins=np.arange(63, step=1))
    plt.title(f'Time {time:.2f}')


def animate_tree_length_histograms(from_time, to_time, cadence):
    img_directory = f'{RUN_DIRECTORY}/distances'
    check_clean_directory(path = img_directory)
    times = np.arange(from_time, to_time, cadence)
    for i, t in custom_tqdm(enumerate(times), total=len(times)):
        plot_tree_length_at_time(time=t)
        plt.savefig(f'{RUN_DIRECTORY}/distances/movie-{i:>04}.png')
        plt.close()

    stitch_movie(image_directory=img_directory, 
                 filename='distance_test', 
                 out_directory=RUN_DIRECTORY, 
                 delete_images=True)


if __name__ == '__main__':
    animate_tree_length_histograms(from_time=80, to_time=100, cadence=0.01)
    # plot_tree_length_at_time(time=50)