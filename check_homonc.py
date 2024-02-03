from amuse.lab import units, constants
from amuse.ic.plummer import new_plummer_model
from amuse.ext.orbital_elements import new_binary_from_orbital_elements as make_binary
from amuse.ext.orbital_elements import generate_binaries
from amuse.datamodel import Particles, Particle
from amuse.lab import nbody_system
from amuse.lab import write_set_to_file, read_set_from_file
from amuse.community.ph4.interface import Ph4

import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm

import matplotlib.pyplot as plt
from amuse.plot import scatter, xlabel, ylabel


def BracketCycler():
    # options = ['[]', '()', '<>', '/|']
    options = ['[]', '[]', '[]', '[]']
    i = 0
    while True:
        yield options[i % len(options)]
        i += 1


def find_composite_multiples(plummer, min_hardness_kt=1):
    """
    Succesively combine binaries in hierarchical configurations
    return a nested list
    Is this easier from posvel?
    """
    binaries_to_plot = []

    def collapse_binary(binary):
        compos = binary.center_of_mass()
        comvel = binary.center_of_mass_velocity()
        commass = binary.mass.sum()
        comparticle = Particles(1)
        comparticle.mass = commass
        comparticle.velocity = comvel
        comparticle.position = compos
        comparticle.id = 1000
        comparticle.color = 'red'

        return comparticle

    def replace_a_binary(plummer, hardness_kt):
        hardness_prefactor =  INITIAL_KE / plummer.kinetic_energy()
        binaries = plummer.get_binaries(G=nbody_system.G, hardness=2/3*hardness_kt * hardness_prefactor)
        if len(binaries) == 0:
            return None
        
        def pick_binary(bins):
            max_id = np.argmax([binary[0].hardness for binary in bins])
            return bins[max_id]
        
        # for target in binaries:
        plummer.color = 'black'
        target = pick_binary(binaries)
        comparticle = collapse_binary(target)
        plummer.add_particle(comparticle)
        target.color = 'blue'
        plummer.remove_particles(target)
        plummer.add_particle(target)
        plot_plummer(plummer)
        plummer.remove_particles(target)
        return target.id, target[0].hardness*3/2

    while True:
        binary_replaced = replace_a_binary(plummer, hardness_kt=min_hardness_kt)
        
        if not binary_replaced:
            return binaries_to_plot

        # (id1, id2), hardness_kt = binary_replaced
        binaries_to_plot.append(binary_replaced)
        
def plot_plummer(plummer):
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')

    x = plummer.x.value_in(nbody_system.length)
    y = plummer.y.value_in(nbody_system.length)
    z = plummer.z.value_in(nbody_system.length)

    ax.scatter(x, y, z, color=plummer.color)
    maxdist = np.percentile(np.abs(np.concatenate([x, y, z])), q=90)
    ax.set_xlim(-maxdist, maxdist)
    ax.set_ylim(-maxdist, maxdist)
    ax.set_zlim(-maxdist, maxdist)

    plt.show()

def plot_multiples_for_plummer(plummer):
    multiples = find_composite_multiples(plummer=plummer.copy(), min_hardness_kt=2)
    scatter(plummer.x, plummer.y, alpha=0.1, color='grey')
    decomp_ids = []
    for multiple in multiples:
        for id in multiple[0]:
            if id == 1000:
                continue
            decomp_ids.append(int(id))

    # decomp_ids = list(set(decomp_ids))
    print(decomp_ids)
    # binarymembers = []
    for id in decomp_ids:
        # binarymembers.append(plummer[int(id)])
        member = plummer[id]
        scatter(member.x, member.y, alpha=1, color='black')

    plt.show()

SNAPSHOT_DIRECTORY = 'snapshots_n16_homonc/'

def snapshot_at_time(time: float) -> Particles:
    filename = SNAPSHOT_DIRECTORY + f'test_snapshots_time{time:.2f}'
    return read_set_from_file(filename=filename, format='csv')

initial_plummer = snapshot_at_time(0.00)
INITIAL_KE = initial_plummer.kinetic_energy()

if __name__ == '__main__':
    
    time = 9.94
    plummer = snapshot_at_time(time)
    find_composite_multiples(plummer, min_hardness_kt=1)
# plummer.id = plummer.id.astype(str)
# plummer.id = [str(identifier) for identifier in plummer.id]
# ids = find_composite_multiples(plummer.copy(), min_hardness_kt=0.1)
    
    # scatter(plummer.x, plummer.y)
    # plt.show()


# if __name__ == '__main__':
    # inspect_interaction_at_time()

    