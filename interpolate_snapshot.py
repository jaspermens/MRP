from helpers import snapshot_at_time, custom_tqdm, check_clean_directory
import numpy as np
from amuse.community.ph4.interface import Ph4
from amuse.lab import nbody_system, Particles


SNAPSHOT_DIRECTORY = 'output/n64/snapshots'

def evolve_until_next(time: float):
    plummer_current = snapshot_at_time(snapshot_directory=SNAPSHOT_DIRECTORY, time=time)
    gravity = Ph4()
    gravity.particles.add_particles(plummer_current)
    gravity.evolve_model(.01 | nbody_system.time)
    return gravity.particles

def compare_interpolated_result(time: float):
    plummer_next = snapshot_at_time(snapshot_directory=SNAPSHOT_DIRECTORY, time=time+.01)
    plummer_interpd = evolve_until_next(time=time)
    print(np.max((plummer_next.position - plummer_interpd.position).value_in(nbody_system.length)))


if __name__ == '__main__':
    compare_interpolated_result(time=93.28)