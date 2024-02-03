from amuse.ic.plummer import new_plummer_model
from amuse.lab import nbody_system
from amuse.lab import write_set_to_file
from amuse.community.ph4.interface import Ph4

from helpers import custom_tqdm, check_clean_directory
import numpy as np

RNGSEED = 1234
# RNGSEED = 123585
# RNGSEED = 1235
# RNGSEED = 12345

RNG = np.random.default_rng(seed=RNGSEED)
N_BODIES = 16
MIN_HARDNESS_KT = 1
OUTPUT_DIRECTORY = f'output/n{N_BODIES}'

from decompose_multiples import find_composite_multiples

def build_history(t_end: float = 20, 
                  step_size: float = 0.01, 
                  output_directory: str = 'output/test') -> None:
    """
    Performs a simulation and outputs the decomposition at each timestep
    (also stores snapshots)
    """
    gravity = Ph4()

    # initialise a plummer sphere using a set seed
    plummer = new_plummer_model(N_BODIES, random=RNG)
    
    # set the star id's to str types explicitly for concatenation later
    plummer.id = [f'{i}' for i in range(N_BODIES)] 
    
    snapshot_filename = f'{output_directory}/snapshots.log'
    history_filename = f'{output_directory}/history.txt'
    
    check_clean_directory(output_directory)

    outfile = open(history_filename, 'a+')
    doc = f" N={n_stars} - {min_hardness_kt}kT "
    outfile.write(f"{doc:=^80}\n")

    initial_ke = plummer.kinetic_energy()

    gravity.particles.add_particles(plummer)
    channel_out = gravity.particles.new_channel_to(plummer)
    t_start = 0
    for step in custom_tqdm(np.arange(t_start, t_end, step_size),
                            total=t_end/step_size):
        
        # update the cluster
        gravity.evolve_model(step * 1| nbody_system.time)
        channel_out.copy()
        
        # save snapshot
        write_set_to_file(plummer, f"{snapshot_directory}/snapshot_time{step:.2f}", "csv", overwrite=True)
        
        # decompose and filter out single stars
        ids = find_composite_multiples(plummer.copy(), min_hardness_kt=MIN_HARDNESS_KT, initial_ke = initial_ke)        
        justmultiples = [id for id in ids if len(id) > 5]
        
        # output the decomposition representation
        print(f'time: {step:.2f} - decomposition: {justmultiples}', file=outfile)
    
    print("Done!")
    gravity.stop()
    outfile.close()


if __name__ in '__main__':
    build_history(t_end=200, step_size=0.01, output_directory=OUTPUT_DIRECTORY)

