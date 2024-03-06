from amuse.ic.plummer import new_plummer_model
from amuse.lab import nbody_system
from amuse.lab import write_set_to_file
from amuse.community.ph4.interface import Ph4

import numpy as np
import os

from decompose_multiples import find_composite_multiples
from helpers import running_on_alice, get_final_snapshot_filename


def is_hb_hard_enough(plummer, target_hardness_kt0, initial_ke) -> bool:
    # just see if there's binaries with >10kT_0 hardness:
    hardness_prefactor = initial_ke / plummer.kinetic_energy()
    min_hardness_kt_adjusted = hardness_prefactor * 2/3 * target_hardness_kt0
    binaries = plummer.get_binaries(G=nbody_system.G, hardness=min_hardness_kt_adjusted)
    return len(binaries) > 0

    
def do_run(n_stars: int, 
           snapshot_frequency: float, 
           output_directory: str, 
           min_hardness_kt: float,
           rng_seed=None, 
           stop_delay: int = 10,
           t_end: float = 200,
           ) -> None:

    if rng_seed:
        rng = np.random.default_rng(rng_seed)
    else:
        rng = np.random.default_rng()

    gravity = Ph4()
    plummer = new_plummer_model(n_stars, random=rng)
    plummer.id = [f'{i:>02}' for i in range(n_stars)]

    snapshot_filename = f'{output_directory}/snapshots.log' # just one file because fuck it
    final_snapshot_filename = get_final_snapshot_filename(n_stars=n_stars, run_id=rng_seed) # rng_seed is the run_id. Knew that would come in handy!
    
    initial_ke = plummer.kinetic_energy()
    gravity.particles.add_particles(plummer)
    channel_out = gravity.particles.new_channel_to(plummer)

    t_start = 0
    patience = 0
    for step in np.arange(t_start, t_end, 1/snapshot_frequency):
        # update the cluster
        model_time = step * 1 | nbody_system.time
        gravity.evolve_model(model_time)
        channel_out.copy()
        
        # save snapshot
        write_set_to_file(plummer.savepoint(model_time), snapshot_filename, "amuse", append_to_file=True)
        
        # see if we're done
        hb_hard_enough = is_hb_hard_enough(plummer=plummer, 
                                           target_hardness_kt0=10, 
                                           initial_ke=initial_ke
                                           )
    
        if hb_hard_enough:
            patience += 1
        else:
            patience = 0
            
        if patience > stop_delay:
            break

    write_set_to_file(plummer.savepoint(model_time), final_snapshot_filename, "amuse")

    gravity.stop()


def make_buncha_data(n_runs: int, n_stars: int, snapshot_frequency: int, min_hardness_kt: float) -> None:
    main_output_directory = f'output/n{n_stars}'

    if running_on_alice():
        main_output_directory = f'/home/s2015242/data1/{main_output_directory}'
    else:
        main_output_directory = f'/home/jasper/studie/MRP/{main_output_directory}'

    if not os.path.exists(main_output_directory):
        print(os.getcwd())
        raise FileNotFoundError(f"Output directory {main_output_directory} not found")

    # TODO: prototype binary hardness history -> traverse history backwards and keep track of hardness of... hardest binary?
    
    # TODO: refactor this so either everything happens in pp, or the updated decomp stuff happens here. No use doing the decomps if I'm just going to redo them

    get_i = lambda: len(os.listdir(main_output_directory))

    def claim_directory():
        i = get_i()
        output_directory = f'{main_output_directory}/run_{i:>04}'
        try:
            os.mkdir(output_directory)
        except FileExistsError:
            # if the file already exists, then try another one. 
            # things get kind of messy when you're running 12 runs in parallel
            return claim_directory()
        
        return output_directory, i

    for _ in range(n_runs):
        output_directory, i = claim_directory()

        do_run(n_stars=n_stars, 
               snapshot_frequency=snapshot_frequency, 
               output_directory=output_directory, 
               min_hardness_kt=min_hardness_kt, 
               rng_seed=i,
               )



if __name__ in '__main__':
    ...
    # import time
    # test_get_max_h_in_decomp_list()
    # sleep_time = np.random.randint(low=0, high=100)
    # time.sleep(sleep_time/10)
    # print(sleep_time/10)
    # make_buncha_data(n_runs=1, n_stars=16, snapshot_frequency=128, min_hardness_kt=1)
