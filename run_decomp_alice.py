from amuse.ic.plummer import new_plummer_model
from amuse.lab import nbody_system
from amuse.lab import write_set_to_file
from amuse.community.ph4.interface import Ph4

import numpy as np
import os

from decompose_multiples import find_composite_multiples
from helpers import check_clean_directory, running_on_alice

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
    plummer.id = [f'{i}' for i in range(n_stars)]

    snapshot_filename = f'{output_directory}/snapshots.log' # just one file because fuck it
    history_filename = f'{output_directory}/history.txt'

    outfile = open(history_filename, 'a+')
    # outfile.write(f"="*20 + f" N={n_stars} - {min_hardness_kt}kT" + "="*20 + "\n")
    doc = f" N={n_stars} - {min_hardness_kt}kT "
    outfile.write(f"{doc:=^80}\n")

    initial_ke = plummer.kinetic_energy()
    gravity.particles.add_particles(plummer)
    channel_out = gravity.particles.new_channel_to(plummer)

    t_start = 0
    for step in np.arange(t_start, t_end, 1/snapshot_frequency):
        # update the cluster
        model_time = step * 1 | nbody_system.time
        gravity.evolve_model(model_time)
        channel_out.copy()
        
        # save snapshot
        write_set_to_file(plummer.savepoint(model_time), snapshot_filename, "amuse", append_to_file=True)
        
        # decompose and filter out single stars
        ids = find_composite_multiples(plummer.copy(), min_hardness_kt=min_hardness_kt, initial_ke=initial_ke)        
        justmultiples = [id for id in ids if len(id) > 5]

        max_inner_hardness = get_max_hardness_in_decomp_list(justmultiples)
        outfile.write(f'\ntime: {step:.2f} - decomposition: {justmultiples}')
        
        if max_inner_hardness < 10:
            patience = 0
            continue
            
        if patience > stop_delay:
            gravity.stop()
            outfile.close()
            return step
            
        patience += 1

    gravity.stop()
    outfile.close() 


def get_innermost_hardness(decomp: str) -> float:
    decomp_cleaned = decomp\
        .replace('[', '')\
        .replace(']', '')\
        .replace(',', '')
    decomp_array = np.array(decomp_cleaned.split())

    hardnesses = decomp_array[np.where(['.' in d for d in decomp_array])].astype(float)
    return hardnesses[-1]

def get_max_hardness_in_decomp_list(decomp_list: list) -> float:
    if len(decomp_list) == 0:
        return 0.
    binding_energies = [get_innermost_hardness(decomp=decomp) for decomp in decomp_list]
    return np.max(binding_energies)


def make_buncha_data(n_runs: int, n_stars: int, snapshot_frequency: int, min_hardness_kt: float) -> None:
    main_output_directory = f'output/n{n_stars}'

    if running_on_alice():
        main_output_directory = f'/home/s2015242/data1/{main_output_directory}'
    
    if not os.path.exists(main_output_directory):
        raise FileNotFoundError("Output directory not found")
    # check_clean_directory(main_output_directory)

    # TODO: prototype plots -> make sure we can parse the stuff efficiently
    # TODO: prototype binary hardness history -> traverse history backwards and keep track of hardness of... hardest binary?
    # TODO: snapshot cadence naar 2^-3 (of maybe of 2^-4)?? idk. Allebei prima, data is goedkoop I guess. Variabele!
    # EERST op alice aan zetten, dan ondertussen aan het werk. Anders zonde en moet je vet lang wachten op resultaten :)
    # TODO: when stop run? -> parse history somehow I guess
    # TODO: history parser that does non-hierarchical systems

    for _ in range(n_runs):
        i = len(os.listdir(main_output_directory))
        output_directory = f'{main_output_directory}/run_{i}'

        if os.path.exists(output_directory):
            time.sleep(np.random.randint(low=0, high=100)/10)
            i = len(os.listdir(main_output_directory))
            output_directory = f'{main_output_directory}/run_{i}'

            if os.path.exists(output_directory):
                raise FileExistsError(f"Output directories run_{i} and run_{i+1} already exists! Seems like someone did a fucky...")
        
        os.mkdir(output_directory)

        do_run(n_stars=n_stars, 
               snapshot_frequency=snapshot_frequency, 
               output_directory=output_directory, 
               min_hardness_kt=min_hardness_kt, 
               rng_seed=i)


# def test_diff_step_sizes():
#     snapshot_freqs = [2**3, 10, 2**6, 100, 2**7]
#     n_bodies = 16
#     output_dir = 'output/test'
#     min_hardness_kt = 3
#     rng_seed = 123
#     t_ccs = np.zeros(5)
#     for seed in range(10):
#         for i, sf in enumerate(snapshot_freqs):
#             t_cc = do_run(n_stars=n_bodies, snapshot_frequency=sf, output_directory=output_dir, min_hardness_kt=min_hardness_kt, rng_seed=seed)
#             print(f'frequency {sf} - t_cc {t_cc}')
#             t_ccs[i] += t_cc

#     print(f'averages: {t_ccs/10}')

# def test_read_snapshot():
#     from amuse.lab import read_set_from_file
#     snapshot_file = 'output/test/snapshots.log'
#     snapshots = read_set_from_file(snapshot_file, 'amuse')
#     snapshots_list = [s for s in snapshots.history]
#     snapshot = snapshots_list[-1]
#     ids = find_composite_multiples(snapshot.copy(), initial_ke=snapshots_list[0].kinetic_energy(), min_hardness_kt=1)
#     print(ids)


if __name__ in '__main__':
    import time
    sleep_time = np.random.randint(low=0, high=100)
    time.sleep(sleep_time/10)
    # print(sleep_time/10)
    make_buncha_data(n_runs=10, n_stars=16, snapshot_frequency=128, min_hardness_kt=1)
