from amuse.ic.plummer import new_plummer_model
from amuse.lab import nbody_system
from amuse.lab import write_set_to_file
from amuse.community.ph4.interface import Ph4

import numpy as np
import os

from decompose_multiples import find_composite_multiples
from helpers import running_on_alice


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
    history_filename = f'{output_directory}/history.txt'

    outfile = open(history_filename, 'a+')
    # outfile.write(f"="*20 + f" N={n_stars} - {min_hardness_kt}kT" + "="*20 + "\n")
    doc = f" N={n_stars} - {min_hardness_kt}kT "
    outfile.write(f"{doc:=^80}\n")

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
        
        # decompose and filter out single stars
        ids = find_composite_multiples(plummer.copy(), min_hardness_kt=min_hardness_kt, initial_ke=initial_ke)        
        justmultiples = [id for id in ids if len(id) > 5]

        max_inner_hardness = get_max_hardness_in_decomp_list(justmultiples)
        outfile.write(f'\ntime: {step:.3f} - decomposition: {justmultiples}')
        
        if max_inner_hardness < 10:
            patience = 0
            continue
            
        if patience > stop_delay:
            break
            
        patience += 1

    outfile.write(f'\nEnd time:\n{model_time.value_in(nbody_system.time):.2f}')
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

    binding_energies = [get_innermost_hardness(decomp=binary_slice) for full_decomp in decomp_list for binary_slice in full_decomp.split("],")]
    return np.max(binding_energies)


def test_get_max_h_in_decomp_list():
    decomp_list = ['22.8[ 31, 50]', '4.0[ 29, 3.0[ 48, 5.9[ 04, 4.4[ 24, 15]]]]']
    decomp_list = ['1.5[ 19, 1.2[ 22.8[ 31, 50], 9.2[ 48, 3.3[ 04, 4.0[ 29, 3.3[ 24, 15]]]]]]']
    print(get_max_hardness_in_decomp_list(decomp_list=decomp_list))

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
    
    # TODO: test multicore speed => seems to work just fine! Not great, but fine! 4 processes did it 3x as fast. 

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
    # import time
    test_get_max_h_in_decomp_list()
    # sleep_time = np.random.randint(low=0, high=100)
    # time.sleep(sleep_time/10)
    # print(sleep_time/10)
    # make_buncha_data(n_runs=1, n_stars=16, snapshot_frequency=128, min_hardness_kt=1)
