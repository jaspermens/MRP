import numpy as np

from amuse.lab import Particles
from helpers import get_final_snapshot, read_history_csv, get_run_ids_for_n_stars, custom_tqdm, get_output_path
from heritage_from_history import get_hardest_binary_from_decomp
from amuse.lab import nbody_system
from amuse.ext.orbital_elements import orbital_elements_from_binary


def get_final_binary_in_run(n_stars: int, run_id: int):
    final_snapshot = get_final_snapshot(n_stars=n_stars, run_id=run_id)
    final_binary_ids = get_final_binary_ids_in_run(n_stars=n_stars, run_id=run_id)

    primary_id, secondary_id = final_binary_ids
    final_binary = Particles()
    final_binary.add_particle(final_snapshot[final_snapshot.id == primary_id])
    final_binary.add_particle(final_snapshot[final_snapshot.id == secondary_id])

    return final_binary


def get_final_binary_ids_in_run(n_stars: int, run_id: int):
    *_, decomps = read_history_csv(n_stars=n_stars, run_id=run_id)
    final_decomp = decomps[-1]

    _, final_binary_ids = get_hardest_binary_from_decomp(decomp=final_decomp)

    return final_binary_ids


def get_final_bindist_ecc_for_run(n_stars: int, run_id: int) -> float:
    final_binary = get_final_binary_in_run(n_stars=n_stars, run_id=run_id)
    distance_to_cc = final_binary.center_of_mass().length().value_in(nbody_system.length)
    eccentricity = orbital_elements_from_binary(final_binary)[3]

    return distance_to_cc, eccentricity


def get_final_binary_properties_for_n_stars(n_stars: int):
    output_path = get_output_path()
    npz_filename = f'{output_path}/n{n_stars}_finalbinary_properties.npy'
    try:
        binary_props = np.load(file=npz_filename)
        print(f'File found! re-using the {len(binary_props)} binary properties from {npz_filename}...')
        return binary_props.T
    
    except FileNotFoundError:
        print(f'File {npz_filename} not found- redoing')
        pass

    run_ids = get_run_ids_for_n_stars(n_stars=n_stars)
    
    binary_props = np.zeros((len(run_ids), 2), dtype=float)
    for i,run_id in custom_tqdm(enumerate(run_ids), total=len(run_ids)):
        binary_props[i] = get_final_bindist_ecc_for_run(n_stars=n_stars, run_id=run_id) 

    np.save(file=npz_filename, arr=binary_props)
    return binary_props.T


if __name__ == '__main__':
    ...