import numpy as np

from amuse.lab import Particles
from helpers import get_final_snapshot, read_history_csv, get_run_ids_for_n_stars, custom_tqdm, get_output_path
from heritage_from_history import get_hardest_binary_from_decomp
from amuse.lab import nbody_system
from amuse.ext.orbital_elements import orbital_elements_from_binary



def get_final_binary_ids_in_run(n_stars: int, run_id: int):
    *_, hbhs,  decomps = read_history_csv(n_stars=n_stars, run_id=run_id)
    if hbhs[-1] < 10:
        return None
    
    final_decomp = decomps[-1]

    _, final_binary_ids = get_hardest_binary_from_decomp(decomp=final_decomp)

    return final_binary_ids


def get_final_binary_in_run(n_stars: int, run_id: int):
    final_binary_ids = get_final_binary_ids_in_run(n_stars=n_stars, run_id=run_id)
    if final_binary_ids is None:
        return None
    
    final_snapshot = get_final_snapshot(n_stars=n_stars, run_id=run_id)

    primary_id, secondary_id = final_binary_ids
    final_binary = Particles()
    final_binary.add_particle(final_snapshot[final_snapshot.id == primary_id])
    final_binary.add_particle(final_snapshot[final_snapshot.id == secondary_id])

    return final_binary, final_snapshot

def get_cluster_core_pos(snapshot: Particles):
    core_pos, _, _ = snapshot.densitycentre_coreradius_coredens(number_of_neighbours=5)
    return core_pos

def get_final_bindist_ecc_for_run(n_stars: int, run_id: int) -> float:
    final_binary, final_snapshot = get_final_binary_in_run(n_stars=n_stars, run_id=run_id)
    if final_binary is None:
        return -1, -1
    
    core_pos = get_cluster_core_pos(final_snapshot)

    distance_to_cc = (final_binary.position - core_pos).length().value_in(nbody_system.length)
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
    for i, run_id in custom_tqdm(enumerate(run_ids), total=len(run_ids)):
        print(run_id)
        distance, ecc = get_final_bindist_ecc_for_run(n_stars=n_stars, run_id=run_id) 
        binary_props[i] = np.array([distance, ecc])

    np.save(file=npz_filename, arr=binary_props)
    return binary_props.T


def test_densitycentre_finding():
    from amuse.datamodel.particle_attributes import particle_potential
    def myown_bound_subset(plummer):
        distances_from_center = plummer.position.lengths().value_in(nbody_system.length)

        farenough_stars = plummer[distances_from_center > np.percentile(distances_from_center, q=10)]
        for star in farenough_stars:
            total_energy = (particle_potential(plummer, star, G=nbody_system.G) + 0.5 * (star.velocity**2).sum()).value_in(nbody_system.length**2 * nbody_system.time**(-2))
            # if ke > 0:
            print(f'{star.id} - {total_energy}')
        
        print(len(plummer.bound_subset(G=nbody_system.G)))
        print(plummer.bound_subset(G=nbody_system.G).id)
        kes = np.array([(particle_potential(plummer, star, G=nbody_system.G) + 0.5 * (star.velocity**2).sum()).value_in(nbody_system.length**2 * nbody_system.time**(-2)) for star in farenough_stars])
        print(len(kes[kes>=0]))
        # print(farenough_stars.id)

    n_stars = 64
    run_id = 2

    final_snapshot = get_final_snapshot(run_id=run_id, n_stars=n_stars)
    myown_bound_subset(final_snapshot)
    # print(final_snapshot.densitycentre_coreradius_coredens(number_of_neighbours=5))

    # bound_subset = final_snapshot.bound_subset(G=nbody_system.G)
    # print(bound_subset.id)
    # print(bound_subset.densitycentre_coreradius_coredens(number_of_neighbours=5))

    # print(get_final_bindist_ecc_for_run(n_stars=n_stars, run_id=run_id))


def test_dc():
    run_id = 10
    n_stars = 16

    print(get_final_binary_in_run(n_stars=n_stars, run_id=run_id)[0].position.length().value_in(nbody_system.length))
    print(get_final_bindist_ecc_for_run(run_id=run_id, n_stars=n_stars))

if __name__ == '__main__':
    test_dc()