from helpers import read_snapshot_file, read_history_csv, get_run_directory
from helpers import snapshot_at_time
from decompose_multiples import find_composite_multiples
from amuse.lab import nbody_system, Particles
import matplotlib.pyplot as plt

def get_before_after_snapshots() -> tuple:
    n_stars = 16
    run_id = 1
    all_snapshots = read_snapshot_file(run_directory=get_run_directory(n_stars=n_stars, run_id=run_id))
    times = [2.875, 2.883]
    times = [7.992, 8.000]
    
    times = [5.734, 5.742]
    hardest_binary_ids = ['15', '06']
    other_star_ids = ['05', '12', '10', '00', '13', '04']

    times = [0.359, 0.367]
    hardest_binary_ids = ['05', '14']
    other_star_ids = ['05', '14', '02', '00', '15', '10', '09', '04']

    snapshot_before = snapshot_at_time(snapshots=all_snapshots, time=times[0], snapshot_frequency=128)
    snapshot_after = snapshot_at_time(snapshots=all_snapshots, time=times[1], snapshot_frequency=128)
    print(snapshot_before.get_timestamp())
    print(snapshot_after.get_timestamp())
    initial_ke = all_snapshots[0].kinetic_energy()
    return snapshot_before, snapshot_after, initial_ke, hardest_binary_ids, other_star_ids


def get_decomp_for_snapshot(snapshot, initial_ke):
    ids = find_composite_multiples(plummer=snapshot.copy(), initial_ke=initial_ke, min_hardness_kt=1)
    decomp = [id for id in ids if len(id) > 5]
    return decomp



def print_distances_to_binary(snapshot, hardest_binary_ids, other_star_ids):
    snapshot = snapshot.copy()
    primary_id, secondary_id = hardest_binary_ids
    primary = snapshot[snapshot.id == primary_id]
    secondary = snapshot[snapshot.id == secondary_id]
    hardest_binary = Particles()
    hardest_binary.add_particle(primary)
    hardest_binary.add_particle(secondary)

    compos = hardest_binary.center_of_mass()
    snapshot.position -= compos
    distances = []
    for star_id in other_star_ids:
        star = snapshot[snapshot.id == star_id]
        distance = star.position.length().value_in(nbody_system.length)
        distances.append(f'{distance:.2f}')

    print(distances)
    print(other_star_ids)

    
def compare_distances():
    snapshot_before, snapshot_after, initial_ke, hardest_binary, other_star_ids = get_before_after_snapshots()
    print("="*100)
    print(get_decomp_for_snapshot(snapshot=snapshot_before, initial_ke=initial_ke))
    print_distances_to_binary(snapshot_before, hardest_binary_ids=hardest_binary, other_star_ids=other_star_ids)
    print("-"*100)
    print(get_decomp_for_snapshot(snapshot=snapshot_after, initial_ke=initial_ke))
    print_distances_to_binary(snapshot=snapshot_after, hardest_binary_ids=hardest_binary, other_star_ids=other_star_ids)
    print("="*100)

    


if __name__ == '__main__':
    compare_distances()