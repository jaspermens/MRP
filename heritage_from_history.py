import numpy as np
from helpers import get_run_directory, read_history_csv


def read_history_deprecated(n_stars: int = 16, run_id: int = 0):
    run_directory = get_run_directory(n_stars=n_stars, run_id=run_id)
    history_fn = f'{run_directory}/history.txt'

    with open(history_fn, 'r') as file:
        # split into lines and skip the header and footer:
        history = file.read().split('\n')[2:-2] 

    return history


def time_decomp_from_line_deprecated(line: str):
    time, decomp = line.split(' - ')
    decomp_raw = decomp.split(':')[1]
    time = float(time.split(':')[1])

    return time, decomp_raw


def extract_binary_info_from_slice(decomp_slice: str) -> tuple:
    decomp_cleaned = decomp_slice\
        .replace('[', '')\
        .replace(']', '')\
        .replace(',', '')\
        .replace("'", '')
    
    if len(decomp_cleaned) <= 1:
        return 0, [None, None]
    
    decomp_array = np.array(decomp_cleaned.split())
    hardness, star1, star2 = decomp_array[-3:]
    # hardnesses = decomp_array[np.where(['.' in d for d in decomp_array])].astype(float)
    return hardness, [star1, star2]


def extract_binaries_from_decomp(decomp: str):
    binary_slices = [b for binary in decomp.split("', '") for b in binary.split('],')]

    hardnesses = np.zeros_like(binary_slices, dtype=object)
    components = np.zeros_like(binary_slices, dtype=object)

    for i, binary_slice in enumerate(binary_slices):
        hardnesses[i], components[i] = extract_binary_info_from_slice(binary_slice)

    return hardnesses.astype(float), components


def get_hardest_binary_from_decomp(decomp: str):
    hardnesses, binaries = extract_binaries_from_decomp(decomp=decomp)

    hardest_binary_index = np.argmax(hardnesses)
    hardest_binary = binaries[hardest_binary_index]
    hardest_hardness = hardnesses[hardest_binary_index]

    return hardest_hardness, hardest_binary


def test_hardest_binary_from_decomp():
    decomp1 = "['1.9[ 03, 1.8[ 01, 3.3[ 11, 3.5[ 14, 6.9[ 00, 4.6[ 08, 10.5[ 12, 8.2[ 13, 12.3[ 10, 2.2[ 12.7[ 05, 15], 1.0[ 04, 3.8[ 06, 1.4[ 02, 4.7[ 07, 09]]]]]]]]]]]]]]']"
    decomp2 = "['1.0[ 11, 1.2[ 3.3[ 00, 03], 1.5[ 09, 1.0[ 13, 5.1[ 05, 8.2[ 06, 3.6[ 07, 1.3[ 01, 12]]]]]]]]']"
    print(get_hardest_binary_from_decomp(decomp=decomp2))


def get_theseus_binary_from_decomp(target_binary: list, decomp: str):
    """
    ew ew ew ew please rewrite
    there HAS to be a numpy thing for this...
    """
    hardnesses, binaries = extract_binaries_from_decomp(decomp=decomp)

    member_containing_hardnesses = []
    member_containing_binaries = []
    target_primary, target_secondary = target_binary
    
    for hardness, binary in zip(hardnesses, binaries):
        if (target_primary not in binary) and (target_secondary not in binary):
            continue
        member_containing_hardnesses.append(hardness)
        member_containing_binaries.append(binary)

    if len(member_containing_binaries) == 0:
        return 0, target_binary
    
    hardest_id = np.argmax(member_containing_hardnesses)
    hardest_member_hardness = member_containing_hardnesses[hardest_id]
    hardest_member_binary = member_containing_binaries[hardest_id]

    return hardest_member_hardness, hardest_member_binary


def compile_heritage(n_stars: int = 16, run_id: int = 0):
    """Returns time and data arrays for the hardness history of the final hard binary"""
    times, *_, decomps = read_history_csv(n_stars=n_stars, run_id=run_id)
    times_reversed = times[::-1]
    decomps_reversed = decomps[::-1]
    final_time, final_decomp = times_reversed[0], decomps_reversed[0]
    assert final_time > 0

    final_hardness, final_components = get_hardest_binary_from_decomp(final_decomp)
    # maybe this should be: 
    # if None in final_components:
    # this current version also cuts out the DNF's
    if final_hardness < 10:
        return [], np.array([], dtype=float), []
    history_hardnesses = [final_hardness]
    history_times = [final_time]
    history_binaries = [final_components]

    binary = final_components
    patience = 0
    for time, decomp in zip(times_reversed, decomps_reversed):
        if patience == 10:
            break
        history_times.append(time)
        hardness, binary = get_theseus_binary_from_decomp(binary, decomp)
        
        if hardness == 0:
            patience += 1
        
        history_hardnesses.append(hardness)
        history_binaries.append(binary)

    return history_times, np.array(history_hardnesses).astype('float'), history_binaries



if __name__ in '__main__':
    print(compile_heritage(n_stars=16, run_id=1))

    # test_hardest_binary_from_decomp()
    # test_decomp_str_from_line(time=6.10)