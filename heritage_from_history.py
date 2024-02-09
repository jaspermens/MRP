import numpy as np
from helpers import get_run_directory


def read_history(n_stars: int = 16, run_id: int = 0):
    run_directory = get_run_directory(n_stars=n_stars, run_id=run_id)
    history_fn = f'{run_directory}/history.txt'

    with open(history_fn, 'r') as file:
        # split into lines and skip the header and footer:
        history = file.read().split('\n')[2:-2] 

    return history


def time_decomp_from_line(line: str):
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
    binary_slices = decomp.split('],')

    hardnesses = np.zeros_like(binary_slices, dtype=object)
    components = np.zeros_like(binary_slices, dtype=object)

    for i, binary_slice in enumerate(binary_slices):
        hardnesses[i], components[i] = extract_binary_info_from_slice(binary_slice)

    return hardnesses, components


def get_hardest_binary_from_decomp(decomp: str):
    hardnesses, binaries = extract_binaries_from_decomp(decomp=decomp)

    hardest_binary_index = np.argmax(hardnesses)
    hardest_binary = binaries[hardest_binary_index]
    hardest_hardness = hardnesses[hardest_binary_index]

    return hardest_hardness, hardest_binary


def get_theseus_binary_from_decomp(target_binary: list, decomp: str):
    """
    ew ew ew ew please rewrite
    there HAS to be a numpy thing for this...
    """
    hardnesses, binaries = extract_binaries_from_decomp(decomp=decomp)

    member_containing_hardnesses = []
    member_containing_binaries = []
    for star in target_binary:
        for bi, binary in enumerate(binaries):
            if star not in binary:
                continue
            member_containing_hardnesses.append(hardnesses[bi])
            member_containing_binaries.append(binaries[bi])

    if len(member_containing_binaries) == 0:
        return 0, target_binary
    
    hardest_id = np.argmax(member_containing_hardnesses)
    hardest_member_hardness = member_containing_hardnesses[hardest_id]
    hardest_member_binary = member_containing_binaries[hardest_id]
    return hardest_member_hardness, hardest_member_binary


def compile_heritage(n_stars: int = 16, run_id: int = 0):
    """Returns time and data arrays for the hardness history of the final hard binary"""
    history_backwards = read_history(n_stars=n_stars, run_id=run_id)[::-1]

    final_time, decomp = time_decomp_from_line(history_backwards[0])
    final_hardness, final_components = get_hardest_binary_from_decomp(decomp)
    
    hardnesses = [final_hardness]
    times = [final_time]
    binaries = [final_components]

    binary = final_components
    patience = 0
    for line in history_backwards[1:]:
        if patience == 10:
            break
        time, decomp = time_decomp_from_line(line)

        times.append(time)
        hardness, binary = get_theseus_binary_from_decomp(binary, decomp)
        
        if hardness == 0:
            patience += 1
        
        hardnesses.append(hardness)
        binaries.append(binary)

    return times, np.array(hardnesses).astype('float'), binaries



if __name__ in '__main__':
    print(compile_heritage(run_id=0))
    # test_decomp_str_from_line(time=6.10)