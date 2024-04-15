import numpy as np
from amuse.lab import Particles, nbody_system, constants
from helpers import custom_tqdm, read_snapshot_file, get_run_directory
from decompose_multiples import find_composite_multiples
import os
from amuse.io.base import IoException

# ok so now we're going to traverse backward. 
# we get the final binary from the final snapshot
# we get the initial temperature from the first snapshot (could re-init using the seed too, but I doubt that's faster than just accessing the snapshot)
# for each snapshot, we get time, temperature, fhb_hardness, Edot per star, r per star, decomp for like eight stars/iterations (so we can see scrambles but never decomp all 128 stars)
# write to file and we're done!

def get_binaries_amuse(particles,hardness=10,G = constants.G):
    """
    returns the binaries in a particleset. binaries are selected according to a hardness criterion [hardness=10]
    This function returns the binaries as a list of i,j particles. Triple detection is not done.
    
    >>> from amuse import datamodel
    >>> m = [1,1,1] | units.MSun
    >>> x = [-1,1,0] | units.AU
    >>> y = [0,0,1000] | units.AU
    >>> z = [0,0,0] | units.AU
    >>> vx = [0,0,0] | units.kms
    >>> vy = [1.,-1.,0] | units.kms
    >>> vz = [0,0,0] | units.kms
    >>> particles = datamodel.create_particle_set( mass=m,x=x,y=y,z=z,vx=vx,vy=vy,vz=vz )
    >>> binaries = particles.get_binaries()
    >>> print len(binaries)
    1
    
    """
    n=len(particles)
    total_Ek=(0.5*particles.mass*(particles.vx**2+particles.vy**2+particles.vz**2)).sum()
    average_Ek=total_Ek/particles.mass.sum()
    max_mass=particles.mass.amax()
    limitE=hardness*average_Ek

    a=np.argsort(particles.x.number)

    binaries=[]

    for i in range(n-1):
        j=i+1
        while j<n and (particles.x[a[j]]-particles.x[a[i]])<2*G*max_mass/limitE:
            r2=(particles.x[a[j]]-particles.x[a[i]])**2+ \
               (particles.y[a[j]]-particles.y[a[i]])**2+ \
               (particles.z[a[j]]-particles.z[a[i]])**2 
            v2=(particles.vx[a[j]]-particles.vx[a[i]])**2+ \
               (particles.vy[a[j]]-particles.vy[a[i]])**2+ \
               (particles.vz[a[j]]-particles.vz[a[i]])**2 
            r=r2**0.5
            eb=G*(particles.mass[a[i]]+particles.mass[a[j]])/r-0.5*v2
            if eb > limitE:
                binary=particles[[a[i],a[j]]].copy()
                binary.hardness=eb/average_Ek
                binaries.append(binary)
            j+=1  

    return binaries

def get_binary_with_star_amuseversion(star_id: int, particles: Particles, min_hardness: float):
    binaries = particles.get_binaries(G=nbody_system.G, hardness=min_hardness)
    candidates = []
    for binary in binaries:
        for comp in binary:
            if int(comp.id) == star_id:
                candidates.append(binary)

    hardnesses = [c.hardness[0] for c in candidates]
    binary = candidates[np.argmax(hardnesses)]
    
    return binary.id, binary.hardness[0], binary        

## Theseus binary + hardness
def get_binary_with_star(star_id: int, particles: Particles, min_hardness: float):
    total_Ek=(0.5*particles.mass*(particles.vx**2+particles.vy**2+particles.vz**2)).sum()
    average_Ek=total_Ek/particles.mass.sum()
    
    dists = (particles.position - particles[star_id].position).lengths()

    closest_stars = np.argsort(dists)
    assert star_id == closest_stars[0]
    
    star_id = closest_stars[0]
    n_neighbors = 3

    closest_n_ids = closest_stars[1:n_neighbors+1] 
    hardnesses = (nbody_system.G * 2 * particles[0].mass / dists[closest_n_ids] - 0.5*(particles[closest_n_ids].velocity - particles[star_id].velocity).lengths_squared()) / average_Ek
    hardest_companion = np.argmax(hardnesses)
    hardest_hardness = hardnesses[hardest_companion]
    hardest_id = closest_n_ids[hardest_companion]

    if hardest_hardness < min_hardness:
        return 0, 0, 0
    # because is it really still a binary if there's an unbound star inside?

    binary_ids = (star_id, hardest_id)
    binary = particles[[star_id, hardest_id]].copy()
    
    # binary = particles[[star_id, closest_id]].copy()
    # binary.hardness = hardest_hardness

    return binary_ids, hardest_hardness, binary

def theseus_binary_in_snapshot(snapshot: Particles, 
                               last_binary: Particles,    #  technically the next binary. used as an educated guess
                               patience: int,             #  patience counter for demo systems
                               min_hardness_kt0: float,   #  min hardness in units of kT0 (so not converted yet)
                               hardness_prefactor         # conversion factor to current kT
                               ):
    """ 
    gets the ancestor binary in the next/previous snapshot. 
    should be faster than n squared at least most of the time
    if we cannot find any binaries, we're probably in some kind of demo system, 
    or we're at the end/start. So, let's have a (generous) patience counter.
    """

    max_hardness = 0
    for component_id in last_binary.id:
        _, candidate_hardness, candidate_binary = get_binary_with_star(star_id=int(component_id), 
                                                                                          particles=snapshot, 
                                                                                          min_hardness=min_hardness_kt0 * hardness_prefactor)
        if candidate_hardness > max_hardness:
            max_hardness = candidate_hardness
            binary = candidate_binary

    hardness = max_hardness/hardness_prefactor  

    if max_hardness == 0:
        patience += 1
        binary = last_binary
    
    else:
        patience = 0

    return binary, hardness, patience


## Edots
def edots_for_snapshot(snapshot: Particles, binary: Particles):
    primary_pos, secondary_pos = binary.position.lengths().value_in(nbody_system.length)
    other_pos = snapshot.position.lengths().value_in(nbody_system.length)

    primary_pf = powerfunc(component_pos=primary_pos, other_pos=other_pos) 
    secondary_pf = powerfunc(component_pos=secondary_pos, other_pos=other_pos) 

    cmvel = binary.center_of_mass_velocity().length().value_in(nbody_system.length / nbody_system.time)
    primary_vel, secondary_vel = binary.velocity.lengths().value_in(nbody_system.length / nbody_system.time)
    
    return - (primary_pf * (primary_vel - cmvel) + secondary_pf * (secondary_vel - cmvel))

def powerfunc(component_pos: np.ndarray, other_pos: np.ndarray):
    return - (component_pos - other_pos) / (np.linalg.norm(component_pos - other_pos)**3)


def theseus_binary_generator(snapshots: list[Particles],
                             min_hardness_kt0: float,
                             initial_ke,
                             final_binary: Particles,
                             ):

    binary = final_binary
    for snapshot in snapshots:
        # to convert the kT0 hardness to current temperature hardness:
        hardness_prefactor = 2/3 * initial_ke / snapshot.kinetic_energy()

        binary, hardness, patience = theseus_binary_in_snapshot(snapshot=snapshot,
                                                                last_binary=binary,
                                                                patience=patience,
                                                                min_hardness_kt0=min_hardness_kt0,
                                                                hardness_prefactor=hardness_prefactor)
        hardness_kt0 = hardness / hardness_prefactor
        
        yield binary, hardness_kt0

## logistics
def claim_run(run_id: int, n_stars: int):
    run_directory = get_run_directory(n_stars=n_stars, run_id=run_id)
    claimed_flag = f'{run_directory}/claimed.txt'

    os.mknod(claimed_flag) # will raise a FileExistsError if the file already exists
    # if os.path.exists(f'{run_directory}/decompositions.txt'):
    #     raise FileExistsError  # this is just a double check- technically not necessary 


def redo_decomp(n_stars: int, start_num: int = 0):
    run_id = start_num
    while True:
        try:
            # try to claim the run
            claim_run(run_id=run_id, n_stars=n_stars)
            # if successful, redo the decomposition (and put the result in the claimed file)
            redo_decomp_for_run(run_id=run_id, n_stars=n_stars)
            run_id += 1

        except FileExistsError:
            # if the file already exists, then a run is either done or in progress
            # move on to the next one
            run_id += 1

        except FileNotFoundError:
            # then the run directory doesn't exist, so we're done (or something is borked)
            print('end of the line! No more runs to do!')
            break


## main function really
def redo_decomp_for_run(run_id: int, n_stars: int, min_hardness_kt0: float = 0.1) -> None:
    run_directory = get_run_directory(n_stars=n_stars, run_id=run_id)

    snapshots = read_snapshot_file(run_directory=run_directory)
    
    header = "time - T/T0 - binary - hardness - nearest neighbors - distances to com - energy exchange rates - patience"
    
    pickle_array = np.zeros(shape=(len(snapshots), len(header.split(' - '))), dtype=object)
    all_distances = np.zeros(shape=(len(snapshots), n_stars), dtype=float)
    all_edots = np.zeros_like(all_distances)
    all_hardnesses = np.zeros(shape=(len(snapshots)), dtype=float)
    csvlines = np.zeros(shape=(len(snapshots), len(header.split(' - '))), dtype='U2000')

    initial_ke = snapshots[0].kinetic_energy()
    
    final_hardness_prefactor = 2/3 * initial_ke / snapshots[-1].kinetic_energy()
    try:
        final_binary = snapshots[-1].get_binaries(hardness=10*final_hardness_prefactor, G=nbody_system.G)[0]
    except IndexError:
        np.save(f"{run_directory}/distances", all_distances)
        np.save(f"{run_directory}/edots", all_edots)
        np.save(f"{run_directory}/hardnesses", all_hardnesses)
        np.savetxt(f"{run_directory}/postprocess3.txt", csvlines, fmt="%s  %s  %s  %3s  %s  %s  %s  %s", header=header)
        np.save(f"{run_directory}/postprocess3", pickle_array, allow_pickle=True)
        return
    
    current_binary = final_binary
    patience = 0


    for i, snapshot in enumerate(snapshots[::-1]):
        time = snapshot.get_timestamp().value_in(nbody_system.time)
        line_time = f"{time:.3f}"
        
        current_ke = snapshot.kinetic_energy()
        temperature = current_ke/initial_ke

        line_temperature = f"{temperature:.3f}"
        hardness_prefactor = 2/3 * initial_ke / snapshot.kinetic_energy()

        # grab the binary info
        current_binary, current_hardness, patience = theseus_binary_in_snapshot(snapshot=snapshot, 
                                                                                last_binary=current_binary, 
                                                                                patience=patience, 
                                                                                min_hardness_kt0=min_hardness_kt0, 
                                                                                hardness_prefactor=hardness_prefactor
                                                                                )
        line_binary_ids = str(current_binary.id).replace("'", "")
        line_hardness = f"{current_hardness:.3f}"
        line_patience = patience

        # grab edots, dists, neighbors:
        dists = (current_binary.center_of_mass() - snapshot.position).lengths().value_in(nbody_system.length)

        closest_star_ids = [f"{name:>02}" for name in dists.argsort()[:6]]

        edots = edots_for_snapshot(snapshot=snapshot, binary=current_binary)

        current_binary_ids = current_binary.id.astype(int)
        edots[current_binary_ids] = 0

        line_edots = str(edots.astype("U3")).replace("\n", "").replace("'", "")
        line_comdists = str(dists.astype("U3")).replace("\n", "").replace("'", "")
        line_friends = str(closest_star_ids).replace("'", "")

        line = np.array([line_time, 
                         line_temperature, 
                         str(line_binary_ids), 
                         line_hardness, 
                         line_friends, 
                         line_comdists, 
                         line_edots, 
                         line_patience,
                         ], dtype="U2000")

        csvlines[-i-1, :] = line
        all_distances[-i-1, :] = dists
        all_edots[-i-1, :] = edots
        all_hardnesses[-i-1] = current_hardness

        pickle_array[-i-1, :] = np.array([time, 
                                       temperature, 
                                       current_binary.id, 
                                       current_hardness, 
                                       closest_star_ids, 
                                       dists, 
                                       edots, 
                                       patience], dtype=object)

    np.save(f"{run_directory}/distances", all_distances)
    np.save(f"{run_directory}/edots", all_edots)
    np.save(f"{run_directory}/hardnesses", all_hardnesses)
    np.savetxt(f"{run_directory}/postprocess3.txt", csvlines, fmt="%s  %s  %s  %3s  %s  %s  %s  %s", header=header)
    np.save(f"{run_directory}/postprocess3", pickle_array, allow_pickle=True)

def test_theseus_binary():
    n_stars = 16
    run_id = 45

    run_directory = get_run_directory(n_stars=n_stars, run_id=run_id)

    snapshots = read_snapshot_file(run_directory=run_directory)
    initial_ke = snapshots[0].kinetic_energy()
    final_hardness_prefactor = 2/3 * initial_ke / snapshots[-1].kinetic_energy()
    final_binary = snapshots[-1].get_binaries(hardness=10*final_hardness_prefactor, G=nbody_system.G)[0]
    b = final_binary
    p = 0
    ps = [p]
    hs = [final_binary.hardness[0] / final_hardness_prefactor]
    kes = [snapshots[-1].kinetic_energy()/initial_ke]
    swaptimes = []
    for snapshot in snapshots[::-1]:
        ke = snapshot.kinetic_energy()
        hardness_prefactor = 2/3 * initial_ke / ke
        newbinary, h, p = theseus_binary_in_snapshot(snapshot=snapshot, last_binary=b, patience=p, min_hardness_kt0=1, hardness_prefactor=hardness_prefactor)
        if str(newbinary.id) != str(b.id):
            print(f"SWAP! t = {snapshot.get_timestamp()} - from {b.id} to {newbinary.id}")
            swaptimes.append(snapshot.get_timestamp().value_in(nbody_system.time))
        if p != 0:
            print(f"PATIENT! t = {snapshot.get_timestamp()} - p = {p}")
        b = newbinary
        kes.append(ke/initial_ke)
        ps.append(p)
        hs.append(h)
    
    np.save(file='patiences', arr=np.array(ps))
    np.save(file='hardnesses', arr=np.array(hs))
    np.save(file='kes', arr=np.array(kes))
    import matplotlib.pyplot as plt

    ps = np.load(file='patiences.npy')
    hs = np.load(file='hardnesses.npy')
    x_axis = np.arange(0, len(ps)/128, 1/128)
    plt.plot(x_axis, ps[::-1])
    plt.vlines(swaptimes, ymin=0, ymax=max(ps), colors='black', ls='--')
    plt.plot(x_axis, kes[::-1])
    plt.plot(x_axis, hs[::-1])
    plt.show()


def hist_dists():
    run_id = 17
    n_stars = 16
    
    run_directory = get_run_directory(run_id=run_id, n_stars=n_stars)  
    dists = np.load(f"{run_directory}/distances.npy")
    density_proxies = np.array([dists[:, i] for i in range(16)]).T

    import matplotlib.pyplot as plt
    plt.plot(dists)
    # plt.plot(np.load(f"{run_directory}/hardnesses.npy"))
    plt.show()

if __name__ == '__main__':
    # TODO: handle DNFS
    # TODO: handle escaped stars?
    # TODO: solidify header
    # TODO: edot and r formatting -> 4 digits?
    # TODO: decomp? or like argsort of distance and edots, closest 5 ish neighbors (to com???) for legibility
    # TODO: also write distances and edots out to npz while we're at it right?
    # test_theseus_binary()
    # hist_dists()
    for i in custom_tqdm(range(50), total=50):
        redo_decomp_for_run(run_id=i, n_stars=16, min_hardness_kt0=0.1)