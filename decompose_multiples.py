from amuse.datamodel import Particles, Particle
from amuse.lab import nbody_system

import numpy as np


def find_composite_multiples(plummer, initial_ke, min_hardness_kt=1, for_plot: bool=False):
    """
    Succesively combines binaries in hierarchical configurations
    Returns: list[str] id's of the stars, encoding binary hierarchy
    """
    binary_replaced = True
    while binary_replaced:
        binary_replaced = replace_a_binary(plummer, initial_ke=initial_ke, hardness_kt=min_hardness_kt, for_plot=for_plot)

    return plummer.id


def replace_a_binary(plummer: Particles, initial_ke, hardness_kt: float = 1, for_plot: bool=False) -> bool | tuple:
    """
    Replaces the hardest binary in the particleset with a COM particle
    Only considers binaries with a hardness > hardness_kt
    RETURNS: boolean to indicate if a binary was replaced
      note: the particleset is edited in place, so use a copy!
    """ 
    def hardest_binary(bins):
        """Returns the hardest binary in the set of binaries"""
        max_id = np.argmax([binary[0].hardness for binary in bins])
        return bins[max_id]
    
    # convert the target hardness (in kT at initial time) 
    # to units of 3/2kT at current time:
    hardness_prefactor =  initial_ke / plummer.kinetic_energy()
    min_hardness = hardness_prefactor * 2/3 * hardness_kt

    binaries = plummer.get_binaries(G=nbody_system.G, hardness=min_hardness)

    # if there's no binaries, we're done
    if len(binaries) == 0:
        return None if for_plot else False
    
    # replace the hardest binary by a COM particle and remove the components
    target = hardest_binary(binaries)
    comparticle = collapse_binary(target, hardness_pf = hardness_prefactor, for_plot=for_plot)
    plummer.remove_particles(target)
    plummer.add_particle(comparticle)

    # return this information for plotting purposes 
    # (not used to make the history)
    if for_plot:
        return target.id, target[0].hardness * 3/2 / hardness_prefactor
    else:
        return True

def collapse_binary(binary, hardness_pf: float, for_plot:bool) -> Particle:
    """
    Returns a center-of-mass particle for the given binary
    """
    compos = binary.center_of_mass()
    comvel = binary.center_of_mass_velocity()
    commass = binary.mass.sum()
    comparticle = Particles(1)
    comparticle.mass = commass
    comparticle.velocity = comvel
    comparticle.position = compos

    if for_plot:
        comparticle.id = 10000
        return comparticle
    
    # convert hardness to initial kT units:
    hardness = binary[0].hardness * 3/2 / hardness_pf 
    
    # try to put the single stars to the left for easier reading
    if len(str(binary[1].id)) > len(str(binary[0].id)):
        comparticle.id = f"{hardness:.3f}[ {binary[0].id}, {binary[1].id}]"

    else:
        comparticle.id = f"{hardness:.3f}[ {binary[1].id}, {binary[0].id}]"

    return comparticle


if __name__ in '__main__':
    ...
    # build_history(t_end=200, step_size=0.01, output_directory='output/n64')
    # test_composite_multiples()
    # with open('decomp_log_n16.txt', 'a+') as outputfile:
        # outputfile.write(f'\n' + "="*20 + f" new run - {MIN_HARDNESS_KT}kT - seed {RNGSEED} " + "="*20 + "\n")
        # build_history(t_end=50, 
                    #   outputfile=outputfile
                    #   )

    # build_history(t_end=1)


    



