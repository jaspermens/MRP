from amuse.datamodel import Particles
from amuse.lab import nbody_system
from amuse.lab import read_set_from_file

import matplotlib.pyplot as plt
import numpy as np

import matplotlib.pyplot as plt
from amuse.plot import scatter
from decompose_multiples import find_composite_multiples
from helpers import get_run_directory, snapshot_at_time, read_snapshot_file

N_STARS = 16
RUN_NUMBER = 0
SNAPSHOT_DIRECTORY = get_run_directory(n_stars=N_STARS, run_id=RUN_NUMBER)

SNAPSHOTS = read_snapshot_file(run_directory=SNAPSHOT_DIRECTORY)
initial_plummer = SNAPSHOTS[0]
INITIAL_KE = initial_plummer.kinetic_energy()


def plot_multiples_for_plummer(plummer):
    """Don't use this- it's trash"""
    multiples = find_composite_multiples(plummer=plummer.copy(), min_hardness_kt=1, initial_ke=INITIAL_KE, for_plot=False)
    scatter(plummer.x, plummer.y, alpha=0.1, color='grey')
    decomp_ids = []
    nstars = len(plummer)
    for multiple in multiples:
        for id in multiple[0]:
            if id > nstars:
                continue
            decomp_ids.append(int(id))

    # decomp_ids = list(set(decomp_ids))
    print(decomp_ids)
    # binarymembers = []
    for id in decomp_ids:
        # binarymembers.append(plummer[int(id)])
        member = plummer[id]
        scatter(member.x, member.y, alpha=1, color='black')

    plt.show()


import tkinter as tk
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk 
import matplotlib
matplotlib.use('TkAgg')


def pos_all(plummer):
    return (plummer.x.value_in(nbody_system.length), 
            plummer.y.value_in(nbody_system.length),
            plummer.z.value_in(nbody_system.length))


def plot(param):
    global fig, ax1, ax2, window, canvas
    
    ax1.clear()
    ax2.clear()
    time = plot_time_slider.get()
    viewsize = get_view_size()
        
    plummer = snapshot_at_time(snapshots=SNAPSHOTS, time=time)
    # plummer = test_plummer()

    def get_focus_star_pos():
        focus_star_id = get_focus_star_id()
        if focus_star_id == -1:
            focus_x, focus_y, focus_z = 0, 0, 0
        else:
            focus_star = plummer[plummer.id == focus_star_id]
            focus_x = focus_star.position.x.value_in(nbody_system.length)
            focus_y = focus_star.position.y.value_in(nbody_system.length)
            focus_z = focus_star.position.z.value_in(nbody_system.length)

        return focus_x, focus_y, focus_z
    
    focus_x, focus_y, focus_z = get_focus_star_pos()
    
    # def dim_to_alpha(dim):
    #     absmax = np.max(np.abs(dim))
    #     neg_one_to_one_ish = np.array(dim) / absmax
    #     zero_to_one_to_zero = 1 - np.abs(neg_one_to_one_ish)
    #     return zero_to_one_to_zero

    def dim_to_alpha_centered(dim, center):
        scale = get_view_size()
        dist_from_center = np.abs(dim - center)
        dist_from_center_norm = dist_from_center / scale
        one_to_negative_something = 1 - dist_from_center_norm
        one_to_zero = np.maximum(one_to_negative_something, 0)
        return one_to_zero**(1.5)

        
    # dim_to_alpha = lambda dim: (np.array(dim) - np.min(dim)) / (np.max(dim) - np.min(dim))

    comx, comy, comz = 0, 0, 0
    all_x, all_y, all_z = pos_all(plummer)
    # member_x, member_y, member_z, hardnesses, ids = pos_members(plummer)

    # ax1.scatter(all_x, all_y, c='grey', marker=np.array([f'${iden}$' for iden in range(len(all_x))]), alpha=dim_to_alpha(all_z - comz))
    # ax2.scatter(all_z, all_y, c='grey', marker=np.array([f'${iden}$' for iden in range(len(all_x))]), alpha=dim_to_alpha(all_x - comz))
    for i, (x, y, z) in enumerate(zip(all_x, all_y, all_z)):
        ax1.scatter(x, y, c='grey', marker=f'${i}$', alpha=dim_to_alpha_centered(all_z, center=focus_z)[i])
        ax2.scatter(z, y, c='grey', marker=f'${i}$', alpha=dim_to_alpha_centered(all_x, center=focus_x)[i])
    
    ax1.scatter(comx, comy, c='red', alpha=dim_to_alpha_centered(all_z, center=focus_z), marker="*")
    ax2.scatter(comz, comy, c='red', alpha=dim_to_alpha_centered(all_x, center=focus_x), marker="*")
    
    # for i, (x, y, z, hardness, name) in enumerate(zip(*pos_members(plummer))):
    #     ax1.scatter(x, y, c=hardness, cmap='viridis', vmin=0, vmax=20, marker=f"${name}$")
    #     ax2.scatter(z, y, c=hardness, cmap='viridis', vmin=0, vmax=20, marker=f"${name}$")
        
    # ax1.scatter(member_x, member_y, c=hardnesses, cmap='viridis', vmin=0, vmax=20)
    # ax2.scatter(member_z, member_y, c=hardnesses, cmap='viridis', vmin=0, vmax=20)

    ax1.set_aspect('equal')
    ax2.set_aspect('equal')

    ax1.set_xlim(focus_x - viewsize, focus_x + viewsize)
    ax2.set_xlim(focus_z - viewsize, focus_z + viewsize)
    ax1.set_ylim(focus_y - viewsize, focus_y + viewsize)
    ax2.set_ylim(focus_y - viewsize, focus_y + viewsize)
    
    fig.suptitle(f'Time {time}')
    ax1.set_xlabel('x')
    ax2.set_xlabel('z')
    ax1.set_ylabel('y')
    ax2.set_ylabel('y')

    canvas.draw()


def get_focus_star_id():
    return int(focus_id_slider.get())
    
def get_view_size():
    return view_size_slider.get()

window = tk.Tk()
window.title('Decomp Scatter')
window.geometry('2000x1200')

plot_time_slider = tk.Scale(master = window,  
                     command = plot,
                     from_ = 190,
                     to = 200, 
                     orient="horizontal",
                     resolution = 0.01,
                     tickinterval=1,
                    #  height = 10,  
                    #  width = 400, 
                     length=1200,
                    #  text = "Time",
                     ) 

plot_time_slider.set(1)
plot_time_slider.pack()

focus_id_slider = tk.Scale(master=window,
                          command = plot,
                          from_ = -1,
                          to = N_STARS,
                          orient='horizontal',
                          resolution=1,
                          tickinterval=N_STARS/16,
                          length=1200)

focus_id_slider.set(0)
focus_id_slider.pack()

view_size_slider = tk.Scale(master=window,
                          command = plot,
                          from_ = .01,
                          to = 5,
                          orient='horizontal',
                          resolution=.01,
                          tickinterval=1,
                          length=1200)

view_size_slider.set(3)
view_size_slider.pack()


fig, (ax1, ax2)= plt.subplots(1, 2, figsize=(24,10), dpi=150, layout='constrained')

canvas = FigureCanvasTkAgg(fig, master=window)
    
canvas.get_tk_widget().pack()

toolbar = NavigationToolbar2Tk(canvas, window)
toolbar.update()
canvas.get_tk_widget().pack()

window.mainloop()


if __name__ == '__main__':
    ...


    