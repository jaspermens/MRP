import matplotlib.pyplot as plt

plt.style.use('~/Templates/not_ugly.mplstyle')

def color_for_n(n_stars: int):
    match n_stars:
        case 16: return 'red'
        case 14: return 'firebrick'
        case 12: return 'maroon'
        case 32: return 'black'
        case 64: return 'xkcd:marine blue'
        case 72: return 'xkcd:electric blue'
        case 80: return 'xkcd:azure'
        case 128: return 'xkcd:aqua blue'
        case _: return 'mediumorchid'
    