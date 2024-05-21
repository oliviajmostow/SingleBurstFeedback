
import numpy as np
import matplotlib.pyplot as plt
from readData.DataLoader import DataLoader
from matplotlib.colors import LogNorm
import matplotlib as mpl
import pandas as pd
import cd_function
import shared_data
def main():
        spacing = 15
        radius = np.logspace(-1.2,1, spacing)
        fig, ax = shared_data.set_plot_params()
        runs=['../m10_ic0_sg/sg_mass10-output','../sizetest_9/mb_mass10-output','../sizetest_8/mb_mass10-output','../sizetest_7/mb_mass10-output']
        plot_radius = 1/2*(radius[:-1] +radius[1:])
        density_sg = cd_function.calc_density(runs[0], spacing=spacing)
        density_9 = cd_function.calc_density(runs[1], spacing=spacing)
        density_8 = cd_function.calc_density(runs[2], spacing=spacing)
        density_7 = cd_function.calc_density(runs[3], spacing=spacing)
        #colorbar shenanigans
        norm = mpl.colors.Normalize(vmin=1, vmax=5)
        cmap = mpl.cm.ScalarMappable(norm=norm, cmap=mpl.cm.jet)
        cmap.set_array([])
        densities = [density_9, density_8, density_7]
        sizes = [1e7,2.5e7, 5e7]
        Z = [[0,0],[0,0]]
        levels = range(1, 6, 1)
        mymap = mpl.colors.LinearSegmentedColormap.from_list('mycolors',['blue','red'])
        CS3 = plt.contourf(Z, levels, cmap=mymap)
        ax.cla()
        y = densities[1:4]
        z = sizes[1:4]
        for y1,z1 in zip(y,z):
            r = (float(z1) - 1e7)/(4e7)
            g=0
            b=1-r
            ax.plot(plot_radius, y1, color=(r,g,b))
        fig.colorbar(CS3, ax=ax).set_label('Mass Expelled Per Burst ($10^7$$M_{\odot}$/h)', size=16)
        ax.plot(plot_radius, density_sg, color='black', label='steady')
        ax.set_ylim((1000, density_sg.max()*5))
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.axvline(x=.076, linestyle='--', color='black', label='softening length')
        ax.axvline(x=2.8*.076, linestyle='--', color='gray', label='2.8 x softening')
        plt.xlabel('Radius (kpc)')
        ax.legend()
        plt.ylabel('Density ($M_{\odot}$/kpc$^3$)')
        fig.savefig('density.pdf', bbox_inches='tight')
        return


if __name__=='__main__':
	main()
 
