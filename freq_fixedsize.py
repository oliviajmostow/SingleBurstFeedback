
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
        fig, ax = shared_data.set_plot_params(nrows=1, ncols=2)
        runs=['../m10_ic0_sg/sg_mass10-output','../sizetest_9/mb_mass10-output','../ft_s9_20/mb_mass10-output','../ft_s9_10/mb_mass10-output','../ft_s9_5/mb_mass10-output','../ft_s9_1/sb_mass10-output']
#        runs=['../m10_ic0_sg/sg_mass10-output','../sizetest_8/mb_mass10-output','../ft_s8_20/mb_mass10-output','../ft_s8_10/mb_mass10-output','../ft_s8_5/mb_mass10-output','../ft_s8_1/sb_mass10-output']
        plot_radius = 1/2*(radius[:-1] +radius[1:])
        density_sg = cd_function.calc_density(runs[0], spacing=spacing)
        density_30 = cd_function.calc_density(runs[1], spacing=spacing)
        density_20 = cd_function.calc_density(runs[2], spacing=spacing)
        density_10 = cd_function.calc_density(runs[3], spacing=spacing)
        density_5 = cd_function.calc_density(runs[4], spacing=spacing)
        density_1 = cd_function.calc_density(runs[5], spacing=spacing)
        #colorbar shenanigans
        norm = mpl.colors.Normalize(vmin=1, vmax=5)
        cmap = mpl.cm.ScalarMappable(norm=norm, cmap=mpl.cm.jet)
        cmap.set_array([])
        densities = [density_30, density_20, density_10, density_5, density_1]
        burst_nums = [31, 20, 10, 5, 1]
        Z = [[0,0],[0,0]]
        levels = range(1, 32, 5)
        mymap = mpl.colors.LinearSegmentedColormap.from_list('mycolors',['blue','red'])
        CS3 = plt.contourf(Z, levels, cmap=mymap)
        ax[1].cla()
        y = densities
        z = burst_nums
        for y1,z1 in zip(y,z):
            r = (float(z1) - 1)/(30)
            g=0
            b=1-r
            ax[0].plot(plot_radius, y1, color=(r,g,b))
        fig.colorbar(CS3, ax=ax[1]).set_label('Number of Bursts', size=16)
        ax[0].plot(plot_radius, density_sg, color='black', label='steady')
        ax[0].set_ylim((1000, density_sg.max()*5))
        ax[0].set_xscale('log')
        ax[0].set_yscale('log')
        ax[0].axvline(x=.076, linestyle='--', color='black', label='softening length')
        ax[0].axvline(x=2.8*.076, linestyle='--', color='gray', label='2.8 x softening')



        runs=['../m10_ic0_sg/sg_mass10-output','../sizetest_8/mb_mass10-output','../ft_s8_20/mb_mass10-output','../ft_s8_10/mb_mass10-output','../ft_s8_5/mb_mass10-output','../ft_s8_1/sb_mass10-output']
        plot_radius = 1/2*(radius[:-1] +radius[1:])
        density_sg = cd_function.calc_density(runs[0], spacing=spacing)
        density_30 = cd_function.calc_density(runs[1], spacing=spacing)
        density_20 = cd_function.calc_density(runs[2], spacing=spacing)
        density_10 = cd_function.calc_density(runs[3], spacing=spacing)
        density_5 = cd_function.calc_density(runs[4], spacing=spacing)
        density_1 = cd_function.calc_density(runs[5], spacing=spacing)
        #colorbar shenanigans
        densities = [density_30, density_20, density_10, density_5, density_1]
        burst_nums = [31, 20, 10, 5, 1]
        y = densities
        z = burst_nums
        for y1,z1 in zip(y,z):
            r = (float(z1) - 1)/(30)
            g=0
            b=1-r
            ax[1].plot(plot_radius, y1, color=(r,g,b))
        ax[1].plot(plot_radius, density_sg, color='black', label='steady')
        ax[1].set_ylim((1000, density_sg.max()*5))
        ax[1].set_xscale('log')
        ax[1].set_yscale('log')
        ax[1].axvline(x=.076, linestyle='--', color='black', label='softening length')
        ax[1].axvline(x=2.8*.076, linestyle='--', color='gray', label='2.8 x softening')
        ax[1].legend()
        ax[0].set_title('$10^7$$M_{\odot}$ Expelled Per Burst')
        ax[1].set_title('$2.5*10^7$$M_{\odot}$ Expelled Per Burst')
        ax[0].set_xlabel('Radius (kpc)')
        ax[0].set_ylabel('Density ($M_{\odot}$/kpc$^3$)')
        fig.savefig('density.pdf', bbox_inches='tight')
        return


if __name__=='__main__':
	main()
 
