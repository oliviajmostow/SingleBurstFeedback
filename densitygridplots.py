
import numpy as np
import matplotlib.pyplot as plt
from readData.DataLoader import DataLoader
from matplotlib.colors import LogNorm
import pandas as pd
import shared_data
from m1_plot_values import earlydensities
from m1_plot_values import middensities
from m1_plot_values import latedensities
from m1_plot_values import plot_radius
def main():
    fig, ax = plt.subplots(2,4)
    for i in range(4):
        ax[0,i].plot(earlydensities[i], plot_radius, color='blue', label='early')
        ax[0,i].plot(middensities[i], plot_radius, color='black', label='middle')
        ax[0,i].plot(latedensities[i], plot_radius, color='orange', label='late')

    for i in range(4):
        ax[1,i].plot(earlydensities[i+4], plot_radius, color='blue', label='early')
        ax[1,i].plot(middensities[i+4], plot_radius, color='black', label='middle')
        ax[1,i].plot(latedensities[i+4], plot_radius, color='orange', label='late')
    for ax in ax.flat:
        ax.legend()
        ax.set_ylim((1000, latedensities[0].max()*10))
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.axvline(x=2.8*.038, linestyle='--', color='gray', label='2.8*softening')
        ax.set(xlabel='Radius (kpc)', ylabel='Density ($M_{\odot}$/kpc$^3$)')
        ax.label_outer()
    fig.savefig('grid_density.pdf', bbox_inches='tight')
    
    return


if __name__=='__main__':
	main()

