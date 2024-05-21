
import numpy as np
import matplotlib.pyplot as plt
from readData.DataLoader import DataLoader
from scipy import stats
import shared_data
import cd_function
def main():
        filepath='../sizetest_7/mb_mass10-output'
        spacing = 20
        vd = cd_function.calc_vd(filepath, spacing)
        ax.plot(vd.plot_radius, vd.dispersion, color='black', label='s7')
        ax.set_xscale('log')
#	ax.set_xlim((10**-1), (10))
        plt.xlabel('Radius (kpc)')
        ax.legend()
#	ax.set_ylim(0,50)
        plt.ylabel('Velocity Dispersion (km/s)')
        fig.savefig('veldisp.pdf', bbox_inches='tight')
        return


if __name__=='__main__':
	main()

