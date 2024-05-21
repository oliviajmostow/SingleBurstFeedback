import numpy as np
import matplotlib.pyplot as plt
from readData.DataLoader import DataLoader
from matplotlib.colors import LogNorm
import pandas as pd
from  astropy.cosmology import FlatLambdaCDM, z_at_value
import shared_data
plt.rcParams.update({'font.size': 14})
def main():
        alpha_list = []
        spacing = 15
        upper = 2
        lower = 0
        for i in range(45, 128):
                cat = DataLoader('../ft_s8_20/mb_mass10-output', i, 1, ['Masses','Coordinates','SubhaloPos'], sub_idx=0)
                particle_radius = np.sqrt(np.sum((np.square(cat['Coordinates']-cat['SubhaloPos'])), axis=1))
                radius = np.logspace(-1,1,spacing)
                radius_normal = np.logspace(-1, 1, spacing)
                radius_shifted = radius_normal[1:]
                plot_radius = 1/2*(radius_normal[0:-1] + radius_shifted)
                radius2 = radius[1:]
                volume = (4/3) * np.pi * (radius)**3
                volume2 =(4/3) * np.pi * (radius2)**3
                shell_volume  = volume2 -  volume[0:-1]
                masses = cat['Masses'] * (1e10)
                hist, bin_edges = np.histogram(particle_radius, radius, weights=masses)
                density = hist/shell_volume
                logp = np.log10(density)
                logr = np.log10(radius)
                dlogr = logr[upper] - logr[lower]
                dlogp = logp[upper] - logp[lower]
                alpha = dlogp/dlogr
                alpha_list.append(alpha)
        df = pd.read_csv('../ExpansionList_1', sep=' ', header=None)
        scale_factors = df[2].array[45:]
        z = (1/scale_factors - 1)
        cosmo = FlatLambdaCDM(H0=69.09, Om0=.3, Tcmb0=2.725)
        t = cosmo.lookback_time(z)
        fig, ax = shared_data.set_plot_params()
#	ax.plot(scale_factors, alpha_list, color='blueviolet', label='high frequency')
        burst_times = [.0526, .058, .063, .069, .0769, .085, .092, .102, .111, .120, .134, .145, .160, .170, .190, .210, .225, .245, .265, .295, .320, .350, .380, .410, .440, .470, .50, .55, .6, .65, .7, .75, .8, .85, .9, .95, 1]
        lf_bursts = [.058, .069, .085, .102, .120, .145, .170, .210, .245, .295, .350, .410, .470, .55, .65, .75, .85,.95]
        burst_times_30 = [.085, .092, .102, .111, .120, .134, .145, .160, .170, .190, .210, .225, .245, .265, .295, .320, .350, .380, .410, .440, .470, .50, .55, .6, .65, .7, .75, .8, .85, .9, .95]
        burst_times_20 = [.092, .111, .120, .145, .160, .190, .210, .245, .265, .320, .350, .410, .440, .50, .55, .65, .7, .8, .85, .95]
        burst_times_10 = [.092, .120, .160, .210, .265, .350, .440, .55, .7, .85]
        for i in burst_times_20:
            ax.axvline(x=(i), color='silver', ls='--')
        ax.plot(scale_factors, alpha_list, color='black')
        ax.set_xlabel('Scale Factor')
        ax.set_ylabel('Inner Slope')
        ax.set_xlim(.35, 1)
        ax.set_ylim(-2, 1.5)
        ax.legend()
        fig.savefig('alpha_time.pdf', bbox_inches='tight')
        return


if __name__=='__main__':
	main()

