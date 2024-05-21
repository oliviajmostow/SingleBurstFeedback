import numpy as np
import matplotlib.pyplot as plt
from readData.DataLoader import DataLoader
from matplotlib.colors import LogNorm
import pandas as pd
import shared_data
from  astropy.cosmology import FlatLambdaCDM
plt.rcParams.update({'font.size': 14})
def main():
	alpha_list = []
	spacing = 20
	upper = 2
	lower = 0
	for i in range(45, 128):
		cat = DataLoader('../mag2_red2_g2_up/sb_mass8-output', i, 1, ['Masses','Coordinates','SubhaloPos'], sub_idx=0)
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
	last_Gyr = t.value<=1
	final = np.array(alpha_list)[last_Gyr]
	final_med = np.median(final)
	print(final_med)
	fig, ax = shared_data.set_plot_params()
#	ax.plot(scale_factors, alpha_list, color='blueviolet', label='high frequency')
	burst_times = [.0526, .058, .063, .069, .0769, .085, .092, .102, .111, .120, .134, .145, .160, .170, .190, .210, .225, .245, .265, .295, .320, .350, .380, .410, .440, .470, .50, .55, .6, .65, .7, .75, .8, .85, .9, .95, 1]
	lf_bursts = [.058, .069, .085, .102, .120, .145, .170, .210, .245, .295, .350, .410, .470, .55, .65, .75, .85,.95]
	lf_bursts_shifted = np.array(lf_bursts) + .01
	new_hf_bursts = [.051, .056, .063, .069, .076, .084, .093, .103, .113, .125, .138, .154, .168, .186, .207, .228, .251, .280, .307, .336, .381, .422, .485]
	ax.plot(t, alpha_list, color='dodgerblue')
	ax.axvline(x=cosmo.lookback_time(1).value, color='silver', ls = '--', label='burst')
	ax.set_xlim(0, 13)
	ax2 = ax.twiny()
	newlabel = [12, 6, 4, 3, 2, 1]
	convert = lambda x : cosmo.lookback_time(x).value
	newpos = [convert(z) for z in newlabel]
	ax2.set_xticks(newpos)
	ax2.set_xticklabels(newlabel)
	ax.invert_xaxis()
	ax2.invert_xaxis()
#	for i in new_hf_bursts:
#		ax.axvline(x=cosmo.lookback_time(1/i-1).value, color='silver', ls='--')
	ax2.set_xlabel('Redshift')
	ax2.set_xlim(ax.get_xlim())
	ax.set_xlabel('Time (Gyr)')
	ax.set_ylabel('Inner Slope')
	ax.plot(t, alpha_list, color='dodgerblue')
#	ax.set_ylim(-2, .5)
	ax.legend()
	print(radius[upper])
	print(radius[lower])
	fig.savefig('alpha_time.pdf', bbox_inches='tight')
	return


if __name__=='__main__':
	main()

