
import numpy as np
import matplotlib.pyplot as plt
from readData.DataLoader import DataLoader
from matplotlib.colors import LogNorm
import pandas as pd
plt.rcParams.update({'font.size': 14})
def main():
	cat = DataLoader('../sg_m8_hr/sg_mass8-output', 126, 1, ['Coordinates','SubhaloPos'], sub_idx=0)
	sub = cat['Coordinates'] - cat['SubhaloPos']
	sub2 = np.square(sub)
	sub_sum = np.sum(sub2, axis=1)
	particle_radius = np.sqrt(sub_sum)
	radius = np.logspace(-1,1,10)
	radius_normal = np.logspace(-1, 1,10)
	radius_shifted = radius_normal[1:]
	plot_radius = 1/2*(radius_normal[0:-1] + radius_shifted)
	radius2 = radius[1:]
	volume = (4/3) * np.pi * (radius)**3

	volume2 =(4/3) * np.pi * (radius2)**3
	shell_volume  = volume2 -  volume[0:-1]
	volume_p1 = (4/3) * np.pi * (radius[0])**3
	vol_list =  list(shell_volume)
	shell_volume = np.array(vol_list)
	hist, bin_edges = np.histogram(particle_radius, radius)
	fig, ax = plt.subplots()
	ax.hist(particle_radius, bins=radius, fill=False)
	plt.xlabel('Radius (kpc)')
	plt.ylabel('Particles')
	ax.set_xscale('log')
	ax.set_xlim(.1, 3)
	ax.set_ylim(0, 1000)
	ax.axvline(x=.22, ls='--')
	fig.savefig('radius_hist.pdf', bbox_inches='tight')
	return


if __name__=='__main__':
	main()

