
import numpy as np
import matplotlib.pyplot as plt
from readData.DataLoader import DataLoader
from  astropy.cosmology import Planck13, z_at_value
import pandas as pd
import shared_data
from astropy import units as u
plt.rcParams.update({'font.size': 12})
def main():
	halo_mass = []
	for i in range(39, 128):
		cat = DataLoader('../m8_dmo_g2/dmo_mass8-output', i, 1, ['Masses'], sub_idx=0)
		halo_mass.append(np.sum(cat['PartType1/Masses']))
	df = pd.read_csv('../ExpansionList_1', sep=' ', header=None)
	halo_mass = np.array(halo_mass)
	halo_mass = halo_mass * (1e10)
	halo_mass_hr = []
	for i in range(39,128):
		cat = DataLoader('../m8_dmo_hr/dmo_mass8-output', i, 1, ['Masses'], sub_idx=0)
		halo_mass_hr.append(np.sum(cat['PartType1/Masses']))
	halo_mass_hr = np.array(halo_mass_hr)
	halo_mass_hr = halo_mass_hr*1e10
	
#	want the corresponding lookback times for each snapshot, can do based on redshift or probably scale factor too
	scale_factors = df[2].array
	z = (1/scale_factors - 1)[39:]
	print(len(z))
	cosmo = Planck13
	t = cosmo.lookback_time(z)

	cmf = []
	for i in range(0, 89):
		fractional_mass = halo_mass[i]/halo_mass[-1]
		cmf.append(fractional_mass)

	fig, ax = shared_data.set_plot_params()
	ax.plot(t,halo_mass, color='mediumseagreen', label='galaxy 1')
	ax.plot(t, halo_mass_hr, color='black', label='galaxy 2')
	ax.set_xlabel('Lookback Time (Gyr)')
	ax.set_ylabel('Halo Mass')
	ax2 = ax.twiny()
	newlabel = [12, 6, 4, 3, 2, 1]
	#convert = lambda x: z_at_value(Planck13.age, x)
	convert = lambda x : cosmo.lookback_time(x).value
	newpos = [convert(z) for z in newlabel]
	ax2.set_xticks(newpos)
	ax2.set_xticklabels(newlabel)
	ax.invert_xaxis()
	ax.set_yscale('log')
	ax2.invert_xaxis()
	ax.legend()
	ax2.set_xlabel('Redshift')
	ax2.set_xlim(ax.get_xlim())
	fig.savefig('halomass.pdf', bbox_inches='tight')
		

	return


if __name__=='__main__':
	main()

