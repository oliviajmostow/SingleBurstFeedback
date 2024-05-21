
import numpy as np
import matplotlib.pyplot as plt
from readData.DataLoader import DataLoader
from  astropy.cosmology import Planck13
import pandas as pd
import shared_data
def main():
	seed_mass = []
	halo_mass = []
	minsnap=32
	early_indexes = [0,1,1,0,0,0,0,0,0,0,0,0,1,1,1,0,1]
#	for i in range(22, 39):
#		cat = DataLoader('../sg_m8/sg_mass8-output', i, 5, ['Masses'], fof_idx=early_indexes[i-22])
#		seed_mass.append(cat['PartType5/Masses'][0])
	for i in range(127, 128):
                cat = DataLoader('../m10_destructive_run/mb_mass10-output', i, [1,5], ['Masses'], sub_idx=0)
                print(i)
                seed_mass.append(cat['PartType5/Masses'][0])
                print(cat['PartType5/Masses'][0])
	df = pd.read_csv('../ExpansionList_1', sep=' ', header=None)
	seed_mass = np.array(seed_mass)
	seed_mass = seed_mass * (1e10)
	for i in range(0, 128):
		cat = DataLoader('../717_m10_hf_stop5/mb_mass10-output', i, 1, ['Masses'],fof_idx=0)
		halo_mass.append(np.sum(cat['PartType1/Masses']))
	halo_mass =np.array(halo_mass) * (1e10)
#	for i in range(45,128):
#		cat2 = DataLoader('../fourthrun/fourth-output', i,1, ['Masses', 'GroupMassType'])
#		print(cat2['GroupMassType'][0,5])
#		if cat2['GroupMassType'][0,5] == 0:
#			print(i)
	
#	want the corresponding lookback times for each snapshot, can do based on redshift or probably scale factor too
	scale_factors = df[2].array
	z = 1/scale_factors - 1
	z2 = z[minsnap:]
	cosmo = Planck13
	t = cosmo.lookback_time(z)
	t2 = cosmo.lookback_time(z2)
#saving params to use in shmr
	M_10 = 11.590
	M_11 = 1.195
	N_10 = .0351
	N_11 = -.0247
	B_10 = 1.376
	B_11 = -.826
	G_10 = .608
	G_11 = .329
#redshift-dependent params
	x = z/(z+1)
	log_Mz = M_10 + M_11 * x
	M = 10**log_Mz
	N = N_10 + N_11 * x
	B = B_10 + B_11 * x
	G = G_10 + G_11 * x
	stellar_mass = halo_mass * 2 * N* np.power((halo_mass/M)**(-B) + (halo_mass/M)**(G), -1)
	fig, ax = shared_data.set_plot_params()
	ax.plot(t2, seed_mass, color='blue')
#	ax.plot(t, halo_mass)
#	ax.plot(t, stellar_mass)
	ax2 = ax.twiny()
	ax.set_xlim(0, 14)
	newlabel = [12, 6, 4, 3, 2, 1]
	convert = lambda x : cosmo.lookback_time(x).value
	newpos = [convert(z) for z in newlabel]
	ax2.set_xticks(newpos)
	ax2.set_xticklabels(newlabel)
	ax.invert_xaxis()
	ax2.invert_xaxis()
	ax2.set_xlabel('Redshift')
	ax2.set_xlim(ax.get_xlim())
#	ax.set_yscale('log')
	ax.set_xlabel('Lookback Time (Gyr)')
	ax.set_ylabel('Baryon Mass')
	fig.savefig('mass_time.pdf', bbox_inches='tight')
	fig, ax = plt.subplots()
	ax.plot(halo_mass, stellar_mass)
	ax.set_yscale('log')
	ax.set_xscale('log')
	ax.set_xlabel('Halo Mass ($M_{\odot}$)')
	ax.set_ylabel('Stellar Mass ($M_{\odot}$)')
	fig.savefig('shmr.pdf')
	fig, ax = plt.subplots()
	ax.plot(scale_factors, stellar_mass)
	ax.set_xlabel('Scale Factor')
	ax.set_ylabel('Stellar Mass ($M_{\odot}$)')
	ax.set_yscale('log')
	fig.savefig('sf_sm.pdf')
	fig, ax = shared_data.set_plot_params()
	ax.plot(t2, seed_mass)
	ax.set_yscale('log')
	fig.savefig('seed.pdf')
	return


if __name__=='__main__':
	main()

