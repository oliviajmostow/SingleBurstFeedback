
import numpy as np
import matplotlib.pyplot as plt
from readData.DataLoader import DataLoader
from  astropy.cosmology import FlatLambdaCDM, z_at_value
import pandas as pd
import shared_data
from astropy import units as u
plt.rcParams.update({'font.size': 12})
def main():
        halo_mass = []
        early_indexes = [0,1,1,0,0,0,0,0,0,0,0,0,1,1,1,0,1]
#	for i in range(0,22):
#		cat = DataLoader('../eighthrun/eighth-output', i, 1, ['Masses'], fof_idx=0)
#		halo_mass.append(np.sum(cat['PartType1/Masses']))
#	for i in range(22, 39):
#		cat = DataLoader('../eighthrun/eighth-output', i, 1, ['Masses'], fof_idx=early_indexes[i-22])
#		halo_mass.append(np.sum(cat['PartType1/Masses']))
        for i in range(30, 128):
            cat = DataLoader('../m10_dmo/dmo_mass10-output', i, 1, ['Masses'], fof_idx=0)
            halo_mass.append(np.sum(cat['PartType1/Masses']))
        df = pd.read_csv('../ExpansionList_1', sep=' ', header=None)
        halo_mass = np.array(halo_mass)
        halo_mass = halo_mass * (1e10)
        print(halo_mass[-1])
#	for i in range(45,128):
#		cat2 = DataLoader('../fourthrun/fourth-output', i,1, ['Masses', 'GroupMassType'])
#		print(cat2['GroupMassType'][0,5])
#		if cat2['GroupMassType'][0,5] == 0:
#			print(i)
	
#	want the corresponding lookback times for each snapshot, can do based on redshift or probably scale factor too
        scale_factors = df[2].array
        z = (1/scale_factors - 1)[30:]
#	z2 = z[45:]
        print(len(z))
        cosmo = FlatLambdaCDM(H0=69.09, Om0=.3, Tcmb0=2.725)
        t = cosmo.lookback_time(z)
#	t2 = cosmo.lookback_time(z2)
#saving params to use in shmr
        M_10 = 11.590
        M_11 = 1.195
        N_10 = .0351
        N_11 = -.0247
        B_10 = 1.62
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
        print(stellar_mass[-1])
#	fig, ax = plt.subplots()
#	print(halo_mass[0])
#	print(halo_mass[-1])
#	print(stellar_mass[-1])
	#ax.plot(scale_factors, halo_mass)
#	ax.set_yscale('log')
#	ax.set_xscale('log')
#	ax.set_xlabel('Scale Factor ($M_{\odot}$)')
#	ax.set_ylabel('Stellar Mass ($M_{\odot}$)')
#	fig.savefig('stellarmass.pdf')
        cmf = []
        for i in range(0, 98):
            fractional_mass = halo_mass[i]/halo_mass[-1]
            cmf.append(fractional_mass)
        fig, ax = shared_data.set_plot_params()
        ax.plot(t, cmf, color='dodgerblue')
        ax.set_xlabel('Lookback Time (Gyr)')
        ax.set_ylabel('Halo CMF')
        ax2 = ax.twiny()
        newlabel = [12, 6, 4, 3, 2, 1]
	#convert = lambda x: z_at_value(Planck13.age, x)
        convert = lambda x : cosmo.lookback_time(x).value
        newpos = [convert(z) for z in newlabel]
        ax2.set_xticks(newpos)
        ax2.set_xticklabels(newlabel)
        ax.invert_xaxis()
        ax2.invert_xaxis()
        ax2.set_xlabel('Redshift')
        ax2.set_xlim(ax.get_xlim())
#	ax.axvline(x=cosmo.lookback_time(7).value, ls='--', color='gray')
#	ax.axvline(x=cosmo.lookback_time(3).value, ls='--', color='gray')
#	ax.axvline(x=cosmo.lookback_time(2).value, ls='--', color='gray')
        fig.savefig('halocmf.pdf', bbox_inches='tight')
		
        return


if __name__=='__main__':
	main()

