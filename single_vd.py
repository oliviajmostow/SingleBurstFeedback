
import numpy as np
import matplotlib.pyplot as plt
from readData.DataLoader import DataLoader
from scipy import stats
import shared_data
def main():
	time = 127
	colors = ['mediumblue', 'dodgerblue', 'lightskyblue', 'pink', 'red','black']
	labels = ['z=10','z=8','z=6','z=4','z=2','sg','multiple bursts','1e10 halo, dmo']
	num_bins = 50
	i = -1
	runs= ['../mag3_red10_hr/sb_mass8-output','../mag3_red8_hr/sb_mass8-output', '../mag3_red6_hr/sb_mass8-output', '../mag3_red4_hr/sb_mass8-output', '../mag3_red2_hr/sb_mass8-output', '../sg_m8_hr/sg_mass8-output','../717_m10_hf_stop5/mb_mass10-output', '../eighthrun/eighth-output']
	fig, ax = shared_data.set_plot_params()
	cat = DataLoader(runs[i],time, 1, ['Masses','Coordinates','SubhaloPos','Velocities'], sub_idx=0)
	#disp = np.sqrt(np.sum(np.square(cat['Velocities']), axis=1))
	velocities = cat['Velocities']*np.sqrt(cat.time)
	particle_radius = np.sqrt(np.sum((np.square(cat['Coordinates']-cat['SubhaloPos'])), axis=1))*cat.time
	std_x, b, bn = stats.binned_statistic(particle_radius, velocities[:,0],'std', bins=num_bins)
	std_y, b, bn = stats.binned_statistic(particle_radius, velocities[:,1],'std', bins=num_bins)
	std_z, b, bn = stats.binned_statistic(particle_radius, velocities[:,2],'std', bins=num_bins)
	#binned stat function takes in the variable for the bins (radius), values (disp), and computes mean by default
	#std_velocity has the standard deviation for x, y, z; now need to take total magnitude
	dispersion = np.sqrt(np.square(std_x) + np.square(std_y) + np.square(std_z))
	plot_radius = 1/2*(b[0:-1] + b[1:])
	print(cat.time)
	ax.plot(plot_radius, dispersion,color=colors[i], label=labels[i])
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

