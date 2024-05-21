
import numpy as np
import matplotlib.pyplot as plt
from readData.DataLoader import DataLoader
from scipy import stats
import shared_data
def main():
	runs = ['../mag1_red8_g2/sb_mass8-output', '../mag1_red2_g2/sb_mass8-output']
	cat = DataLoader(runs[0], 127, 1, ['Masses','Coordinates','SubhaloPos','Velocities'], sub_idx=0)
	#disp = np.sqrt(np.sum(np.square(cat['Velocities']), axis=1))
	velocities = cat['Velocities']*np.sqrt(cat.time)
	print(cat.time)
	particle_radius = np.sqrt(np.sum((np.square(cat['Coordinates']-cat['SubhaloPos'])), axis=1))
	std_x, b, bn = stats.binned_statistic(particle_radius, velocities[:,0],'std', bins=100)
	std_y, b, bn = stats.binned_statistic(particle_radius, velocities[:,1],'std', bins=100)
	std_z, b, bn = stats.binned_statistic(particle_radius, velocities[:,2],'std', bins=100)
	#binned stat function takes in the variable for the bins (radius), values (disp), and computes mean by default
	fig, ax = shared_data.set_plot_params()
	#std_velocity has the standard deviation for x, y, z; now need to take total magnitude
	dispersion = np.sqrt(np.square(std_x) + np.square(std_y) + np.square(std_z))
	plot_radius = 1/2*(b[0:-1] + b[1:])
	ax.plot(plot_radius, dispersion,color='black', label='burst at z=8')
	ax.set_xscale('log')
	plt.xlabel('Radius (kpc)')
	plt.ylabel('Velocity Dispersion (km/s)')
	cat = DataLoader(runs[1], 127, 1, ['Masses','Coordinates','SubhaloPos','Velocities'], sub_idx=0)
	velocities = cat['Velocities']*np.sqrt(cat.time)
	particle_radius = np.sqrt(np.sum((np.square(cat['Coordinates']-cat['SubhaloPos'])), axis=1))
	std_x, b, bn = stats.binned_statistic(particle_radius, velocities[:,0],'std', bins=100)
	std_y, b, bn = stats.binned_statistic(particle_radius, velocities[:,1],'std', bins=100)
	std_z, b, bn = stats.binned_statistic(particle_radius, velocities[:,2],'std', bins=100)
	dispersion = np.sqrt(np.square(std_x) + np.square(std_y) + np.square(std_z))
	plot_radius = 1/2*(b[0:-1] + b[1:])
	ax.plot(plot_radius, dispersion,color='orchid', label='burst at z=2')
	ax.legend()
	fig.savefig('veldisp.pdf', bbox_inches='tight')
	return


if __name__=='__main__':
	main()

