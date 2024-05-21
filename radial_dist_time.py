
import numpy as np
import matplotlib.pyplot as plt
from readData.DataLoader import DataLoader
from matplotlib.colors import LogNorm
import pandas as pd
from astropy.cosmology import Planck13
plt.rcParams.update({'font.size': 14})
def main():
	cat = DataLoader('../stop5_hf_m10/mb_mass10-output', 49, 1, ['Masses','ParticleIDs','Coordinates','GroupPos','SubhaloMass','SubhaloPos', 'GroupMassType'], sub_idx=0)
	particle_radius = np.sqrt(np.sum((np.square(cat['Coordinates'] - cat['SubhaloPos'])),axis=1))
	print(cat.time)
	radius = np.logspace(-1,1, 21)
	hist, bin_edges = np.histogram(particle_radius, radius)
	print(bin_edges)
	inner = particle_radius < 1
	inner_particles = cat['ParticleIDs'][inner]
	avg_radii = []
	radius_time = []
	sf = []
	for i in range(50, 128):
		cat = DataLoader('../stop5_hf_m10/mb_mass10-output', i, 1, ['Masses','ParticleIDs','Coordinates','GroupPos','SubhaloMass','SubhaloPos', 'GroupMassType'], sub_idx=0)
		particle_radius = np.sqrt(np.sum((np.square(cat['Coordinates'] - cat['SubhaloPos'])),axis=1))
		for j in range(len(cat['ParticleIDs'])):
			if cat['ParticleIDs'][j] in inner_particles:
				radius_time.append(particle_radius[j])
		avg_radii.append(np.array(radius_time).mean())
		sf.append(cat.time)
	fig, ax = plt.subplots()
	sf = np.array(sf)
	cosmo = Planck13
	t = cosmo.lookback_time((1/sf - 1))
	ax.plot(t, avg_radii)
	ax.invert_xaxis()
	ax.set_xlabel('t/Gyr')
	ax.set_ylabel('Average Radius')
	fig.savefig('orbit_time.pdf')

	return


if __name__=='__main__':
	main()

