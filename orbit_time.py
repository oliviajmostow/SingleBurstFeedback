
import numpy as np
import matplotlib.pyplot as plt
from readData.DataLoader import DataLoader
from matplotlib.colors import LogNorm
import pandas as pd
plt.rcParams.update({'font.size': 14})
def main():
	cat = DataLoader('../eighthrun/eighth-output', 49, 1, ['Masses','ParticleIDs','Coordinates','GroupPos','SubhaloMass','SubhaloPos', 'GroupMassType'], sub_idx=0)
	particle_radius = np.sqrt(np.sum((np.square(cat['Coordinates'] - cat['SubhaloPos'])),axis=1))
	print(cat.time)
	radius = np.logspace(-1,1, 21)
	hist, bin_edges = np.histogram(particle_radius, radius)
	print(bin_edges)
	ingroup1 = particle_radius < bin_edges[0]
	ingroup2 = (particle_radius > bin_edges[1]) & (particle_radius < bin_edges[2])
	ingroup3 = (particle_radius > bin_edges[2]) & (particle_radius < bin_edges[3])
	ingroup4 = (particle_radius > bin_edges[3]) & (particle_radius < bin_edges[4])
	ingroup5 = (particle_radius > bin_edges[4]) & (particle_radius < bin_edges[5])
	ingroup6 = (particle_radius > bin_edges[5]) & (particle_radius < bin_edges[6])
	ingroup7 = (particle_radius > bin_edges[6]) & (particle_radius < bin_edges[7])
	ingroup8 = (particle_radius > bin_edges[7]) & (particle_radius < bin_edges[8])
	ingroup9 = (particle_radius > bin_edges[8]) & (particle_radius < bin_edges[9])
	ingroup10 = (particle_radius > bin_edges[9]) & (particle_radius < bin_edges[10])
	ingroup11 = (particle_radius > bin_edges[10]) & (particle_radius < bin_edges[11])
	ingroup12 = (particle_radius > bin_edges[11]) & (particle_radius < bin_edges[12])
	ingroup13 = (particle_radius > bin_edges[12]) & (particle_radius < bin_edges[13])
	ingroup14 = (particle_radius > bin_edges[13]) & (particle_radius < bin_edges[14])
	ingroup15 = (particle_radius > bin_edges[14]) & (particle_radius < bin_edges[15])
	ingroup16 = (particle_radius > bin_edges[15]) & (particle_radius < bin_edges[16])
	ingroup17 = (particle_radius > bin_edges[16])& (particle_radius < bin_edges[17])
	ingroup18 = (particle_radius > bin_edges[17])&(particle_radius < bin_edges[18])
	ingroup19 = (particle_radius > bin_edges[18])&(particle_radius < bin_edges[19])
	ingroup20 = (particle_radius > bin_edges[19])&(particle_radius < bin_edges[20])
	group_list = [ingroup1, ingroup2, ingroup3, ingroup4, ingroup5, ingroup6, ingroup7, ingroup8, ingroup9, ingroup10, ingroup11, ingroup12, ingroup13, ingroup14, ingroup15, ingroup16, ingroup17, ingroup18, ingroup19, ingroup20]
	empty_lists = [ [] for x in range(20) ]
	group_averages = [[] for x in range(20)]
	particle_lists_list = []
	print(group_list)
	for i in group_list:
		particles = cat['ParticleIDs'][i]
		particle_lists_list.append(particles)
	for i in range(50, 128):
		print(i)
		cat = DataLoader('../stop5_hf_m10/mb_mass10-output', i, 1, ['Masses','ParticleIDs','Coordinates','GroupPos','SubhaloMass','SubhaloPos', 'GroupMassType'], sub_idx=0)
		particle_radius = np.sqrt(np.sum((np.square(cat['Coordinates'] - cat['SubhaloPos'])),axis=1))
		for j in range(len(cat['ParticleIDs'])):
			for z in range(len(particle_lists_list)):
				if cat['ParticleIDs'][j] in particle_lists_list[z]:
					empty_lists[z].append(particle_radius[j])
		for k in range(len(empty_lists)):
			x = empty_lists[k].mean()
			group_averages[k].append(empty_lists[k])
	print(group_averages[0])
			
	return


if __name__=='__main__':
	main()

