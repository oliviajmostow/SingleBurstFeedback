
import numpy as np
import matplotlib.pyplot as plt
from readData.DataLoader import DataLoader
from matplotlib.colors import LogNorm
import pandas as pd
def main():
	cat = DataLoader('../717_hf_stop5/mb_mass10-output', 127, 1, ['Masses','Coordinates','Potential','ParticleIDs'])
	x = cat['Coordinates'][:,0]
	y = cat['Coordinates'][:,1]
	center_potential = cat['Potential'].min()
	fig, ax = plt.subplots()
#	ax.hist2d(x, y, bins=100, norm=LogNorm())
	print(center_potential)
	df = pd.DataFrame()
	df['x_coord']=cat['Coordinates'][:,0]
	df['y_coord']=cat['Coordinates'][:,1]
	df['z_coord']=cat['Coordinates'][:,2]
	df['Potential']=cat['Potential']	
	df['Masses']=cat['Masses']
	df['id']=cat['ParticleIDs']
	#density associated with a given particle's surronding sphere w/ a fixed number of particles- need shell density of a given radius, groupby radial distance?
	df = df.set_index('id')
	center_data = df[df.Potential==center_potential]
	#ity'rint(df.head())
	#print(center_data)
	#print('(' + str(center_data.x_coord.loc[25340.0]) + ',' + str(center_data.y_coord.loc[25340.0]) + ',' + str(center_data.z_coord.loc[25340.0]) + ')')
	x_0 = center_data.x_coord.iloc[0]
	y_0 = center_data.y_coord.iloc[0]
	z_0 = center_data.z_coord.iloc[0]
	df['radial_distance'] = ((df['x_coord']-x_0)**2+(df['y_coord']-y_0)**2+(df['z_coord']-z_0)**2)**.5
#	ax.plot(x_0, y_0, 'o')
	cat2 = DataLoader('../fourthrun/fourth-output', 95, 5, ['Coordinates'], sub_idx=0)
	print(cat2['Coordinates'])
	bh_x = cat2['Coordinates'][0,0]
	bh_y = cat2['Coordinates'][0,1]
	rel_pos = cat['Coordinates'] - cat2['Coordinates']
	bh_radius = np.sqrt(np.sum((np.square(rel_pos)), axis=1))
	cut = bh_radius <= 50
	inside = cat['Coordinates'][cut]
	x_coord = inside[:,0]
	y_coord = inside[:,1]
	ax.hist2d(x_coord, y_coord, bins=100, norm=LogNorm())
	ax.plot(bh_x, bh_y, 'o')
#	ax.set_xlim(47800, 48200)
#	ax.set_ylim(51000, 52000)
	fig.savefig('img.pdf')
	return


if __name__=='__main__':
	main()

