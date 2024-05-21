
import numpy as np
import matplotlib.pyplot as plt
from readData.DataLoader import DataLoader
from matplotlib.colors import LogNorm
import pandas as pd
import shared_data
def main():
	cat = DataLoader('../tenthrun/tenth-output', 74, 1, ['Masses','Coordinates','Potential','ParticleIDs'], fof_idx=0)
	print(cat.time)
	x = cat['Coordinates'][:,0]
	y = cat['Coordinates'][:,1]
	center_potential = cat['Potential'].min()
	fig, ax = plt.subplots()
	ax.hist2d(x, y, bins=50, norm=LogNorm())
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
	inside = df[df.radial_distance<=100]
	x_coord = inside.x_coord
	y_coord = inside.y_coord
#	h=ax.hist2d(x_coord, y_coord, bins=100, norm=LogNorm())
#	plt.colorbar(h[3], ax=ax, label='density')
#	ax.plot(x_0, y_0, 'o')
	ax.set_xticks([])
	ax.set_yticks([])
	fig.savefig('img.pdf')
	return


if __name__=='__main__':
	main()

