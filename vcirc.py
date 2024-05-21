
import numpy as np
import matplotlib.pyplot as plt
from readData.DataLoader import DataLoader
from matplotlib.colors import LogNorm
import pandas as pd
from halotools.empirical_models import NFWProfile
import shared_data
def main():
	cat = DataLoader('../sg_m8/sg_mass8-output',127, 1, ['Masses','Coordinates','GroupPos','SubhaloMass','SubhaloPos', 'GroupMassType'], sub_idx=0)
	sub = cat['Coordinates'] - cat['SubhaloPos']
	print('Scale Factor = ' + str(cat.time))
	sub2 = np.square(sub)
	sub_sum = np.sum(sub2, axis=1)
	particle_radius = np.sqrt(sub_sum)
#remove particles with radius 0 to avoid errors dividing
	nonzero = particle_radius != 0
	particle_radius2 = particle_radius[nonzero]
#conversion parameters
	a = cat.time
	h = .6909
	G = 6.67e-11
#convert to m and physical not comoving coordinates
	particle_radius_m = particle_radius2 * (3.1e19) * a/h
#masses from 1e10 solar mass to kg			
	masses = cat['Masses'] * (1e10)*(1.99e30)/h
	masses = masses[nonzero]
	m_enclosed = []
	for r in particle_radius_m:
		inner_radii = particle_radius_m < r
		m_enc = masses[inner_radii].sum()
		m_enclosed.append(m_enc)
	df = pd.DataFrame()
	df['radius'] = particle_radius_m
	df['m_enclosed'] = np.array(m_enclosed)
	df = df.sort_values(by='radius', ascending=False)
	radii_sorted = df['radius'].array
	mass_sorted = df['m_enclosed'].array
	v_0 = np.sqrt(2*G*mass_sorted[0]/radii_sorted[0])
	print(radii_sorted[0:5])
	print('///////Results From Calculating Directly////////')
	print('V_0 = ' + str(v_0))
	M_r = mass_sorted/radii_sorted
	U = -1*G*(mass_sorted/radii_sorted)
	dU = (U[1:] - U[:-1])
	delta_u = (np.sum(dU))
	print('Delta U = ' + str(delta_u))
	v_f = np.sqrt((2*delta_u) + (np.square(v_0)))
	print('V_f = ' + str(v_f))
	m_out = 10e42/(.5*(np.square(v_f))) / (1.99e30)
	print('M_out/M_formed = ' + str(m_out))
	print('////////Results Assuming NFW/////////////')
#try assuming an NFW profile
	model = NFWProfile()
#need max radius in Mpc/h and max mass in M_sun/h
	max_radius = particle_radius.max()/1000 *a/h
	max_mass = cat['Masses'].sum() * 1e10
	radius = np.logspace(-4, np.log10(max_radius), 1000)
	total_mass = np.zeros(1000) + max_mass
	conc = 11/8
	mass_enc = model.enclosed_mass(radius, total_mass, conc)
	df2 = pd.DataFrame()
	df2['radius'] = radius
	df2['mass_enc'] = mass_enc
	df2 = df2.sort_values(by='radius', ascending=False)
	radii = df2['radius'].array
	mass = df2['mass_enc'].array
#unit conversions because result is in M_sun/h
	mass = mass*(1.99e30)/.702
	radii = radii*3.1e22*a/.702
	v_0 = np.sqrt(2*G*mass[0]/radii[0])
	print(max_radius*3.1e22)
	print(particle_radius.max()*a/h)
	print(radii[-1])
	print('V_0 = ' + str(v_0))
	U = -1*G*(mass/radii)
	dU = (U[1:] - U[:-1])
	delta_u = np.sum(dU)
	print('Delta U = '+ str(delta_u))
	v_f = np.sqrt((2*delta_u) + (np.square(v_0)))
	print('V_f = ' + str(v_f))
	m_out = 10e42/(.5*(np.square(v_f))) / (1.99e30)
	print('M_out/M_formed = ' + str(m_out))
	return


if __name__=='__main__':
	main()

