
import numpy as np
import matplotlib.pyplot as plt
from readData.DataLoader import DataLoader
from matplotlib.colors import LogNorm
import pandas as pd
def main():
	cat = DataLoader('../hr_sg_mass9/hr_sg_mass9-output', 127, 1, ['Masses','Coordinates','SubhaloPos','SubhaloMass'], sub_idx=0)
	#get radius from center and coordinates
	sub = cat['Coordinates'] - cat['SubhaloPos']
	sub2 = np.square(sub)
	sub_sum = np.sum(sub2, axis=1)
	particle_radius = np.sqrt(sub_sum)
	log_part_radius = np.log(particle_radius)
	#evenly spaced values in lin space to put the particles into bins based on radius
#	logr = np.linspace(1e-1, log_part_radius.max(),125)
#	plot_radius = 1/2*(logr[0:-1] + logr[1:]) #midpoints of all the bins
	#shell volume of space between adjacent values for bins to use in computing density
#	volume = (4/3) * np.pi * (10**logr)**3
#	volume2 =(4/3) * np.pi * ((10**logr)[1:])**3
#	shell_volume  = volume2 -  volume[0:-1]
	#convert to solar masses
	masses = cat['Masses'] * (1e10)
	#compute density from weighted sum inside each bin/volume of shell
#	hist, bin_edges = np.histogram(log_part_radius, logr, weights=masses)
#	density = hist/shell_volume
#	logp = np.log(density)
#	plot_logr = plot_radius #right now these are the same length
	#calculating derivative of density w.r.t. radius, using plot radius to have same shape
#	delta_logp = logp[1:] - logp[0:-1]
#	delta_logr = logr[1:] - logr[0:-1]
#	plot_logr_midpoints = 1/2*(plot_logr[1:]+plot_logr[0:-1])
#	delta_logr = plot_logr_midpoints[1:] - plot_logr_midpoints[:-1]
#	dp_dr = delta_logp/delta_logr
#	fig, ax = plt.subplots()
#	ax.plot((10**plot_logr_midpoints), dp_dr)
#	fig.savefig('log_firstorderslope.pdf')
#	

#	density calculation from other file
	radius = np.logspace(-1,1,50)
	logr = np.log10(radius)
	radius_normal = np.logspace(-1,1,50)
	radius_shifted = radius_normal[1:]
	plot_radius = 1/2* (radius_normal[0:-1] + radius_shifted)
	radius2 = radius[1:]
	volume = (4/3) * np.pi * (radius)**3
	volume2 = (4/3) * np.pi * (radius2)**3
	shell_volume  = volume2 -  volume[0:-1]
	hist, bin_edges = np.histogram(particle_radius, radius, weights=masses)
	density = hist/shell_volume
	fig, ax = plt.subplots()
	ax.plot(plot_radius, density)
	ax.set_xscale('log')
	ax.set_yscale('log')
	logp = np.log10(density)
	print(density[:5])
	print(logp[:5])
	fig.savefig('testplot.pdf')
	fig, ax = plt.subplots()
	plot_logr = np.log10(plot_radius)
	ax.plot(plot_logr, logp)
	fig.savefig('testplot2.pdf')
	cat2 = DataLoader('../hr_sg_mass10/hr_sg_mass10-output', 127, 1, ['Masses','Coordinates','SubhaloPos','SubhaloMass'], sub_idx=0)
	sub = cat2['Coordinates'] - cat['SubhaloPos']
	sub2 = np.square(sub)
	sub_sum = np.sum(sub2, axis=1)
	particle_radius = np.sqrt(sub_sum)
	log_part_radius10 = np.log(particle_radius)
	masses10 = cat2['Masses']*1e10
	radius10 = np.logspace(-1,1,50)
	logr10 = np.log10(radius10)
	plot_radius10 = 1/2*(radius10[:-1] + radius10[1:])
	shell_volume = 4/3 * np.pi * ((radius10[:-1]**3) - (radius10[1:]**3))
	hist, bin_edges = np.histogram(particle_radius, radius10, weights=masses10)
	density10 = hist/shell_volume
	logp10 = np.log10(density10)
	logr10 = np.log10(radius10)
	plot_logr10 = np.log10(plot_radius10)
	
	#second order approximation (formula from numpy calculated manually)
	f_xminus = logp[0:-2] #starts one before f_x, ends one before
	f_x = logp[1:-1] #excludes first
	f_xplus = logp[2:] #starts one after, also needs to end one after to be same length
	plr2 = logr[1:-1]#the midpoints of the bins that are included in the approximation, first and last are cut off because can't go back/forth one
	h_s = plot_logr[2:] - plot_logr[1:-1]
	h_d = plot_logr[1:-1] - plot_logr[0:-2]
	dp_dr = ((np.square(h_s))*(f_xplus) + (np.square(h_d) - np.square(h_s))*f_x - (np.square(h_d))*f_xminus)/(h_s * h_d * (h_d + h_s))
#	fig, ax = plt.subplots()
	#dp_dr has length two less than original, which means it should be plotted against plr2
	plr2_midpoints = 1/2* (plr2[1:]+plr2[:-1])
#	ax.plot((10**(plr2_midpoints)), dp_dr)
#	ax.set_xscale('log')
#	ax.set_xlabel('r')
#	ax.set_ylabel('dlog(p)/dlog(r)')
#	ax.set_xlim(0,500)
#	fig.savefig('slope2_sg9.pdf')
	

	new_bins = plr2[::5]
	hist, bin_edges = np.histogram(plr2_midpoints, new_bins, weights = dp_dr)
	hist2, bin_edges = np.histogram(plr2_midpoints, new_bins)
	average_slopes = hist/hist2
	fig, ax = plt.subplots()
	ax.plot(10**(1/2*(new_bins[1:]+new_bins[:-1])), average_slopes)
	ax.set_xlabel('Radius (kpc)')
	ax.set_ylabel('\u03B1')
	ax.set_xscale('log')
	ax.set_ylim(-3, 1)
#	fig.savefig('slope2_sg9_smooth.pdf')
#	same thing but for the 10 solar mass data
	f_xminus = logp10[0:-2]
	f_x = logp10[1:-1]
	plr2_10 = logr10[1:-1]
	h_s = plot_logr10[2:]-plot_logr[1:-1]
	h_d = plot_logr[1:-1] - plot_logr[0:-2]
	dp_dr_10 = ((np.square(h_s))*(f_xplus) + (np.square(h_d) - np.square(h_s))*f_x - (np.square(h_d))*f_xminus)/(h_s * h_d * (h_d + h_s))
	plr2_midpoints_10 = 1/2*(plr2_10[1:] + plr2_10[:-1])
	new_bins = plr2_10[::5]
	hist, bin_edges = np.histogram(plr2_midpoints_10, new_bins, weights = dp_dr_10)
	hist2, bin_edges = np.histogram(plr2_midpoints_10, new_bins)
	average_slopes_10 = hist/hist2
	ax.plot(10**(1/2*(new_bins[1:]+new_bins[:-1])), average_slopes_10)
	fig.savefig('slope_sg.pdf')
	
#	fourth order derivative approximation (using even spacing formula)
	h = logr.max()/50
	pr4 = logr[2:-2]
	f_x = logp[2:-2]
	f_xplus2 = logp[4:]
	f_xminus2 = logp[:-4]
	f_xminus = logp[1:-3]
	f_xplus = logp[3:-1]
	dp_dr = (-f_xplus2 + 8*f_xplus -8*f_xminus + f_xminus2)/(12*h)
#	dp_dr has length 4 less than original, plot against pr4
	fig, ax = plt.subplots()
	ax.plot(10**(1/2*(pr4[1:]+pr4[:-1])), dp_dr)
	ax.set_xlabel('r')
	ax.set_ylabel('dp/dr')
	ax.set_xscale('log')
	fig.savefig('slope_sg10.pdf', bbox_inches='tight')
	new_bins = pr4[::4]
	pr4_midpoints = 1/2* (pr4[1:]+pr4[:-1])
	hist, bin_edges = np.histogram(pr4_midpoints, new_bins, weights=dp_dr)
	hist2, bin_edges = np.histogram(pr4_midpoints, new_bins)
	average_slopes = hist/hist2
	fig, ax = plt.subplots()
	ax.plot(10**(1/2*(new_bins[1:]+new_bins[:-1])), average_slopes)
	ax.set_xlabel('Radius (kpc)')
	ax.set_ylabel('\u03B1')
#	ax.set_xlim(0, 200)
	ax.set_xscale('log')
	ax.set_ylim(-3, 1)
	for item in ([ax.title, ax.xaxis.label, ax.yaxis.label]): item.set_fontsize(18)
	fig.savefig('slope_sg10_smooth.pdf', bbox_inches='tight')
	return


if __name__=='__main__':
	main()

