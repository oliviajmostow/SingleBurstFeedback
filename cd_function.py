
import numpy as np
from readData.DataLoader import DataLoader
from scipy import stats

def calc_density(filepath, spacing, snap=127):
        cat = DataLoader(filepath, snap, 1, ['Masses','Coordinates','GroupPos','SubhaloMass','SubhaloPos', 'GroupMassType'], sub_idx=0)
        particle_radius = np.sqrt(np.sum(np.square(cat['Coordinates']-cat['SubhaloPos']), axis=1))
        radius = np.logspace(-1.2,1, spacing)
        plot_radius = 1/2*(radius[:-1] +radius[1:])
        shell_volume = (4/3) * np.pi * ((radius[1:])**3-radius[:-1]**3)
        masses = cat['Masses'] * (1e10)
        hist, bin_edges = np.histogram(particle_radius, radius, weights=masses)
        density = hist/shell_volume
        return density

def calc_vd(filepath, spacing, snap=127):
    cat = DataLoader(filepath, snap, 1, ['Velocities','Masses','Coordinates','GroupPos','SubhaloMass','SubhaloPos', 'GroupMassType'], sub_idx=0)
    velocities = cat['Velocities']*np.sqrt(cat.time)
    particle_radius = np.sqrt(np.sum((np.square(cat['Coordinates']-cat['SubhaloPos'])), axis=1))*cat.time
    std_x, b, bn = stats.binned_statistic(particle_radius, velocities[:,0],'std', bins=num_bins)
    std_x, b, bn = stats.binned_statistic(particle_radius, velocities[:,0],'std', bins=num_bins)
    std_z, b, bn = stats.binned_statistic(particle_radius, velocities[:,2],'std', bins=num_bins)
    dispersion = np.sqrt(np.square(std_x) + np.square(std_y) + np.square(std_z))
    plot_radius = 1/2*(b[0:-1] + b[1:])
    return density, plot_radius


