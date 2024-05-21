import numpy as np
import matplotlib.pyplot as plt
from readData.DataLoader import DataLoader
from matplotlib.colors import LogNorm
import pandas as pd
import shared_data
from  astropy.cosmology import FlatLambdaCDM
df = pd.read_csv('../ExpansionList_1', sep=' ', header=None)
scale_factors = df[2].array[45:]
z = (1/scale_factors - 1)
cosmo = FlatLambdaCDM(H0=69.09, Om0=.3, Tcmb0=2.725)
t = cosmo.lookback_time(z)
last_Gyr = t.value<=1
for x in range(9):
    print(cosmo.lookback_time(x).value)
def alpha(filepath, upper, lower, spacing):
    alpha_list =[]
    for i in range(121, 128):
        cat = DataLoader(filepath, i, 1, ['Masses','Coordinates','SubhaloPos'], sub_idx=0)
        particle_radius = np.sqrt(np.sum((np.square(cat['Coordinates']-cat['SubhaloPos'])), axis=1))
        radius = np.logspace(-.973,1,spacing)
        plot_radius = 1/2*(radius[0:-1] + radius[1:])
        shell_volume = (4/3) * np.pi * (radius[1:]**3 - radius[0:-1]**3)
        masses = cat['Masses'] * (1e10)
        hist, bin_edges = np.histogram(particle_radius, radius, weights=masses)
        density = hist/shell_volume
        logp = np.log10(density)
        logr = np.log10(radius)
        dlogr = logr[upper] - logr[lower]
        dlogp = logp[upper] - logp[lower]
        alpha = dlogp/dlogr
        alpha_list.append(alpha)
    final = np.array(alpha_list)
    final_med = np.median(final)
    print(radius[upper])
    print(radius[lower])
    return final_med
alpha('../m10_ic0_sg/sg_mass10-output', upper=1, lower=0, spacing=20)



