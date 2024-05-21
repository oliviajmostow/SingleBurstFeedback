
import numpy as np
import matplotlib.pyplot as plt
from readData.DataLoader import DataLoader
from matplotlib.colors import LogNorm
import pandas as pd
import shared_data
def main():
        spacing = 15
        latedensities = []
        earlydensities = []
        middensities = []
        runs = ['../tt1_m1_r8/sb_mass8-output','../tt1_m1_r7/sb_mass8-output','../tt1_m1_r6/sb_mass8-output','../tt1_m1_r5/sb_mass8-output','../tt1_m1_r4/sb_mass8-output','../tt1_m1_r3/sb_mass8-output','../tt1_m1_r2/sb_mass8-output','../tt1_m1_r1/sb_mass8-output']
        for i in range(len(runs)):
            cat = DataLoader(runs[i], 127, 1, ['Masses','Coordinates','GroupPos','SubhaloMass','SubhaloPos', 'GroupMassType'], sub_idx=0)
            sub = cat['Coordinates'] - cat['SubhaloPos']
            sub2 = np.square(sub)
            sub_sum = np.sum(sub2, axis=1)
            particle_radius = np.sqrt(sub_sum)
            radius = np.logspace(-1.35,1, spacing)
            radius_normal = np.logspace(-1.35, 1, spacing)
            radius_shifted = radius_normal[1:]
            plot_radius = 1/2*(radius_normal[0:-1] + radius_shifted)
            radius2 = radius[1:]
            volume = (4/3) * np.pi * (radius)**3
            volume2 =(4/3) * np.pi * (radius2)**3
            shell_volume  = volume2 -  volume[0:-1]
            volume_p1 = (4/3) * np.pi * (radius[0])**3
            vol_list =  list(shell_volume)
            shell_volume = np.array(vol_list)
            masses = cat['Masses'] * (1e10)
            hist, bin_edges = np.histogram(particle_radius, radius, weights=masses)
            density = hist/shell_volume
            earlydensities.append(density)
        runs = ['../tt2_m1_r8/sb_mass8-output','../tt2_m1_r7/sb_mass8-output','../tt2_m1_r6/sb_mass8-output','../tt2_m1_r5/sb_mass8-output','../tt2_m1_r4/sb_mass8-output','../tt2_m1_r3/sb_mass8-output','../tt2_m1_r2/sb_mass8-output','../tt2_m1_r1/sb_mass8-output']
        for i in range(len(runs)):
            cat = DataLoader(runs[i], 127, 1, ['Masses','Coordinates','GroupPos','SubhaloMass','SubhaloPos', 'GroupMassType'], sub_idx=0)
            sub = cat['Coordinates'] - cat['SubhaloPos']
            sub2 = np.square(sub)
            sub_sum = np.sum(sub2, axis=1)
            particle_radius = np.sqrt(sub_sum)
            radius = np.logspace(-1.35,1, spacing)
            radius_normal = np.logspace(-1.35, 1, spacing)
            radius_shifted = radius_normal[1:]
            plot_radius = 1/2*(radius_normal[0:-1] + radius_shifted)
            radius2 = radius[1:]
            volume = (4/3) * np.pi * (radius)**3
            volume2 =(4/3) * np.pi * (radius2)**3
            shell_volume  = volume2 -  volume[0:-1]
            masses = cat['Masses'] * (1e10)
            hist, bin_edges = np.histogram(particle_radius, radius, weights=masses)
            density = hist/shell_volume
            latedensities.append(density)
        runs = ['/blue/paul.torrey/oliviamostow/Projects/mag1_red8_g2_up/sb_mass8-output','../mag1_red7_g2_up/sb_mass8-output','../mag1_red6_g2_up/sb_mass8-output','../mag1_red5_g2_up/sb_mass8-output','../mag1_red4_g2_up/sb_mass8-output','../mag1_red3_g2_up/sb_mass8-output','../mag1_red2_g2_up/sb_mass8-output','../mag1_red1_g2_up/sb_mass8-output']
        for i in range(len(runs)):
            cat = DataLoader(runs[i], 127, 1, ['Masses','Coordinates','GroupPos','SubhaloMass','SubhaloPos', 'GroupMassType'], sub_idx=0)
            sub = cat['Coordinates'] - cat['SubhaloPos']
            sub2 = np.square(sub)
            sub_sum = np.sum(sub2, axis=1)
            particle_radius = np.sqrt(sub_sum)
            radius = np.logspace(-1.35,1, spacing)
            plot_radius = 1/2*(radius_normal[0:-1] + radius[1:])
            volume = (4/3) * np.pi * (radius)**3
            volume2 =(4/3) * np.pi * (radius[1:])**3
            shell_volume  = volume2 -  volume[0:-1]
            masses = cat['Masses'] * (1e10)
            hist, bin_edges = np.histogram(particle_radius, radius, weights=masses)
            density = hist/shell_volume
            middensities.append(density)
        print(earlydensities)
        return


if __name__=='__main__':
	main()

