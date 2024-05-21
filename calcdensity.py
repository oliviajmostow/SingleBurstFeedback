
import numpy as np
import matplotlib.pyplot as plt
from readData.DataLoader import DataLoader
from matplotlib.colors import LogNorm
import pandas as pd
import shared_data
plt.rcParams.update({'font.size': 14})
def main():
        spacing = 15
        colors = ['black', 'blue','red','orange','pink','green','purple','gray','lightblue']
        fig, ax = shared_data.set_plot_params()
        snaps= [127,127,127,127,127,127,127,127,127]
        labels=['size 1','steady','size 2','size 3', 'size 4', 'size 5', 'size 6','size 8','size 9']
        runs=['../m10_ic0_sg/sg_mass10-output','../m10_mb_highfreq_90/mb_mass10-output','../m10_mb_highfreq_90_soft2/mb_mass10-output','../m10_mb_highfreq_90_soft2/mb_mass10-output']
        runs=['../m10_destructive_run/mb_mass10-output','../m10_ic0_sg/sg_mass10-output','../sizetest_2/mb_mass10-output','../sizetest_3/mb_mass10-output','../sizetest_4/mb_mass10-output','../sizetest_5/mb_mass10-output','../sizetest_6/mb_mass10-output','../sizetest_8/mb_mass10-output','../sizetest_9/mb_mass10-output']
        runs=['../m10_ic0_sg/sg_mass10-output','../sizetest_8/mb_mass10-output','../sizetest_9/mb_mass10-output']
        labels=['steady','repeated bursts (larger mass expelled)','repeated bursts (smaller mass expelled)']
        colors=['orange','blue','gray','lightgray','black']
        runs = ['../m10_ic0_sg/sg_mass10-output','../sizetest_9/mb_mass10-output','../ft_s9_20/mb_mass10-output','../ft_s9_10/mb_mass10-output','../ft_s9_5/mb_mass10-output']
        labels = ['steady','30 bursts','20 bursts','10 bursts','5 bursts']
        runs = ['../sg_m8_g2/sg_mass8-output','../mag1_red7_g2_up/sb_mass8-output']
        labels = ['steady growth', 'single burst']
        colors = ['black','red']
        for i in [0]:
                cat = DataLoader(runs[i], snaps[i], 1, ['Masses','Coordinates','GroupPos','SubhaloMass','SubhaloPos', 'GroupMassType'], sub_idx=0)
                sub = cat['Coordinates'] - cat['SubhaloPos']
                sub2 = np.square(sub)
                sub_sum = np.sum(sub2, axis=1)
                particle_radius = np.sqrt(sub_sum)
                radius = np.logspace(-1.2,1, spacing)
                print(cat.time)
                radius_normal = np.logspace(-1.2, 1, spacing)
                radius_shifted = radius_normal[1:]
                plot_radius = 1/2*(radius_normal[0:-1] + radius_shifted)
                radius2 = radius[1:]
                volume = (4/3) * np.pi * (radius)**3
                volume2 =(4/3) * np.pi * (radius2)**3
                shell_volume  = volume2 -  volume[0:-1]
                volume_p1 = (4/3) * np.pi * (radius[0])**3
                masses = cat['Masses'] * (1e10)
                hist, bin_edges = np.histogram(particle_radius, radius, weights=masses)
                density = hist/shell_volume
#		fig, ax = shared_data.set_plot_params()
                ax.plot(plot_radius, density, color=colors[i], label=labels[i])
        ax.set_ylim((1000, density.max()*5))
        ax.set_xscale('log')
        ax.set_yscale('log')
#        ax.axvline(x=.076, linestyle='--', color='black', label='softening length')
        ax.axvline(x=2.8*.076, linestyle='--', color='gray', label='2.8 x softening')
        plt.xlabel('Radius (kpc)')
        ax.legend()
        plt.ylabel('Density ($M_{\odot}$/kpc$^3$)')
        fig.savefig('density.pdf', bbox_inches='tight')
        return


if __name__=='__main__':
	main()

