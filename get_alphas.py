import numpy as np
import matplotlib.pyplot as plt
from readData.DataLoader import DataLoader
from matplotlib.colors import LogNorm
import pandas as pd
import calc_alpha
from  astropy.cosmology import FlatLambdaCDM
plt.rcParams.update({'font.size': 14})
def main():
        runs =['/blue/paul.torrey/oliviamostow/Projects/mag1_red8_g2_up/sb_mass8-output','../mag1_red7_g2_up/sb_mass8-output','../mag1_red6_g2_up/sb_mass8-output','../mag1_red5_g2_up/sb_mass8-output','../mag1_red4_g2_up/sb_mass8-output', '../mag1_red3_g2_up/sb_mass8-output','../mag1_red2_g2_up/sb_mass8-output']
        mag2runs =['../mag2_red8_g2_up/sb_mass8-output','../mag2_red7_g2_up/sb_mass8-output','../mag2_red6_g2_up/sb_mass8-output','../mag2_red5_g2_up/sb_mass8-output','../mag2_red4_g2_up/sb_mass8-output', '../mag2_red3_g2_up/sb_mass8-output','../mag2_red2_g2_up/sb_mass8-output']
        sizetestruns = ['../sizetest_9/mb_mass10-output', '../sizetest_8/mb_mass10-output','../sizetest_7/mb_mass10-output','../m10_ic0_sg/sg_mass10-output']
        fqruns8 = ['../ft_s8_20/mb_mass10-output', '../ft_s8_10/mb_mass10-output','../ft_s8_5/mb_mass10-output','../ft_s8_1/sb_mass10-output']
        fqruns9 = ['../ft_s9_20/mb_mass10-output', '../ft_s9_10/mb_mass10-output','../ft_s9_5/mb_mass10-output','../ft_s9_1/sb_mass10-output']
        mbm8runs = ['../mb_m1_r8_g2/mb_mass8-output','../mb_m1_r7_g2/mb_mass8-output','../mb_m1_r6_g2/mb_mass8-output','../mb_m1_r5_g2/mb_mass8-output','../mb_m1_r4_g2/mb_mass8-output','../mb_m1_r3_g2/mb_mass8-output','../mb_m1_r2_g2/mb_mass8-output','../mb_m1_r1_g2/mb_mass8-output']
#        print(calc_alpha.alpha('../m10_mb_highfreq_90_soft1/mb_mass10-output',upper=5, lower=3, spacing=20))
        fqruns7 = ['../ft_s7_20/mb_mass10-output', '../ft_s7_10/mb_mass10-output','../ft_s7_5/mb_mass10-output','../ft_s7_1/sb_mass10-output']
        alpha_list = []
        for i in range(len(sizetestruns)):
            alpha_list.append(calc_alpha.alpha(sizetestruns[i], upper=5, lower=3, spacing=20))
        print(alpha_list)
        return


if __name__=='__main__':
	main()

