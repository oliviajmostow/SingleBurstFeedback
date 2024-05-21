from readData.DataLoader import DataLoader
import numpy as np
import matplotlib.pyplot as plt
import sys
import shared_data

from matplotlib.colors import LogNorm
from util.calc_hsml import get_gas_density_around_stars as get_dens

def cart(r, t, p):
    x = r*np.cos(t)*np.sin(p)
    y = r*np.sin(t)*np.sin(p)
    z = r*np.cos(p)
    return list(zip(x,y,z))

def main():

    path = "/home/j.rose/Projects/my_MW_DAO/IllustrisTNGcode/new_MW_2023/Hydro/"
    part_types = [1]
    keys = ["Coordinates", "Masses", "GroupPos", "Group_R_Crit200", "GroupMassType"] 

    names, all_runs = shared_data.get_names(sys.argv[1])
    snap_num = shared_data.get_snap(sys.argv[2])
    fig, axi = shared_data.set_plot_params(ncols=3, nrows=2)
    #fig, ax = shared_data.set_plot_params(ncols=10, nrows=5)
    axi = axi.ravel()
    
    ax = [axi[0]]*3 + [axi[1],axi[2]] + [axi[3]]*3 + [axi[4]]
    #ax = axi

    axi[-1].axis("off")

    #call = plt.cm.inferno(np.linspace(0,1,7))
    #colors = [call[0]] + list(call[2:])

    rmax = 200
    rmin = .1
    npoints = 30 #cubed
    soft = 0.305/2  #.305 == medium res

    for i,run in enumerate(all_runs):

        mass_ratio = []
        slopes = []

        cat = DataLoader(path+run, snap_num=snap_num, part_types=part_types, keys=['GroupPos', 'GroupMass', 'Group_R_Crit200', 'GroupMassType'])

        gal_pos = cat['GroupPos'][0]
        gal_r = np.sqrt(np.sum(np.square(cat['GroupPos'] - gal_pos), axis=1))

        cut = (cat['GroupMassType'][:,1]*1e10 > 2e8) & (cat['GroupMassType'][:,2] == 0) & (cat['GroupMassType'][:,4] > 0) & (cat['Group_R_Crit200']/cat.h > 30) # 30)

        #cut[0] = False
        #cut[20:] = False

        print(np.sum(cut))
        #continue

        for j,fof_idx in enumerate(np.arange(len(gal_r))[cut]):
            #break

            cat = DataLoader(path+run, snap_num=snap_num, part_types=part_types, keys=keys, fof_idx=fof_idx)

            r200 = cat['Group_R_Crit200'] / cat.h
            rmin = r200 * 0.01
            rmax = r200 * 0.02

            coords = (cat['Coordinates'] - cat['GroupPos']) / cat.h #* cat.time
            masses = cat['Masses'] *1e10/cat.h

            if len(masses) < 5:
                continue

            #create the log-spaced spherical grid of points to calculate the density around
            #all_r = np.logspace(np.log10(rmin), np.log10(rmax), npoints) #radial stepping
            all_r = np.logspace(-1,1,100)
            t = np.linspace(-np.pi, np.pi, npoints+1)[:-1] #theta span
            p = np.linspace(0, np.pi, npoints+1)[:-1] #phi span 
            p, t = np.meshgrid(p, t) #theta-phi grid
            points = np.array([cart(r,t, p) for r in all_r]) #converted to catesean coords

            #reorganize points
            x = points[:,:,0].ravel()
            y = points[:,:,1].ravel()
            z = points[:,:,2].ravel()

            #calculate dm densities around set points
            h = get_dens(coords[:,0], coords[:,1], coords[:,2], masses, x, y, z)

            #take averages and plot 
            np2 = npoints*npoints
            med = np.array([np.median(h[np2*i:np2*(i+1)]) for i in range(100)])

            mass_ratio.append(np.log10(cat['GroupMassType'][4] / cat['GroupMassType'][1]))

            rmin_idx = np.argmin(np.abs(all_r - rmin))
            rmax_idx = np.argmin(np.abs(all_r - rmax))
            slopes.append((np.log10(np.average(med[rmax_idx-1:rmax_idx+1])) - np.log10(np.average(med[rmin_idx-1:rmin_idx+1]))) / (np.log10(rmax) - np.log10(rmin)))

            #ax[j].plot(all_r, med, label=f'{mass_ratio[-1]:.2f}')

            #ax[j].set_yscale("log")
            #ax[j].set_xscale("log")

            #ax[j].plot([rmin, rmin],[np.min(med), np.max(med)], 'k--', alpha=0.5)
            #ax[j].plot([rmax, rmax],[np.min(med), np.max(med)], 'k--', alpha=0.5)

            #ax[j].set_title(f'{slopes[-1]:.2f}')
            #ax[j].legend()

        
        print(names[i], slopes)

        ax[i].plot(mass_ratio, slopes, '.', label=names[i], markersize=10) #, color=colors[i])
        ax[i].legend()

        ax[i].set_xlim((-5,-.5))
        ax[i].set_ylim((-3,1))

        ax[i].set_xlabel("log$_{10}$ (M$_*$/M$_{halo}$)")
        ax[i].set_ylabel("Average Field Slope")


    if True:
        obsd = shared_data.ObsData("obsdata/core_cusp_fire.csv")
        cols = [['FIRE_high', 'FIRE_low'], ['NFW_high', 'NFW_low']]
        colors = ['b', 'k']
        for k,col in enumerate(cols):
            x_high = obsd.data[col[0]]['x']
            y_high = obsd.data[col[0]]['y']
            p_high = np.poly1d(np.polyfit(x_high, y_high, 10))
            x_low = obsd.data[col[1]]['x']
            y_low = obsd.data[col[1]]['y']
            p_low = np.poly1d(np.polyfit(x_low, y_low, 10))
            xplot = np.linspace(-6,0, 100)

            for ax in axi[:-1]:
                if k == 0:
                    ax.fill_between(xplot, p_high(xplot), p_low(xplot), color=colors[k], alpha=0.3)
                if k == 1:
                    ax.fill_between(xplot, np.average(y_high), np.average(y_low), color=colors[k], alpha=0.3)

    fig.savefig("plots/cores_field.pdf", bbox_inches='tight')

    return

if __name__=="__main__":
    main()
