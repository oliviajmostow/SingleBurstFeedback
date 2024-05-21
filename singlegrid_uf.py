
import numpy as np
import matplotlib.pyplot as plt
from readData.DataLoader import DataLoader
import matplotlib as mpl
import pandas as pd
import cd_function
import shared_data
import matplotlib as mpl
def main():
    a = [-1.41, -2.13, -1.77, -1.71, -1.71, -1.69, -1.39, -2.01, -1.87, -1.57,-1.57,-1.53, -1.46, -1.78,-1.83]
    fig, ax = shared_data.set_plot_params()
    alphas = np.empty([2,7])
    alphas[0] = a[:7]
    alphas[1] = a[7:14]
    #change extent so gridlines surround each square
    extent = (0, alphas.shape[1], alphas.shape[0],0)
    ax.set_yticks([0,1],labels=['.2','1'][::-1])
    ax.set_xticks([0,1,2,3,4,5,6],labels=['13.0','12.9','12.7','12.5','12.1','11.5','10.4'])
    norm = mpl.colors.Normalize(vmin=-2.0, vmax=-.2)
    cmap = mpl.colors.LinearSegmentedColormap.from_list('my_colormap',['blue','red'],256)
    img = ax.imshow(alphas, interpolation='nearest',cmap=cmap, norm=norm)
    #workaround to get gridlines
    plt.hlines(y=np.arange(0,2)+.5, xmin=np.full(2,0)-.5, xmax=np.full(2,7)-.5, color='w')
    plt.vlines(x=np.arange(0,7)+.5, ymin=np.full(7,0)-.5, ymax=np.full(7,2)-.5, color='w')
    for (j,i), label in np.ndenumerate(alphas):
        ax.text(i,j, label, ha='center', va='center',color='w', size=20)
    fig.colorbar(img, ax=ax, norm=norm).set_label(r'$\alpha$')
    ax.set_xlabel('Time of Burst (Gyr)')
    ax.set_ylabel(r'Mass Expelled ($10^{5}M_{\odot}/h$)')
    #ax.grid(True, color='black')
    fig.savefig('cmap.pdf',bbox_inches='tight')
    return


if __name__=='__main__':
	main()
 
