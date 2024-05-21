
import numpy as np
import matplotlib.pyplot as plt
from readData.DataLoader import DataLoader
from matplotlib.colors import BoundaryNorm
import matplotlib as mpl
import pandas as pd
import cd_function
import shared_data
import matplotlib as mpl
def main():
    alphas = [-0.84,-0.66,-0.33,-0.44,-0.43,-1.09,-0.51,-0.57, -0.43,-0.32,-1.09,-1.07,-0.82,-0.94,-0.66, np.nan,np.nan,np.nan,np.nan, -1.05]
    fig, ax = shared_data.set_plot_params()
    m10_alphas = np.empty([4,5])
    m10_alphas[0] = alphas[:5]
    m10_alphas[1] = alphas[5:10]
    m10_alphas[2] = alphas[10:15]
    m10_alphas[3] = alphas[15:]
    #change extent so gridlines surround each square
    extent = (0, m10_alphas.shape[1], m10_alphas.shape[0],0)
    ax.set_yticks([0,1,2,3],labels=['.4','1','2.5','5'][::-1])
    ax.set_xticks([0,1,2,3,4],labels=['1','5','10','20','31'])
    norm = mpl.colors.Normalize(vmin=-1.25, vmax=-0.25)
    cmap = mpl.colors.LinearSegmentedColormap.from_list('my_colormap',['blue','red'],256)
    img = ax.imshow(m10_alphas, interpolation='nearest',cmap=cmap, norm=norm)
    #workaround to get gridlines
    plt.hlines(y=np.arange(0,3)+.5, xmin=np.full(3,0)-.5, xmax=np.full(3,5)-.5, color='w')
    plt.vlines(x=np.arange(0,4)+.5, ymin=np.full(4,0)-.5, ymax=np.full(4,4)-.5, color='w')
    for (j,i), label in np.ndenumerate(m10_alphas):
        ax.text(i,j, label, ha='center', va='center',color='w', size=20)
    #fig.colorbar(img, ax=ax, norm=norm).set_label(r'$\alpha$')
    cb = plt.colorbar(img, ax=ax, norm=norm)
    cb.set_label(r'$\alpha$')
    h, l = cb.ax.get_legend_handles_labels()
    cb.ax.axhline(y=-1.14, color='white',label='smooth growth')
    cb.ax.legend(bbox_to_anchor=[6,-0.05],facecolor='lightgray', framealpha=.5)
    #set coordinates of the labels
    size_x = m10_alphas.shape[0]
    size_y = m10_alphas.shape[1]
    x_end = m10_alphas.shape[1]
    x_start = 0
    y_end = 0
    y_start = m10_alphas.shape[1]
    jump_x = (x_end - x_start)/(2.0*size_x)
    jump_y = (y_end - y_start)/(2.0*size_y)
    x_positions = np.linspace(start=x_start, stop=x_end, num=size_x, endpoint=False)
    y_positions = np.linspace(start=y_start, stop=y_end, num=size_y, endpoint=False)
   # for y_index, y in enumerate(y_positions):
    #    for x_index, x in enumerate(x_positions):
     #       label = m10_alphas[x_index, y_index]
      #      text_x = x + jump_x
       #     text_y = y + jump_y
        #    ax.text(text_x, text_y, label, color='white', ha='center', va='center'))
    ax.set_xticks([0,1,2,3,4],labels=['1','5','10','20','31'])
    ax.set_yticks([0,1,2,3],labels=['.4','1','2.5','5'][::-1])
    ax.set_xlabel('Number of Bursts')
    ax.set_ylabel(r'Mass Expelled ($10^{7}M_{\odot}/h$)')
    #ax.grid(True, color='black')
    fig.savefig('cmap.pdf',bbox_inches='tight')
    return


if __name__=='__main__':
	main()
 
