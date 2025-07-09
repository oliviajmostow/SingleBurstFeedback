from scipy.optimize import curve_fit
import pandas as pd
import numpy as np
import matplotlib, os
from matplotlib import rc
from matplotlib import pyplot as plt

def get_names(run_string):

    all_runs = run_string.strip().split(",")

    names = []
    for run in all_runs:
       
        name = ''

        if 'MW' in run:
            name += run.split("_")[0]
            name += " "

        if 'CDM' in run:
            name += 'CDM'
        else:
            if 'DM' in run:
                name += f'ETHOS$_{run[-4]}$'
            else:
                name += f'ETHOS$_{run.split("/")[0][-1]}$'

        if 'DM' in run and 'CDM' not in run:
            name += " DMO"

        if 'CDMDM' in run:
            name += " DMO"

        names.append(name)

    return names, all_runs

def get_snap(snap_str):
    
    snap_li_str = snap_str.strip().split(",")
    snap_li = [int(el) for el  in snap_li_str]

    if len(snap_li) == 1:
        snap_li = snap_li[0]

    return snap_li

def set_ticks(ax):
    ax.tick_params('both', which='minor', length=4, direction='in', bottom=True, top=True, left=True, right=True)
    ax.tick_params('both', which='major', length=8, direction='in', bottom=True, top=True, left=True, right=True)

    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(2)

    return 

def set_plot_params(nrows=1, ncols=1, figsize=None):

    plt.rc('font',**{'family':'STIXGeneral'})
    plt.rc('text', usetex=True)

    plt.rc('font', size=12)          # controls default text sizes
    plt.rc('axes', titlesize=16)     # fontsize of the axes title
    plt.rc('axes', labelsize=24)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=18)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=18)    # fontsize of the tick labels
    plt.rc('legend', fontsize=16)

    plt.rc('lines', linewidth=2)

    if figsize is None:
        fig, ax = plt.subplots(nrows,ncols, figsize=(ncols*5 + (ncols)*3,nrows*5+(nrows-1)*3))
    else:
         fig, ax = plt.subplots(nrows,ncols, figsize=figsize)

    if type(ax)==type(np.zeros(1)):
        for a in ax.ravel():
            set_ticks(a)
    else:
        set_ticks(ax)

    return fig, ax

def change_plot_params():
    plt.rc('font',**{'family':'STIXGeneral'})
    plt.rc('text', usetex=True)

    plt.rc('font', size=12)          # controls default text sizes
    plt.rc('axes', titlesize=16)     # fontsize of the axes title
    plt.rc('axes', labelsize=24)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=18)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=18)    # fontsize of the tick labels
    plt.rc('legend', fontsize=16)

    plt.rc('lines', linewidth=3)

    return None

def get_fof_idx(name, snap):

    if not os.path.exists(f"plots/tracking/{name}.txt"):
        print(f"This galaxy has not been tracked: {name}")
        return None

    all_idx = np.loadtxt(f"plots/tracking/{name}.txt")
    idx = 127 - snap

    if idx >= len(all_idx):
        return None

    return int(all_idx[idx])

def get_sub_idx(name, snap, sub, start_snap=127):

    if start_snap == 127:
        path = f"plots/tracking/{name}/sub_{sub}.txt"
    else:
        path = f"plots/tracking/{name}_{start_snap}/sub_{sub}.txt"

    if not os.path.exists(path):
        print(f"This galaxy has not been tracked: {name} {sub}")
        return None, None

    all_idx = np.loadtxt(path)
    idx = start_snap - snap

    if idx >= len(all_idx):
        return None, None

    if all_idx[idx][0] < 0:
        return None, None

    return int(all_idx[idx][0]), int(all_idx[idx][1]) #fof, sub

def get_matched_sub(name, sub):

    if name == 'CDM':
        return sub

    all_idx = np.loadtxt('plots/tracking/sub_list.txt')
    idx = int(name.split("_")[1][0])
    
    scut = all_idx[:,0] == sub

    if np.sum(scut) == 0:
        return -1

    return int(all_idx[scut, idx])
    

class ObsData():

    def __init__(self, path, points=False, lines=False, error=False, y_is_log=True, xerror=False, yerror=False):
        self.path = path
        self.y_is_log = y_is_log
        self.error = error
        self.xerror = xerror
        self.yerror = yerror

        self.cols = None
        self.data = self.read_data()
        
        if xerror is not False or yerror is not False:
            self.error = xerror if xerror is not False else yerror 

        if lines:
            self.lines = self.calc_lines()
        if points:
            self.xpoints, self.ypoints = self.get_points()
        if self.error is not False:
            self.xpoints, self.ypoints, self.xerr, self.yerr = self.calc_error()

    def calc_error(self):
        xpoints = []
        ypoints = []
        xerr = []
        yerr = []
        for col in self.data:
            x = np.array(self.data[col]['x'])
            y = np.array(self.data[col]['y'])

            ymax = np.max(y)
            ymin = np.min(y)
            xmax = np.max(x)
            xmin = np.min(x)

            xpoint = np.median(x[(x>xmin) & (x<xmax)])
            ypoint = np.median(y[(y>ymin) & (y<ymax)])

            if len(x)<5:
                if self.xerror is not False:
                    ypoint = np.median(y)
                elif self.yerror is not False:
                    xpoint = np.median(x)
                else:
                    bounds = np.array([xmax, xmin, ymax, ymin])
                    cut = (bounds-xpoint == 0) | (bounds-ypoint == 0)
                    if cut[0]:
                        xmax = xpoint + (xpoint - xmin)
                    if cut[1]:
                        xmin = xpoint - (xmax - xpoint)
                    if cut[2]:
                        ymax = ypoint + (ypoint - ymin)
                    if cut[1]:
                        ymin = ypoint - (ymax - ypoint)                   

            xpoints.append(xpoint)
            ypoints.append(ypoint)
            if self.error == 'same':
                xerr.append((xmax-xpoint) + (xpoint-xmin) /2)
                yerr.append((ymax-ypoint) + (ypoint-ymin) /2)
            elif self.error == 'unique':
                xerr.append([xmax-xpoint, xpoint-xmin])
                yerr.append([ymax-ypoint, ypoint-ymin])
            else:
                raise KeyError(f"{self._error} not supported for error type, please use 'unique' or 'same")
        return xpoints, ypoints, xerr, yerr

    def get_points(self):
        xpoints = [self.data[key]['x'] for key in self.data]
        ypoints = [self.data[key]['y'] for key in self.data]
        return xpoints, ypoints

    def lin(self, x, m, b):
        return m*x + b

    def calc_lines(self):
        lines = {el:[] for el in self.cols}
        for col in self.cols:
            x = self.data[col]['x']
            y = self.data[col]['y']
            guess = [(y.iloc[-1]-y.iloc[0]) / (x.iloc[-1]-x.iloc[0]), 1]
            if self.y_is_log:
                y = np.log10(y)
            popt, pcov = curve_fit(self.lin, x, y, p0=guess) 
            lines[col] = popt
        return lines

    def read_data(self):
      
        df = pd.read_csv(self.path)
        self.cols = [el for el in set(df.columns) if 'named' not in el]

        x = None
        y = None
        col_name = ''
        data = dict()
        for col in df.columns:
            if x is None:
                x = df.loc[1:,col].dropna().astype(float)
                col_name = col
            elif y is None:
                y = df.loc[1:,col].dropna().astype(float)
                data[col_name] = {'x':x, 'y':y}
                x = None
                y = None
                col_name = ''

        return data


                
