
import numpy as np
import matplotlib.pyplot as plt
from readData.DataLoader import DataLoader
from matplotlib.colors import LogNorm
from astropy.cosmology import Planck13
import pandas as pd
def main():
        snap = 120
        cat = DataLoader('../m10_destructive_run/mb_mass10-output', snap, 1, ['GroupPos','SubhaloMass','SubhaloPos', 'GroupMassType'])
        print(cat['GroupPos'][0:5])
        print(cat['GroupMassType'][0:5])
        print(cat.time)
        return


if __name__=='__main__':
	main()

