
import numpy as np
import matplotlib.pyplot as plt
from readData.DataLoader import DataLoader
from matplotlib.colors import LogNorm
from astropy.cosmology import Planck13
import pandas as pd
def main():
        snap = 78
        cat = DataLoader('../mag1_red2_g2_up/sb_mass8-output', snap, 1, ['Masses'])
        print(cat.time)
        return


if __name__=='__main__':
	main()

