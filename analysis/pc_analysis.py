####  STILL PRELIMINARY, NEEDS MORE WORK :: adapted from Mathematica notebook ####

## run after aggregating data ##

import scipy as sp
import numpy as np
import pandas as pd
import sys
import os
import re
import matplotlib.pyplot as plt
from cycler import cycler

from scipy.optimize import curve_fit


def read_input(pc_file_path):
    num_pcs = len([name for name in os.listdir(pc_file_path) if (os.path.isfile(pc_file_path + name) and "permutation" in name)])

    pc_dfs = []
    for i in range(num_pcs):
        pc_dfs.append(pd.read_csv(pc_file_path+'permutation_data_{}.csv'.format(i+1),header=None))
    temp_df = pd.concat(pc_dfs)
    temperatures = temp_df[0].unique()

    pc_df = pd.DataFrame(columns=temp_df.columns)
    pc_df_std = pd.DataFrame(columns=temp_df.columns)
    for temp in temperatures:
        pc_df = pd.concat([pc_df,temp_df.loc[temp_df[0] == temp].mean().to_frame().T])
        pc_df_std = pd.concat([pc_df_std,temp_df.loc[temp_df[0] == temp].std().to_frame().T])
    pc_df.reset_index(drop = True, inplace = True)
    pc_df_std.reset_index(drop = True, inplace = True)

    pc_df.rename(columns={0:'temperature'}, inplace = True)
    pc_df_std.rename(columns={0:'temperature'}, inplace = True)
    pc_df_std.temperature = pc_df.temperature

    f, ax = plt.subplots()
    for col in pc_df.columns[1:7]:
        ax.errorbar(pc_df['temperature'],pc_df[col],yerr=pc_df_std[col]*2,fmt='o',markersize=7,capsize=2)
    plt.show()

    return (pc_df,pc_df_std)

def perm_cycle_fit_func(k, A, mu_hat):
    return A * np.exp(-mu_hat*k)/k

def mu_hat_fit_func(T, A, Tc, nu):
    return A * (T-Tc) ** nu

def find_mu_hats(pc_df, pc_df_std):
    num_rows = pc_df.shape[0]
    mu_hats = []
    for i in range(num_rows):
        popt, _ = curve_fit(perm_cycle_fit_func, np.array(pc_df.iloc[i][2:].index,dtype='int'), pc_df.iloc[i][2:].values, p0 = (0,0), bounds= ([0,0],[np.inf, np.inf]))
#        print(pc_df.iloc[i][0],popt)
        mu_hats.append([pc_df.iloc[i][0], popt[1]])
    return np.array(mu_hats,dtype = 'float')

def fit_mu_hats(mu_hats):
    popt, _ = curve_fit(mu_hat_fit_func, mu_hats[16:,0], mu_hats[16:,1], p0 = (1,3,.5))
    return(popt[1])
