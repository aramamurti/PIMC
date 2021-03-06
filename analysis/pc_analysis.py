####  Calculates critical temperature using permutation cycle data ####

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

    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')


    f, ax = plt.subplots()
    for col in pc_df.columns[1:7]:
        ax.errorbar(pc_df['temperature'],pc_df[col],yerr=pc_df_std[col]*2,fmt='o',markersize=7,capsize=2, label = col)
    plt.legend()
    ax.xaxis.grid(color='gray',linestyle='dashed')
    ax.yaxis.grid(color='gray',linestyle='dashed')

    ax.set_xlabel(r'$T$')
    ax.set_ylabel(r'$\rho_k$')
    plt.show()

    return (pc_df,pc_df_std)

def perm_cycle_fit_func(k, A, mu_hat,alpha):
    return A * np.exp(-mu_hat*k)/k**alpha

def mu_hat_fit_func(T, A, Tc, nu):
    return A * (T-Tc) ** nu

def find_mu_hats(pc_df, pc_df_std, alpha=1.5):
    num_rows = pc_df.shape[0]
    mu_hats = []
    mu_hat_errs = []
    for i in range(num_rows):
        popt, pcov = curve_fit(perm_cycle_fit_func, np.array(pc_df.iloc[i][2:].index,dtype='int'), pc_df.iloc[i][2:].values, p0 = (0,0,alpha), bounds= ([0,0,alpha],[np.inf, np.inf,alpha+.001]))
#        print(pc_df.iloc[i][0],popt)
        mu_hats.append([pc_df.iloc[i][0], popt[1]])
        mu_hat_errs.append([pc_df.iloc[i][0], np.sqrt(np.diag(pcov))[1]])
    return (np.array(mu_hats,dtype = 'float'),np.array(mu_hat_errs,dtype = 'float'))

def fit_mu_hats(mu_hats,mu_hat_errs, start):
    popt, pcov = curve_fit(mu_hat_fit_func, mu_hats[start:,0], mu_hats[start:,1], p0 = (1,3,.5))
    return(popt)

def scatter_mu_hats(mu_hats, mu_hat_errs, fit_vals, name):
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')

    f, ax = plt.subplots()
    ax.errorbar(mu_hats[19:,0],mu_hats[19:,1],yerr = mu_hat_errs[19:,1],fmt='o',markersize=7,capsize=2,c='b')
    xfit = np.linspace(fit_vals[1],np.max(mu_hats[:,0])+.5,100)
    ax.plot(xfit, fit_vals[0]*(xfit-fit_vals[1])**fit_vals[2],c='b')
    ax.xaxis.grid(color='gray',linestyle='dashed')
    ax.yaxis.grid(color='gray',linestyle='dashed')
    ax.set_xlabel(r'$T/T_c$',fontsize=16)
    ax.set_ylabel(r'$\hat{\mu}$',fontsize=16)
    ax.set_xlim(1.5,5.3)

#    plt.savefig(name,dpi=600)
    plt.show()
