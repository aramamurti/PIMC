####  STILL PRELIMINARY, NEEDS MORE WORK ####


import scipy as sp
import numpy as np
import pandas as pd
import sys
import os
import re
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline
from scipy.optimize import curve_fit
from scipy.optimize import brentq


def read_input(wn_file_path):
    num_wns = len([name for name in os.listdir(wn_file_path) if (os.path.isfile(wn_file_path + name) and "winding" in name)])

    wn_dfs = []
    for i in range(num_wns):
        wn_dfs.append(pd.read_csv(wn_file_path+'winding_data_{}.csv'.format(i+1),header=None))
    temp_df = pd.concat(wn_dfs)
    temperatures = temp_df[0].unique()

    wn_df = pd.DataFrame(columns=temp_df.columns)
    wn_df_std = pd.DataFrame(columns=temp_df.columns)
    for temp in temperatures:
        wn_df = pd.concat([wn_df,temp_df.loc[temp_df[0] == temp].mean().to_frame().T])
        wn_df_std = pd.concat([wn_df_std,temp_df.loc[temp_df[0] == temp].std().to_frame().T])
    wn_df.reset_index(drop = True, inplace = True)
    wn_df_std.reset_index(drop = True, inplace = True)

    wn_df.columns = np.array(temp_df.columns,dtype=int) - 1
    wn_df_std.columns = np.array(temp_df.columns,dtype=int) - 1
    wn_df.rename(columns={-1:'temperature'}, inplace = True)
    wn_df_std.rename(columns={-1:'temperature'}, inplace = True)
    wn_df_std.temperature = wn_df.temperature

    return (wn_df,wn_df_std)

def calc_superfluid_frac(wn_df, wn_df_std):
    psp = []
    psp_err = []
    for i in range(len(wn_df)):
        w2 = 0
        w2_err = 0
        for wp in range(0, wn_df.columns[-1]):
            w2 += wn_df.iloc[i][wp] * wp**2
            w2_err += wn_df_std.iloc[i][wp]**2 * wp**4
        w2_err = np.sqrt(w2_err)
        psp.append([wn_df.iloc[i].temperature,wn_df.iloc[i].temperature*w2])
        psp_err.append([wn_df.iloc[i].temperature,wn_df.iloc[i].temperature*w2_err])
    
    psp = np.array(psp)
    psp_err = np.array(psp_err)
    spl = UnivariateSpline(psp[:,0],psp[:,1], w = 1/(np.maximum(np.ones(len(psp_err))*.005,psp_err[:,1])))
    return (psp,psp_err,spl)

def plot_psp(psps, psp_errs,spls,num_particles,name):
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')

    f,ax = plt.subplots()
    ax.xaxis.grid(color='gray',linestyle='dashed')
    ax.yaxis.grid(color='gray',linestyle='dashed')

    ax.set_xlabel(r'$T/T_c$',fontsize=16)
    ax.set_ylabel(r'$\rho_s/\rho$',fontsize=16)
    
    colors = ['r','g','b']
    for i,psp in enumerate(psps):
        ax.errorbar(psp[:,0],psp[:,1],yerr=psp_errs[i][:,1],fmt='o',markersize=7,capsize=2,c=colors[i], label=num_particles[i][:-1])
        ax.plot(psp[:,0],spls[i](psp[:,0]),c=colors[i])
    plt.legend(fontsize=14)
#    plt.savefig(name, dpi=600)
    plt.show()

def find_interesect(psps,spls):
    diffs = []
    temps = []
    for i in range(len(psps)):
        for j in range(i+1,len(psps)):
            diffs.append(spls[i](psps[i][:,0])-spls[j](psps[j][:,0]))
            ind = np.argmax(diffs[-1]<0)
            if(ind > 0 and ind<len(psps[0])-1):
                temps.append((psps[0][ind,0]*abs(diffs[-1][ind-1])+psps[0][ind-1,0]*abs(diffs[-1][ind]))/(abs(diffs[-1][ind])+abs(diffs[-1][ind-1])))
            else:
                temps.append(psps[0][ind,0])
    return np.average(np.array(temps))
