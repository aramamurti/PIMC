# -*- coding: utf-8 -*-
"""
    Created on Wed. Feb 1 11:39:15 2018
    
    @author: adith
    """


####  STILL PRELIMINARY, NEEDS MORE WORK ####


import scipy as sp
import numpy as np
import pandas as pd
import sys
import os
import re
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import Imputer
import xgboost as xgb
from xgboost import XGBClassifier
from xgboost import XGBRegressor
from collections import Counter
from subprocess import check_output
from sklearn.metrics import mean_absolute_error
from scipy.interpolate import UnivariateSpline
from sklearn import svm

from scipy.optimize import curve_fit
from scipy.optimize import brentq


def sigmoid(x, x0, k):
     y = 1 / (1 + np.exp(-k*(x-x0)))
     return y

def sigmoidoff(x, x0, k):
     y = 1 / (1 + np.exp(-k*(x-x0)))
     return y - 0.65

def create_input(folder_path):


    def represents_int(s):
        try:
            int(s)
            return True
        except ValueError:
            return False

    num_trials = len([name for name in os.listdir(folder_path) if (os.path.isdir(folder_path + name) and represents_int(name))])


    fo = open(folder_path+'1/output/overview_0.txt','r')
    p = '';
    for line in fo:
        if "Particles" in line:
            p+=line
    particles = int(re.findall(r'\d+',p)[0])
    fo.close()


    num_temps =[]
    for i in range(num_trials):
        num_temps.append(len([name for name in os.listdir(folder_path+'{}/output/'.format(i+1)) if (os.path.isfile(folder_path+'{}/output/'.format(i+1) + name) and "overview" in name)]))

    num_temps = set(num_temps)
    
    if(len(num_temps) == 1):
        num_temps = num_temps.pop()
    else:
        sys.exit()

    prefix = folder_path+'{}/output/'.format(1)
    temps = []
    for num in range(0,num_temps):
        fo = open(prefix + 'overview_{}.txt'.format(num),'r')
        t = '';
        for line in fo:
            if "Temperature" in line:
                t+=line
        temps.append(float(re.findall(r"[-+]?\d*\.\d+|\d+", t)[0]))
        fo.close()

    for temp_num in range(num_temps):
        for trial in range(num_trials):
            prefix = folder_path+'{}/output/'.format(trial+1)
            temppc = pd.read_csv(prefix+'permutation_data_{}.csv'.format(temp_num),header=None)
            temp = temps[temp_num]
            if(temp/3.3125 > 1):
                ttc = 1
            else:
                ttc = 0
            temppc = temppc.drop(columns = 0)
            temppc['temperature'] = np.full(len(temppc[1]),temp)
            temppc['T/Tc'] = np.full(len(temppc[1]),ttc)
            winding = pd.read_csv(prefix + 'winding_data_{}.csv'.format(temp_num),header=None)
            winding = winding.drop(columns = 0)
            winding.columns = ['wx','wy','wz']
            winding['w2'] = (winding['wx']**2+winding['wy']**2+winding['wz']**2)/3
#            winding = winding.drop(columns = ['wx','wy','wz'])
            if(temp_num == 0 and trial == 0):
                trial_data = pd.concat([temppc,winding],axis = 1)
            else:
                trial_data = pd.concat([trial_data,pd.concat([temppc,winding],axis = 1)])

    print(trial_data.columns)
    trial_data.to_csv(folder_path+'combined_data/trial_data.csv',index=False)

def create_test(folder_path):
    
    def represents_int(s):
        try:
            int(s)
            return True
        except ValueError:
            return False

    num_trials = len([name for name in os.listdir(folder_path) if (os.path.isdir(folder_path + name) and represents_int(name))])

    fo = open(folder_path+'1/output/overview_0.txt','r')
    p = '';
    for line in fo:
        if "Particles" in line:
            p+=line
    particles = int(re.findall(r'\d+',p)[0])
    fo.close()


    num_temps =[]
    for i in range(num_trials):
        num_temps.append(len([name for name in os.listdir(folder_path+'{}/output/'.format(i+1)) if (os.path.isfile(folder_path+'{}/output/'.format(i+1) + name) and "overview" in name)]))

    num_temps = set(num_temps)
    
    if(len(num_temps) == 1):
        num_temps = num_temps.pop()
    else:
        sys.exit()

    prefix = folder_path+'{}/output/'.format(1)
    temps = []
    for num in range(0,num_temps):
        fo = open(prefix + 'overview_{}.txt'.format(num),'r')
        t = '';
        for line in fo:
            if "Temperature" in line:
                t+=line
        temps.append(float(re.findall(r"[-+]?\d*\.\d+|\d+", t)[0]))
        fo.close()

    for temp_num in range(num_temps):
        for trial in range(num_trials):
            prefix = folder_path+'{}/output/'.format(trial+1)
            temppc = pd.read_csv(prefix + 'permutation_data_{}.csv'.format(temp_num),header=None)
            temp = temps[temp_num]
            temppc = temppc.drop(columns = 0)
            temppc['temperature'] = np.full(len(temppc[1]),temp)
            winding = pd.read_csv(prefix + 'winding_data_{}.csv'.format(temp_num),header=None)
            winding = winding.drop(columns = 0)
            winding.columns = ['wx','wy','wz']
            winding['w2'] = (winding['wx']**2+winding['wy']**2+winding['wz']**2)/3
#            winding = winding.drop(columns = ['wx','wy','wz'])
            if(temp_num == 0 and trial == 0):
                test_data = pd.concat([temppc,winding],axis = 1)
            else:
                test_data = pd.concat([test_data,pd.concat([temppc,winding],axis = 1)])

    test_data.to_csv(folder_path+'combined_data/test_data.csv',index=False)

def read_in_data(folder_path,particles):


    main_file_path = folder_path+'combined_data/trial_data.csv'
    data = pd.read_csv(main_file_path)

    for i in range(1,particles+1):
        data[str(i)] = data[str(i)]/i

    print(data.describe())
    print(data.columns)

    return data

def make_model(data):
#    print(check_output(["ls", folder_path]).decode("utf8"))

    for i in range(16,33):
        data = data.drop([str(i)],axis=1)
    data = data.drop(['wx','wy','wz'],axis=1)
    train_data=data.sample(frac=0.8)
    test_data=data.drop(train_data.index)
    train_y = train_data['T/Tc']
    train_X = train_data.drop(['T/Tc','temperature'], axis=1)
    test_y = test_data['T/Tc']
    test_X = test_data.drop(['T/Tc','temperature'], axis=1)


    model = XGBClassifier(n_estimators = 1000,max_depth=10, learning_rate=0.05)
    model.fit(train_X, train_y, early_stopping_rounds=30,
                 eval_set=[(test_X, test_y)], verbose=True)


    predictions = model.predict(test_X)
    print("Mean Absolute Error : " + str(mean_absolute_error(np.array(np.round(predictions)), test_y)))

#    xgb.plot_tree(model)

#    model = svm.SVC(kernel='rbf', gamma=0.7, C=1.0, verbose = True)
#    model.fit(train_X, train_y)


    train_data=data.sample(frac=0.8)
    test_data=data.drop(train_data.index)
    train_y = train_data['temperature']/3.3125
    train_X = train_data.drop(['T/Tc','temperature'], axis=1)
    test_y = test_data['temperature']/3.3125
    test_X = test_data.drop(['T/Tc','temperature'], axis=1)

    model2 = XGBRegressor(n_estimators = 1000,max_depth=10, learning_rate=0.05)
    model2.fit(train_X, train_y, early_stopping_rounds=30,
                 eval_set=[(test_X, test_y)], verbose=True)

    predictions = model2.predict(test_X)
    print("Mean Absolute Error : " + str(mean_absolute_error(np.array(np.round(predictions)), test_y)))

    return (model,model2, train_X.columns)

def visualize_data(data):

    data_index = []
    maxn = 6

    for i in range(2,9):
        data_index.append(np.sort(data[str(i)].unique()))
    plot_data = []
    for i in range(0,maxn):
        for j in range(0,maxn):
            data_int = np.zeros((int(np.max(data_index[i])+1),int(np.max(data_index[j])+1)))
            for k in data_index[i]:
                for l in data_index[j]:
                    data_int[int(k),int(l)] = data['T/Tc'][data.loc[(data[str(i+2)] == k) & (data[str(j+2)] == l)].index].mean()
            plot_data.append(data_int)


    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')

    fig,ax = plt.subplots(ncols=maxn,nrows=maxn,sharey=True,sharex=True)

    for i in range(len(ax)):
        for j in range(len(ax[i])):
            c = ax[i,j].pcolormesh(plot_data[maxn*i+j],vmin=0, vmax=1)
            ax[i,j].annotate('('+str(i+2)+', '+str(j+2)+')',xy= (6,10),)

    ax[0,0].set_xticks([0,2,4,6,8,10])
    ax[0,0].set_xlim((0,12))
    ax[0,0].set_yticks([0,2,4,6,8,10])
    ax[0,0].set_ylim((0,12))

    plt.subplots_adjust(wspace=0, hspace=0)
    plt.subplots_adjust(right=0.8)

    cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    cbar = plt.colorbar(c,cax=cbar_ax)

    plt.show()


def predict_condensed(folder_path, model, particles, tr_data):
    test = pd.read_csv(folder_path+'combined_data/test_data.csv')
    for i in range(1,particles+1):
        test[str(i)] = test[str(i)]/i

    test_X2 = test[tr_data]
    predicted_ttc = model.predict(test_X2)
    final = pd.DataFrame({'temperature': test.temperature, 'TTc': np.array(predicted_ttc,dtype='float')})

    return final

def predict_temperature(folder_path, model, particles, tr_data):
    test = pd.read_csv(folder_path+'combined_data/test_data.csv')
    for i in range(1,particles+1):
        test[str(i)] = test[str(i)]/i

    test_X2 = test[tr_data]
    predicted_ttc = model.predict(test_X2)
    final = pd.DataFrame({'temperature': test.temperature, 'TTc': np.array(predicted_ttc,dtype='float')})

    return final

def scatter_condensed(dfsp):
    for dfp in dfsp:
        plotting = []
        temps = dfp[1].temperature.unique()
        for temp in temps:
            plotting.append([temp,dfp[1].loc[dfp[1].temperature == temp].TTc.mean()])
        plotting = np.array(plotting)
        plt.scatter(plotting[:,0],plotting[:,1],label = dfp[0])
        plt.legend()

    plt.show()
