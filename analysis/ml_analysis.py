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
import seaborn as sns
from scipy.optimize import curve_fit
from scipy.optimize import brentq
from sklearn import svm
from sklearn.model_selection import GridSearchCV

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
            if(temp_num == 0 and trial == 0):
                test_data = pd.concat([temppc,winding],axis = 1)
            else:
                test_data = pd.concat([test_data,pd.concat([temppc,winding],axis = 1)])

    test_data.to_csv(folder_path+'combined_data/test_data.csv',index=False)

def read_in_data(folder_path,particles):


    main_file_path = folder_path+'combined_data/trial_data.csv'
    data = pd.read_csv(main_file_path)

    for i in range(1,particles+1):
        data[str(i)] = data[str(i)]/particles
    
    data['wx'] = data['wx'].abs()
    data['wy'] = data['wy'].abs()
    data['wz'] = data['wz'].abs()

    print(data.describe())
    print(data.columns)

    return data

def make_model(data):
#    print(check_output(["ls", folder_path]).decode("utf8"))

#    for i in range(15,33):
#        data.drop([str(i)],axis=1, inplace=True)
#    data.drop(['wx','wy','wz'],axis=1, inplace=True)


    num_splits = int(np.floor(len(data.loc[data.temperature==data.temperature.unique()[0]])/200))
    train_data_avgd= []
    train_data = data
    for split in range(num_splits):
        split_data = train_data.sample(frac=1/(num_splits-split))
        train_data = train_data.drop(split_data.index)
        train_data_avgd.append(average_data(split_data))
    avg_data = pd.concat(train_data_avgd)
    avg_data.reset_index(drop=True,inplace=True)
    train_data = avg_data.sample(frac=.8)
    test_data=avg_data.drop(train_data.index)
    train_y = train_data['T/Tc']
    train_X = train_data.drop(['T/Tc','temperature'], axis=1)
    test_y = test_data['T/Tc']
    test_X = test_data.drop(['T/Tc','temperature'], axis=1)

#    model = XGBClassifier(n_estimators = 1000,max_depth=8, learning_rate=0.05)
#    model.fit(train_X, train_y, early_stopping_rounds=10,
#                 eval_set=[(test_X, test_y)], verbose=True)
#    xgb.plot_tree(model)

    model = svm.SVC(kernel='rbf', gamma=.7, C=1, verbose = True)
    model.fit(train_X, train_y)
    predictions = model.predict(test_X)
    print("Mean Absolute Error : " + str(mean_absolute_error(np.array(predictions), test_y)))


    train_data = avg_data.sample(frac=.8)
    test_data = avg_data.drop(train_data.index)
    train_y = train_data['temperature']/3.3125
    train_X = train_data.drop(['T/Tc','temperature'], axis=1)
    test_y = test_data['temperature']/3.3125
    test_X = test_data.drop(['T/Tc','temperature'], axis=1)

    model2 = XGBRegressor(n_estimators = 1000,max_depth=8, learning_rate=0.05)
    model2.fit(train_X, train_y, early_stopping_rounds=10,eval_metric='mae',
                 eval_set=[(test_X, test_y)], verbose=True)

    predictions = model2.predict(test_X)
    print("Mean Absolute Error : " + str(mean_absolute_error(np.array(predictions), test_y)))

    return (model,model2, avg_data)

def visualize_correlations(data):

    data_index = []
    maxn = 5
    
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')

    fig,ax = plt.subplots(ncols=maxn,nrows=maxn,figsize = (10,10))

    colors = np.array(['b','r'])
    cdat = np.take(colors, np.array(data['T/Tc'].values,dtype='int'))
    for i in range(len(ax)):
        for j in range(len(ax[i])):
            ax[i,j].scatter(data[str(i+2)],data[str(j+2)],c=cdat,edgecolor = 'black')
            ax[i,j].annotate(xycoords = 'axes fraction', xy= (0.05,.85), s ='('+str(i+2)+','+str(j+2)+')')
    plt.subplots_adjust(wspace=.35, hspace=.35)

    plt.savefig('corr1.pdf', bbox_inches='tight', pad_inches=0)
    plt.show()

    correlation = data.corr()
    plt.figure(figsize=(25,25))
    cbar_ax = fig.add_axes([.905, .3, .05, .3])
    sns.heatmap(correlation, vmin = -1, vmax=1, cbar_ax = cbar_ax, square=True,annot=True,cmap='viridis')

    plt.savefig('corr2.pdf', bbox_inches='tight', pad_inches=0)
    plt.show()

def average_data(data):
    temps = data.temperature.unique()
    avgd_data = pd.DataFrame([],columns = data.columns)
    for temp in temps:
        avgd_data = avgd_data.append(data.loc[data.temperature == temp].mean().to_frame().T)
    return avgd_data


def predict_temperature(folder_path, model, particles, tr_data, num_parts = 1):
    test = pd.read_csv(folder_path+'combined_data/test_data.csv')
    for i in range(1,particles+1):
        test[str(i)] = test[str(i)]/particles
    test['wx'] = test['wx'].abs()
    test['wy'] = test['wy'].abs()
    test['wz'] = test['wz'].abs()
    
    test_data_avgd = []
    for split in range(num_parts):
        split_data = test.sample(frac=1/(num_parts-split))
        test = test.drop(split_data.index)
        test_data_avgd.append(average_data(split_data))
    test = pd.concat(test_data_avgd)
    test.reset_index(drop=True,inplace=True)

    test_X2 = test[tr_data]
    predicted_ttc = model.predict(test_X2)
    final = pd.DataFrame({'temperature': test.temperature, 'TTc': np.array(predicted_ttc,dtype='float')})

    return final

def scatter_temperature(dfsp):
    for dfp in dfsp:
        plotting = []
        temps = np.round(dfp[1].temperature.unique(),1)
        for temp in temps:
            plotting.append([temp,dfp[1].loc[np.round(dfp[1].temperature,1) == temp].TTc.mean()])
        plotting = np.array(plotting)
        plt.scatter(plotting[:,0],plotting[:,1],label = dfp[0])
        plt.legend()

    plt.show()
