# -*- coding: utf-8 -*-
"""
    Created on Mon Jul 27 15:59:10 2015
    
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
from collections import Counter
from subprocess import check_output
from sklearn.metrics import mean_absolute_error
from scipy.interpolate import UnivariateSpline
from sklearn.naive_bayes import GaussianNB

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

def make_model(folder_path,particles):
#    print(check_output(["ls", folder_path]).decode("utf8"))

    main_file_path = folder_path+'combined_data/trial_data.csv'
    data = pd.read_csv(main_file_path)

#    for i in range(1,particles+1):
#        data[str(i)] = data[str(i)].i

    print(data.describe())
    print(data.columns)

    train_data=data.sample(frac=0.8,random_state=200)
    test_data=data.drop(train_data.index)

    train_y = train_data['T/Tc']
    train_X = train_data.drop(['T/Tc','temperature'], axis=1)
    test_y = test_data['T/Tc']
    test_X = test_data.drop(['T/Tc','temperature'], axis=1)

    my_model = XGBClassifier(n_estimators = 1000,max_depth=8)
    my_model.fit(train_X, train_y, early_stopping_rounds=5,
                 eval_set=[(test_X, test_y)], verbose=True)

#    xgb.plot_tree(my_model)

    predictions = my_model.predict(test_X)

    print("Mean Absolute Error : " + str(mean_absolute_error(np.array(np.round(predictions)), test_y)))

    return my_model

def predict_tc(folder_path, my_model, particles):
    test = pd.read_csv(folder_path+'combined_data/test_data.csv')
#    for i in range(1,particles+1):
#        test[str(i)] = test[str(i)]/i

    test_X2 = test.drop(['temperature'],axis=1)
    predicted_ttc = my_model.predict(test_X2)
    final = pd.DataFrame({'temperature': test.temperature, 'TTc': np.array(predicted_ttc,dtype='float')})

    temps = final.temperature.unique()

    average_ttc = np.zeros((len(temps),2),dtype='float')
    for index,temp in enumerate(temps):
        average_ttc[index] = [temp, final.loc[final.temperature == temp].TTc.mean()]

    popt, pcov = curve_fit(sigmoid, average_ttc[:,0], average_ttc[:,1])

    ax.scatter(average_ttc[:,0], sigmoid(average_ttc[:,0],popt[0],popt[1]))
    root = brentq(sigmoidoff, 1, 5, args=(popt[0],popt[1]))
    return root
