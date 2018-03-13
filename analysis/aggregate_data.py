# -*- coding: utf-8 -*-
"""
    Created on Mon Jul 27 15:59:10 2015
    
    @author: adith
    """

import scipy
import numpy
import csv
import os, os.path
import re
import sys

pi = scipy.pi

def represents_int(s):
    try:
        int(s)
        return True
    except ValueError:
        return False

def aggregate_data(folder_path, output_folder):
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
        num_temps.append(len([name for name in os.listdir(folder_path+'{}/output/'.format(i+1)) if (os.path.isfile(folder_path+'{}/output/'.format(i+1)+name) and "overview" in name)]))

    num_temps = set(num_temps)

    if len(num_temps) == 1:
        
        num_temps = num_temps.pop()
        
        enlist = [0]*num_temps
        stdenlist = [0]*num_temps
        permlist = [[0]*particles]*num_temps
        wnx = [[0]*10]*num_temps
        wny = [[0]*10]*num_temps
        wnz = [[0]*10]*num_temps

        wn = [[0]*10]*num_temps

        prefix = folder_path+'{}/output/'.format(1)
        temps = []
        for num in range(0,num_temps):
            fo = open(prefix+'overview_{}.txt'.format(num),'r')
            t = '';
            for line in fo:
                if "Temperature" in line:
                    t+=line
            temps.append(float(re.findall(r"[-+]?\d*\.\d+|\d+", t)[0]))
            fo.close()

        temp_ind = [[] for x in range(num_temps)]
        
        for i in range(num_trials):
            prefix = folder_path+ '{}/output/'.format(i+1)
            for num in range(0,num_temps):
                fo = open(prefix + 'overview_{}.txt'.format(num),'r')
                t = '';
                for line in fo:
                    if "Temperature" in line:
                        t+=line
                temp_ind[num].append(temps.index(float(re.findall(r"[-+]?\d*\.\d+|\d+", t)[0])))
                fo.close()
    
        for trial in range(num_trials):
            for num in range(0,num_temps):
                
                tempen = []
                temppc = []
                tempwnx = [[0] for x in range(10)]
                tempwny = [[0] for x in range(10)]
                tempwnz = [[0] for x in range(10)]

                prefix = folder_path +'{}/output/'.format(trial+1)
                
                if(trial == 0):
                    file_num = num
                else:
                    file_num = temp_ind[num][trial-1]
                fe = open(prefix + 'energy_data_{}.csv'.format(file_num),'r')
                fp = open(prefix + 'permutation_data_{}.csv'.format(file_num),'r')
                fw = open(prefix + 'winding_data_{}.csv'.format(file_num),'r')

                reader = csv.reader(fe)
                next(reader)
                for row in reader:
                    if(len(row)>1):
                        tempen.append( (float)(row[1]))

                reader = csv.reader(fp)
                for row in reader:
                    pcrow = list(map(int,row[1::]))
                    temppc.append(pcrow)
                
                reader = csv.reader(fw)
                for row in reader:
                    tempwnx[abs(int(row[1]))][0] += 1
                    tempwny[abs(int(row[2]))][0] += 1
                    tempwnz[abs(int(row[3]))][0] += 1
                
                fe.close()
                fp.close()
                fw.close()

                tempenarr = numpy.array(tempen)
                enlist[num] = (numpy.average(tempenarr)).tolist()
                stdenlist[num] = ((numpy.std(tempenarr)).tolist())/numpy.sqrt(len(tempenarr))

                temppcarr = numpy.array(temppc)
                pctot =  numpy.sum(temppcarr, 0)
                permlist[num]= (pctot / (float) (numpy.sum(pctot))).tolist()

                tempwxarr = numpy.array(tempwnx)
                tempwyarr = numpy.array(tempwny)
                tempwzarr = numpy.array(tempwnz)

                tempwarr = numpy.zeros(10)
                for i in range(len(tempwarr)):
                    tempwarr[i] = tempwxarr[i]+tempwyarr[i]+tempwzarr[i]

                wnx[num] = tempwxarr/(float)(numpy.sum(tempwxarr))
                wny[num] = tempwyarr/(float)(numpy.sum(tempwyarr))
                wnz[num] = tempwzarr/(float)(numpy.sum(tempwzarr))
                wn[num] = tempwarr/(float)(numpy.sum(tempwarr))

            wnx = numpy.array(wnx)
            wny = numpy.array(wny)
            wnz = numpy.array(wnz)
            wn = numpy.array(wn)

            dFx = [0]*num_temps
            dFy = [0]*num_temps
            dFz = [0]*num_temps

            dF = [0]*num_temps

            for T in range(len(wnx[:,0])):
                for wp in range(len(wnx[T])):
                    dFx[T] += numpy.real(wnx[T,wp][0]*numpy.exp(1j*pi*wp))*1
                    dFy[T] += numpy.real(wny[T,wp][0]*numpy.exp(1j*pi*wp))*1
                    dFz[T] += numpy.real(wnz[T,wp][0]*numpy.exp(1j*pi*wp))*1
                    dF[T] += numpy.real(wn[T,wp]*numpy.exp(1j*pi*wp))*1
                dFx[T] = -temps[T]*numpy.log(dFx[T])
                dFy[T] = -temps[T]*numpy.log(dFy[T])
                dFz[T] = -temps[T]*numpy.log(dFz[T])
                dF[T] =  -temps[T]*numpy.log(dF[T])

            we = open(folder_path+output_folder+'energy_data_{}.csv'.format(trial+1),'w')
            cwe = csv.writer(we)
            cwe.writerow(temps)
            cwe.writerow(enlist)
            cwe.writerow(stdenlist)

            wp = open(folder_path+output_folder+'permutation_data_{}.csv'.format(trial+1),'w')

            cwp = csv.writer(wp)
            for i in range(len(temps)):
                cwp.writerow([temps[i]]+permlist[i])

            ww = open(folder_path+output_folder+'winding_data_{}.csv'.format(trial+1),'w')
            cww = csv.writer(ww)
            for i in range(len(temps)):
                cww.writerow([temps[i]]+wn[i].tolist())


            wf = open(folder_path+output_folder+'free_energy_data_{}.csv'.format(trial+1),'w')
            cwf = csv.writer(wf)
            cwf.writerow(temps)
            cwf.writerow(dFx)
            cwf.writerow(dFy)
            cwf.writerow(dFz)
            cwf.writerow(dF)


            we.close()
            wp.close()
            wf.close()
            ww.close()
