# -*- coding: utf-8 -*-
"""
Created on Mon Mar 16 09:06:36 2015

@author: root
"""

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

from scipy import signal
from pysoundcard import Stream, continue_flag, complete_flag 
import time

import Funktions_Signal as fs
import Funktions_Soundcard as sc
import Sweeps as sweeps


logfil='FindRumImpuls.txt'
def log(s): #write to log
    with open(logfil,'a') as f:
        f.write(s)
        
Tidspunkt = time.strftime('%d-%m_%H-%M-%S')
log('time: '+Tidspunkt+'\n')

input_device, output_device = sc.pysoundcard_test()
input_kanaler = input_device['input_channels']
output_kanaler = output_device['output_channels']


def Get_sweep(SweepType='LogSinus',Tid=1,StartFrek=1,
             SlutFrek=20000):

    if SweepType=='LogSinus':#eksponentiel sinus sweep
        Sweep,Samp = sweeps.SinusSweep(StartFrek,SlutFrek,Tid)
    elif SweepType=='LinSinus':
        Sweep,Samp = sweeps.LinSweep(StartFrek,SlutFrek,Tid)
    elif SweepType=='MLS':
        Sweep = sweeps.mls((2**15)-2, [1,1,1,1,1,1,1,1,1,1,1,1,1,1,0])
    else:
        raise ValueError('Invalid Sweep type')

    Sweep = np.hstack([np.zeros(44100),Sweep,np.zeros(44100*1)])
    Sweep32 = np.vstack(output_kanaler*[Sweep]).T
    return Sweep,Sweep32


def callback(input_data, num_frames, time_info, status_flags):
    """
    Afspiller et sweep og optager responsen.
    """
    global Optagelse,Sweep32,i
    Optagelse = np.append(Optagelse,input_data,axis=0)
    output_data=Sweep32[i*num_frames:(i+1)*num_frames]
    i+=1
    flag = continue_flag if i < i_stop else complete_flag 
    return (output_data, flag)
    i-=1
    return (input_data,continue_flag)

def RumImpuls(Sweep32,block_length=1024):
    global Optagelse,i,i_stop
    i = 0 
    i_stop=len(Sweep32)//block_length # udregner hvor mange hele block længder der  er
    
    Optagelse = np.zeros((1,input_kanaler))    
    
    s = Stream(callback=callback, block_length=block_length)
    s.start()
    #time.sleep((len(Sweep)//44100)+1)
    #s.stop()
    
    while i < i_stop:
        time.sleep(0.1)
    time.sleep(0.2)
    
    return Optagelse

#==============================================================================
# Optagelse
#==============================================================================

SweepType='LogSinus' #'MLS'
Tid=1
StartFrek=1
SlutFrek=20000
block_length=1024

log('SweepType:%s\nTid:%s\nStartFrek:%s\nSlutFrek:%s \n'  
%(SweepType,Tid,StartFrek,SlutFrek))

Sweep,Sweep32 = Get_sweep(SweepType=SweepType,Tid=Tid,
                         StartFrek=StartFrek,SlutFrek=SlutFrek)


Optagelse=RumImpuls(Sweep32,block_length=block_length)

Impuls = fs.FindImpuls(Sweep[:len(Optagelse)],Optagelse[:,0])
ind = Optagelse[:,1]
ud = Optagelse[:,0]
Impuls2 = fs.FindImpuls(ind,ud)

#==============================================================================
# Resultater Gem data
#==============================================================================

gem_data=True

if gem_data:
    SweepSti = "-".join(('Data/Sweep',Tidspunkt))
    OptagelseSti = "-".join(('Data/Optagelse',Tidspunkt))
    #ImpulsSti = string.join(('Data/Impuls',Tidspunkt))
    
    np.savetxt(SweepSti,Sweep)
    np.savetxt(OptagelseSti,Optagelse)
    #np.savetxt(ImpulsSti,Impuls)
log(str(np.argmax(np.abs(Impuls))))
log('\n')    
log(str(np.argmax(np.abs(Impuls2))))
log('\n')  
log('V: -350')#grader
log('\n\n') #write to log finalized'
