# -*- coding: utf-8 -*-
"""
Created on Mon Mar 30 10:06:38 2015

@author: root
"""

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal

from Funktions_Signal import FindImpuls
import pysoundfile

RigtigImpuls ='Lydoptagelser/0000_estimated_ir_ch24_single_mic.wav'

def SinusSweep(StartFrekvens,SlutFrekvens,Varighed):
    """
    Laver et logaritmisk sinus sweep.
    """
    w1 = StartFrekvens
    w2 = SlutFrekvens
    T = Varighed
    SampFrek = 44100
    Tid = np.linspace(0,T,T*SampFrek)
    Sweep = np.sin(1/2*np.pi + 2*np.pi*((w1 * T)/np.log(w2/w1))*(np.exp((Tid/T)*np.log(w2/w1)) - 1)) 
    return Sweep, SampFrek

def LoadFil(Sti):
    Lyddata, Samplerate = pysoundfile.read(Sti)
    return Lyddata[:,0]

def FindOutputNyImpuls(Sweep,Impuls):
    Output = signal.fftconvolve(Sweep,Impuls,mode='full')
    NyImpuls = FindImpuls(Sweep,Output)
    return Output,NyImpuls
      
def Pad(Signal,AntalNul,Padding=False):
    if Padding==True:
        Signal = np.hstack((np.zeros(AntalNul),Signal,np.zeros(AntalNul)))
    return Signal
    
def Forsinkelse(Sweep,Enheder):
    Forsink = np.zeros(Enheder)
    Forsink[Enheder-1] = 1
    return Forsink
    
def Reflektion(Sweep,Enheder,Gain=0.5):
    Reflek = np.zeros(len(Sweep))
    Reflek[100] = 1
    Reflek[Enheder] = 1 * Gain
    return Reflek
    
def LavImpuls(Mode,Sweep,Enheder=5000,Gain=0.5):
    if Mode=='Delay':
        Impuls = Forsinkelse(Sweep,Enheder)
    elif Mode == 'Reflektion':
        Impuls = Reflektion(Sweep,Enheder,Gain)
    elif Mode == 'Rigtig':
        Impuls = LoadFil(RigtigImpuls)
        Impuls = Impuls[:10000]
    return Impuls

def Plotsweep(d1,d2,PlotFrekvens):
    fig = plt.figure()

    gsx,gsy = (4, 3)
    gridspec1 = plt.GridSpec(gsx, gsy)
    gridspec1.update(left=0.05, right=0.47,hspace=0.3)
    ax1 = plt.subplot(gridspec1[:gsx//2,:])
    plt.setp(ax1.get_xticklabels(), visible=False)
    plt.plot(d1[0], d1[1])
    
    ax2 = plt.subplot(gridspec1[gsx//2:,:],sharex=ax1)
    plt.plot(d1[0], d1[2])
    
    gsx,gsy = (4, 3)
    gridspec2 = plt.GridSpec(gsx, gsy)
    gridspec2.update(left=0.53, right=.98,hspace=0.3)
    ax3 = plt.subplot(gridspec2[:gsx//2,:])
    plt.setp(ax3.get_xticklabels(), visible=False)
    plt.plot(d2[0], d2[1])
    if PlotFrekvens==True:
        plt.ylim(0,2)
    ax4 = plt.subplot(gridspec2[gsx//2:,:],sharex=ax3)
    plt.plot(d2[0], d2[2])
#    if PlotFrekvens==False:
#        fig.suptitle(u'Tidsdomænet', size=20)
    if PlotFrekvens==True:
          plt.ylim(0,2)
#        fig.suptitle(u'Frekvensdomænet',size=20)
    plt.setp(ax1, title='Input')
    plt.setp(ax2, title='Output')
    plt.setp(ax3, title='Impuls')
    plt.setp(ax4, title='NyImpuls')
    
def PlotTuples(Sweep,Output,Impuls,NyImpuls,PlotFrekvens=True):
    Sweep = np.hstack((Sweep,np.zeros(len(Output)-len(Sweep))))
    if len(NyImpuls) >= len(Impuls):
        Impuls = np.hstack((Impuls,np.zeros(len(NyImpuls)-len(Impuls))))
    elif len(Impuls) > len(NyImpuls):
        NyImpuls = np.hstack((NyImpuls,np.zeros(len(Impuls)-len(NyImpuls))))

    if PlotFrekvens==True:
        SweepF = np.fft.rfft(Sweep)
        OutputF = np.fft.rfft(Output)
        ImpulsF = np.fft.rfft(Impuls)
        NyImpulsF = np.fft.rfft(NyImpuls)
        XAkse = np.arange(len(SweepF),step=1)
        XAkseImp = np.arange(len(ImpulsF),step=1)
        InOut = (XAkse,np.abs(SweepF),np.abs(OutputF))
        Impulser = (XAkseImp,np.abs(ImpulsF),np.abs(NyImpulsF))
        Plotsweep(InOut,Impulser,True)
    elif PlotFrekvens==False:
        XAkse = np.arange(len(Sweep),step=1)
        XAkseImp = np.arange(len(Impuls),step=1)
        InOut = (XAkse,Sweep,Output)
        Impulser = (XAkseImp,Impuls,NyImpuls)
        Plotsweep(InOut,Impulser,False)
    plt.savefig('Simuleringsplot.pdf')
        

def mls(i, j): 
    """
    Funktionen tager i ind som tids step som skal være lig 
    længden af MLS sekvensen og j som begyndelsesværdien 
    på formen [1,1,0,1] med længden m
    Den genererer en MLS sequens for m = 2,3,4,6,7,15
    """    
    k = 0    
    binseq = [0]*(i+1)
    binseq[0] = j
    while k < i: #laver shift operationen
        binseq[k+1] = list(np.roll(binseq[k],1))        
        binseq[k+1][0] = (binseq[k][-1]+binseq[k][-2])%2
        k = k + 1     
    #udtager sidste element i hver sekvens
    y = []
    for i in binseq:
        y.append(i[-1])
    return np.array(y)
