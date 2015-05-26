# -*- coding: utf-8 -*-
"""
Created on Thu Feb 26 13:38:03 2015

@author: root
"""

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

from scipy import signal


def mls(i, j): 
    """
    Funktionen tager i ind som tids step som skal være lig 
    længden af MLS sekvensen og j som begyndelsesværdien 
    på formen [1,1,0,1] med længden m
    Den genererer en MLS sequens for m = 2,3,4,6,7,15
    
    Parameters
    -------
    i : int
        lenden af MLS sekvensen.
    j : list of int 
        j er begyndelsesværdien på formen [1,1,0,1].
    Returns
    -------
    Sweep : 1D ndarray of floats
        MLS sekvensen i form af et array.

    Eksempel
    -------        
    sweep = mls((2**15)-2, [1,1,1,1,1,1,1,1,1,1,1,1,1,1,0])        
    """
    k = 0    
    
    binseq = [0]*(i+1)
    binseq[0] = j
 
    while k < i: # laver shift operationen
        binseq[k+1] = list(np.roll(binseq[k],1))        
        binseq[k+1][0] = (binseq[k][-1]+binseq[k][-2])%2

        k = k + 1 
        
    # udtager sidste element i hver sekvens
    y = []
    for i in binseq:
        y.append(i[-1])
    
    return np.array(y)

def SinusSweep(StartFrekvens,SlutFrekvens,Varighed):
    """
    Laver et eksponentiel sinus sweep.
    
    Parameters
    -------
    StartFrekvens : float
        Frekvensen som sinussweepet begynder ved.
    SlutFrekvens : float
        Frekvensen som sinussweepet afslutter ved.
    Varighed : float
        Antal sekunder sinussweepet skal vare.
    Returns
    -------
    Sweep : 1D ndarray of floats
        Sinussweepet i form af et array.
    SampFrek : int
        Samplingsfrekvensen til sweepet
    """
    w1 = StartFrekvens
    w2 = SlutFrekvens
    T = Varighed
    SampFrek = 44100
    Tid = np.linspace(0,T,T*SampFrek)
    Sweep = np.sin(2*np.pi*((w1 * T)/np.log(w2/w1))*(np.exp((Tid/T)*np.log(w2/w1)) - 1)) 
    return Sweep, SampFrek

def SinusSweepSubsample(StartFrekvens,SlutFrekvens,Varighed,Delay=0):
    """
    Laver et eksponentiel sinus sweep med et subsample delay.
    
    Parameters
    -------
    StartFrekvens : float
        Frekvensen som sinussweepet begynder ved.
    SlutFrekvens : float
        Frekvensen som sinussweepet afslutter ved.
    Varighed : float
        Antal sekunder sinussweepet skal vare.
    Delay : float
        Et subsample delay i intervallet [0,1].
    Returns
    -------
    Sweep : 1D ndarray of floats
        Sinussweepet i form af et array.
    SampFrek : int
        Samplingsfrekvensen til sweepet
    """
    if Delay <0 or Delay > 1:
        raise 'Not implemented'
    w1 = StartFrekvens
    w2 = SlutFrekvens
    T = Varighed
    SampFrek = 44100
    Tid = np.linspace(0,T,T*SampFrek)-Delay/SampFrek
    Tid[np.where(Tid<0)] = 0
    Sweep = np.sin(2*np.pi*((w1 * T)/np.log(w2/w1))*(np.exp((Tid/T)*np.log(w2/w1)) - 1))
    return Sweep, SampFrek

def LinSweep(StartFrekvens,SlutFrekvens,Varighed):
    """
    Laver et lineært sinus sweep.
    
    Parameters
    -------
    StartFrekvens : float
        Frekvensen som sinussweepet begynder ved.
    SlutFrekvens : float
        Frekvensen som sinussweepet afslutter ved.
    Varighed : float
        Antal sekunder sinussweepet skal vare.
    Returns
    -------
    Sweep : 1D ndarray of floats
        Sinussweepet i form af et array.
    SampFrek : int
        Samplingsfrekvensen til sweepet
    """
    w1 = StartFrekvens
    w2 = SlutFrekvens
    T = Varighed
    SampFrek = 44100
    Tid = np.linspace(0,T,T*SampFrek)
    k = (w2-w1)/T
    Frekvenser = w1+k*Tid/2
    Sweep = np.sin(1/2 * np.pi + 2*np.pi*(Frekvenser)*Tid)
    return Sweep,SampFrek



"""
Her udføres test på sweep for at se om de virker
- kortere sweep
- kun nogle frekvenser
"""    
if __name__ == '__main__':
    pass
    a=mls((2**15)-2, [1,1,1,1,1,1,1,1,1,1,1,1,1,1,0])
#    a =  mls(32766, [1,1,1,1,1,1,1,1,1,1,1,1,1,1,0]) 
    plt.plot(np.fft.rfft(a))
    b=signal.max_len_seq(16)