# -*- coding: utf-8 -*-
"""
Created on Tue Mar  3 10:38:05 2015

@author: Groot
"""

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
from Funktions_Signal import FindImpuls, Punkt2Afstand, Udfold, Interpolation, latexify

def LoadLogFile(fileToOpen) :
    """
    Gem information fra log fil i en liste. 
    
    Parameters
    -------
    fileToOpen : string
        Logfil der ønskes åbnet. 
    Returns
    -------
    DataBlocks : list of dataBlock
        Liste med dataBlock, som indeholder forskellige informationer
        om de forskellige målinger. 
    """
    dataBlocks = []
    dataBlock = []

    with open(fileToOpen) as fp:
        for line in fp:
            if line.startswith('\n') : 
                dataBlocks.append(dataBlock)
                dataBlock = []
            else : 
                dataBlock.append(line.replace('\n', ''))

    return dataBlocks

def PlotIRwithDelay(Sweep, y, mode="tid", color='', lab='') :
    """
    Plot Impulsrespons med computerens forsinkelse.  
    
    Parameters
    -------
    Sweep : 1D ndarray of floats
        Det genererede sweep med forsinkelse fra computeren. 
    y : 1D ndarray of floats
        Det optagede signal fra mikrofonen.
    """
    Sweep = Sweep[:len(y)]
    IRwithDelay = FindImpuls(Sweep, y)

    if mode == 'tid' :
        plt.plot(IRwithDelay, color, label=lab)
        plt.xlim(6000, 11000)
        plt.xlabel('Samples')
        plt.ylabel('Amplitude')

def PlotIRwithoutDelay(Sweep, xpc, y, mode="") :
    """
    Plot Impulsrespons uden computerens forsinkelse.  
    
    Parameters
    -------
    Sweep : 1D ndarray of floats
        Det genererede sweep med forsinkelse fra computeren. 
    y : 1D ndarray of floats
        Det optagede signal fra mikrofonen.
    xpc : 1D ndarray of floats
        Det afspillede signal med forsinkelse fra OSDAC.
    """
    xpc = xpc[:len(y)]
    IRwithoutDelay = FindImpuls(xpc, y)    

    if mode == 'tid' :
        plt.plot(IRwithoutDelay)
        plt.xlim(0, 6000)
        plt.xlabel('Samples')
        plt.ylabel('Amplitude')

def PlotDACOSdelay(Sweep, xpc, mode='', color=None) :
    """
    Plot forsinkelsen i computeren.
    
    Parameters
    -------
    Sweep : 1D ndarray of floats
        Det genererede sweep med forsinkelse fra computeren. 
    xpc : 1D ndarray of floats
        Det afspillede signal med forsinkelse fra OSDAC.
    """
    Sweep = Sweep[:len(xpc)]
    pcImpuls = FindImpuls(Sweep, xpc)

    if mode == 'tid' :
        pcImpulsFrekvens = np.fft.fft(pcImpuls)
        freq = np.fft.fftfreq(len(pcImpulsFrekvens), 1/44100.)
        pcImpulsFrekvens[np.where(np.abs(freq) > 20000)] = 0
        plt.plot(np.fft.ifft(pcImpulsFrekvens), color)
        plt.xlim(0, 10000)
        plt.xlabel('Samples')
        plt.ylabel('Amplitude')
    elif mode == 'frekvens' :
        pcImpulsFrekvens = np.fft.rfft(pcImpuls)
        freq = np.fft.rfftfreq(len(pcImpulsFrekvens) * 2 -1, 1/44100.)
        plt.plot(freq, np.abs(pcImpulsFrekvens), 'r')
        pcImpulsFrekvens[np.where(freq > 20000)] = 0
        plt.plot(freq, np.abs(pcImpulsFrekvens), 'b')
        plt.xlim(0, 22050)
        plt.xlabel('Frekvens')
        plt.ylabel('Amplitude')

def PlotLoudspeakerBlock(mode='') :
    """
    Plotter højtaleren i frekvens eller tidsdomænet. 
    """
    with np.load('polaer.npz') as data :
        hoejtalerFrekvens = data['vinkler'].item()['0']

    if mode == 'tid' :
        freq = np.fft.rfftfreq(len(hoejtalerFrekvens) * 2 -1, 1/44100.)
        delay = np.concatenate((np.zeros(191),np.ones(1),np.zeros(200000)),axis=0)
        delay = delay[:len(hoejtalerFrekvens) * 2 -1]
        delayf = np.fft.rfft(delay)
        hoejtalerFrekvens = hoejtalerFrekvens / delayf
        hoejtalerFrekvens[np.where(freq > 20000)] = 0
        hoejtaler = np.fft.irfft(hoejtalerFrekvens)
        plt.plot(hoejtaler)
        plt.xlim(0, 50)
        plt.xlabel('Samples')
        plt.ylabel('Amplitude')
    elif mode == 'frekvens' :
        freq = np.fft.rfftfreq(len(hoejtalerFrekvens) * 2 -1, 1/44100.)
        plt.plot(freq, np.abs(hoejtalerFrekvens), 'r')
        hoejtalerFrekvens[np.where(freq > 20000)] = 0
        plt.plot(freq, np.abs(hoejtalerFrekvens), 'b')
        plt.xlim(0, 22050)
        plt.ylim(0, 0.4)
        plt.xscale('log')
        plt.xlabel('Frekvens')
        plt.ylabel('Amplitude')

def PlotRoom(mode='tid', HOJ=None) :
    """
    Plotter rummets impulsrespons i med alt affoldet, sp kun rummet er tilbage. 
    """
    rum , freq , H0  = Udfold(y, xpc, HOJ=HOJ, mode='full')

    if mode == 'tid' :
        plt.plot(rum)
        plt.xlabel('Samples')
        plt.ylabel('Amplitude')

def Diracfunction() :
    """
    Plotter knoeckers dirac funktion. 
    """
    dirac = np.concatenate((np.zeros(199),np.ones(1),np.zeros(800)),axis=0)
    plt.plot(dirac)
    plt.xlabel('Samples')
    plt.ylabel('Amplitude')
    plt.xlim(0, len(dirac))

def BallonImpuls() :
    """
    Plotter impuls fra ballon i det lyddøde rum. 
    """
    ballon = np.loadtxt("Data/Ballon")
    plt.ylim(0, 22050)
    plt.xlabel('Samples')
    plt.ylabel('Amplitude')
    plt.xlim(0, len(ballon))

fileToOpen = "Data/Direkte_lyd/Sinus/FindRumImpulsSINUS.txt"
Sweep = np.loadtxt("Data/Direkte_lyd/Sinus/SweepSinus")
dataBlocks = LoadLogFile(fileToOpen)
dataBlock = dataBlocks[1]
time = dataBlock[0].replace("time: ",'')
Optagelse = np.loadtxt("Data/Direkte_lyd/Sinus/Optagelse-" + time)
y = Optagelse[:,0]
xpc = Optagelse[:,1]
Sweep = Sweep[:len(y)]

PlotRoom(mode='tid', HOJ="Speaker0ClipUdenDelay.txt")

##latexify()
##plt.legend()
plt.show()

