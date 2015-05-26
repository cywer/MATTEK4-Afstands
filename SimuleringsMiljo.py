# -*- coding: utf-8 -*-
"""
Created on Mon Mar 30 09:33:29 2015

@author: root
"""

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal

from Funktions_Signal import FindImpuls,latexify
from SimuleringsFunktioner import FindOutputNyImpuls,LavImpuls,Forsinkelse,Reflektion,Pad,Plotsweep,PlotTuples,mls,SinusSweep

PadNuller=10000
Gain = 0.5
Enheder = 5000
MLSInit = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,0]

"""
For alle funktioner:
  Der er 3 mulige modes: Delay,Reflektion og Rigtig.
  Padding kan vælges som True eller False, tilføjer nuller i begge
    sider af signalet, hvis True.
  PlotFrekvens kan vælges som True eller False.
    False giver et plot af input, output, impuls og nyimpuls i tidsdomænet.
    True giver et plot af input, output, impuls og nyimpuls i frekvensdomænet.

Der er følgende funktioner tilgængelige:
  Standard logaritmisk sinus sweep       Sin()
  Standard MLS-signal                    MLS()
  Padded sinus sweep                     SinPadded()
  Padded MLS-signal                      MLSPadded()
  Forkortet sinus sweep                  SinForkortet()
  Forkortet MLS-signal                   MLSForkortet()
  Forlænget sinus sweep                  SinForlaenget()
  Forlænget MLS-signal                   MLSForlaenget()
  Lave frekvenser sinus sweep            SinLavfrek()
  Høje frekvenser sinus sweep            SinHojfrek()
  For høje frekvenser sinus sweep        SinForHojefrek()
  Middelfrekvenser sinus sweep           SinMiddelfrek()
  MLS hvor nogle punkter mangler         MLSMangler()
  Sinus sweep med støj                   SinStoj()
  MLS-signal med støj                    MLSStoj()
"""

def Sin(Mode='Delay',Padding=False,PlotFrekvens=False):
    """Standard logaritmisk sinus sweep"""
    Sweep,SampF = SinusSweep(20,20000,1)
    Sweep = Pad(Sweep,PadNuller,Padding)
    Impuls = LavImpuls(Mode,Sweep)
    Output, NyImpuls = FindOutputNyImpuls(Sweep,Impuls)
    PlotTuples(Sweep,Output,Impuls,NyImpuls,PlotFrekvens)
    
def MLS(Mode='Delay',Padding=False,PlotFrekvens=False):
    """Standard MLS-signal"""
    Sweep = mls(2**15 - 2,MLSInit)
    Sweep = Pad(Sweep,PadNuller,Padding)
    Impuls = LavImpuls(Mode,Sweep)
    Output, NyImpuls = FindOutputNyImpuls(Sweep,Impuls) 
    PlotTuples(Sweep,Output,Impuls,NyImpuls,PlotFrekvens)
    
def SinPadded(Mode='Delay',Padding=True,PlotFrekvens=False):
    """ LogSinus, padded som standard"""
    Sweep,SampF = SinusSweep(20,20000,1)
    Sweep = Pad(Sweep,PadNuller,Padding)
    Impuls = LavImpuls(Mode,Sweep)
    Output, NyImpuls = FindOutputNyImpuls(Sweep,Impuls)
    PlotTuples(Sweep,Output,Impuls,NyImpuls,PlotFrekvens)
    ()
def MLSPadded(Mode='Delay',Padding=True,PlotFrekvens=False):
    """MLS, padded som standard"""
    Sweep = mls(2**15 - 2,MLSInit)
    Sweep = Pad(Sweep,PadNuller,Padding)
    Impuls = LavImpuls(Mode,Sweep)
    Output, NyImpuls = FindOutputNyImpuls(Sweep,Impuls)
    PlotTuples(Sweep,Output,Impuls,NyImpuls,PlotFrekvens)
    
def SinForkortet(Mode='Delay',Padding=False,Tid=0.2,PlotFrekvens=False):
    """LogSin, tid kan varieres. 0.2 sekunder standard"""
    Sweep,SampF = SinusSweep(20,20000,Tid)
    Sweep = Pad(Sweep,PadNuller,Padding)
    Impuls = LavImpuls(Mode,Sweep)
    Output, NyImpuls = FindOutputNyImpuls(Sweep,Impuls)
    PlotTuples(Sweep,Output,Impuls,NyImpuls,PlotFrekvens)
    
def MLSForkortet(Mode='Delay',Padding=True,PlotFrekvens=False,Laengde=2**7 -2):
    """MLS, længde kan varieres. 2**7 -2 standard. Skal være padded"""    
    Sweep = mls(Laengde,MLSInit)
    Sweep = Pad(Sweep,PadNuller,True) #Skal være padded
    Impuls = LavImpuls(Mode,Sweep)
    Output, NyImpuls = FindOutputNyImpuls(Sweep,Impuls)
    PlotTuples(Sweep,Output,Impuls,NyImpuls,PlotFrekvens)

def SinForlaenget(Mode='Delay',Padding=False,Tid=5,PlotFrekvens=False):
    """LogSin, tid kan varieres. 5 sekunder standard"""
    Sweep,SampF = SinusSweep(20,20000,Tid)
    Sweep = Pad(Sweep,PadNuller,Padding)
    Impuls = LavImpuls(Mode,Sweep)
    Output, NyImpuls = FindOutputNyImpuls(Sweep,Impuls)
    PlotTuples(Sweep,Output,Impuls,NyImpuls,PlotFrekvens)
    
def MLSForlaenget(Mode='Delay',Padding=False,Gange=5,PlotFrekvens=False):
    """MLS, signalet forlænges X gange. 5 er standard"""
    Sweep = mls(2**15 - 2,MLSInit)
    Sweep = np.hstack([(Sweep)]*Gange)
    Sweep = Pad(Sweep,PadNuller,Padding)
    Impuls = LavImpuls(Mode,Sweep)
    Output, NyImpuls = FindOutputNyImpuls(Sweep,Impuls)
    PlotTuples(Sweep,Output,Impuls,NyImpuls,PlotFrekvens)

def SinLavfrek(Mode='Delay',Padding=False,StopFrek=20,PlotFrekvens=False):
    """LogSin, stopfrekvensen kan varieres. Standard 500 Hz"""
    Sweep,SampF = SinusSweep(1,StopFrek,1)
    Sweep = Pad(Sweep,PadNuller,Padding)
    Impuls = LavImpuls(Mode,Sweep)
    Output, NyImpuls = FindOutputNyImpuls(Sweep,Impuls)
    PlotTuples(Sweep,Output,Impuls,NyImpuls,PlotFrekvens)
    
def SinHojfrek(Mode='Delay',Padding=False,StartFrek=15000,PlotFrekvens=False):
    """LogSin, startfrekvensen kan varieres. Standard 15000 Hz"""    
    Sweep,SampF = SinusSweep(StartFrek,20000,1)
    Sweep = Pad(Sweep,PadNuller,Padding)
    Impuls = LavImpuls(Mode,Sweep)
    Output, NyImpuls = FindOutputNyImpuls(Sweep,Impuls)
    PlotTuples(Sweep,Output,Impuls,NyImpuls,PlotFrekvens)
    
def SinForHojefrek(Mode='Delay',Padding=False,StartFrek=35000,SlutFrek=40000,PlotFrekvens=False):
    """LogSin, frekvensen over 1/2 startfrekvens. Start og slutfrekvensen 
    kan varieres. Standard er 35000 og 40000"""    
    Sweep,SampF = SinusSweep(StartFrek,SlutFrek,1)
    Sweep = Pad(Sweep,PadNuller,Padding)
    Impuls = LavImpuls(Mode,Sweep)
    Output, NyImpuls = FindOutputNyImpuls(Sweep,Impuls)
    PlotTuples(Sweep,Output,Impuls,NyImpuls,PlotFrekvens)
    
def SinMiddelfrek(Mode='Delay',Padding=False,StartFrek=7500,StopFrek=12500,PlotFrekvens=False):
    """LogSin, frekvenser midt i det hørbare. Start-og stopfrekvens kan 
    varieres. Standard er 7500 og 12500"""    
    Sweep,SampF = SinusSweep(StartFrek,StopFrek,1)
    Sweep = Pad(Sweep,PadNuller,Padding)
    Impuls = LavImpuls(Mode,Sweep)
    Output, NyImpuls = FindOutputNyImpuls(Sweep,Impuls)
    PlotTuples(Sweep,Output,Impuls,NyImpuls,PlotFrekvens)
    
def MLSMangler(Mode='Delay',Padding=False,Mangel=0.25,PlotFrekvens=False):
    """MLS, hvor noget af signalet "mangler"(sættes til 0). Standard er 0.25"""
    Sweep = mls(2**15 - 2,MLSInit)
    Sweep = np.hstack((Sweep[:len(Sweep)*(1-Mangel)],np.zeros(len(Sweep)*Mangel)))
    Sweep = Pad(Sweep,PadNuller,Padding)
    Impuls = LavImpuls(Mode,Sweep)
    Output, NyImpuls = FindOutputNyImpuls(Sweep,Impuls)
    PlotTuples(Sweep,Output,Impuls,NyImpuls,PlotFrekvens)
    
def SinStoj(Mode='Delay',Padding=False,Spredning=0.1,PlotFrekvens=False):
    """LogSin, hvor der tilføjes normalfordelt støj til outputtet.
    Støjens spredning kan varieres, standard er 0.1"""    
    Sweep,SampF = SinusSweep(20,20000,1)
    Sweep = Pad(Sweep,PadNuller,Padding)
    Impuls = LavImpuls(Mode,Sweep)
    Output = signal.fftconvolve(Sweep,Impuls,mode='full')
    Output = Output + np.random.normal(0,Spredning,len(Output))
    NyImpuls = FindImpuls(Sweep,Output)    
    PlotTuples(Sweep,Output,Impuls,NyImpuls,PlotFrekvens)

def MLSStoj(Mode='Delay',Padding=False,Spredning=0.1,PlotFrekvens=False):
    """MLS, hvor der tilføjes normalfordelt støj til outputtet.
    Støjens spredning kan varieres, standard er 0.1"""  
    Sweep = mls(2**15 - 2,MLSInit)
    Sweep = Pad(Sweep,PadNuller,Padding)
    Impuls = LavImpuls(Mode,Sweep)
    Output = signal.fftconvolve(Sweep,Impuls,mode='full')
    Output = Output + np.random.normal(0,Spredning,len(Output))
    NyImpuls = FindImpuls(Sweep,Output)    
    PlotTuples(Sweep,Output,Impuls,NyImpuls,PlotFrekvens)

latexify()