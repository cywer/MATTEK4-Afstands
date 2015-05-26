# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 12:31:59 2015

@author: root
"""  

from __future__ import division
from matplotlib import pyplot as plt
import numpy as np
from Funktions_Signal import KrydsKoreller,ImprovedGaussian,PhaseLinReg,FindReflektion, FindImpuls, Udfold, latexify, Punkt2Afstand, Afstand2Vag, Afstand2Punkt, Interpolation

"""
Samtlige filer med måledata fra præcisionsforsøget åbnes.
"""

fileToOpen =  ['Data/PraecisionsForsoeg/PF1/FindRumImpulsPF1.txt',
               'Data/PraecisionsForsoeg/PF2/FindRumImpulsPF2Papir.txt',
               'Data/PraecisionsForsoeg/PF3/FindRumImpulsPF3Haefte.txt']

def LoadLogFile(fileToOpen) : 
    dataBlocks = []
    dataBlock = []
    with open(fileToOpen) as fp:
        for line in fp:
            if line.startswith('\r\n') or line.startswith('\n') : 
                dataBlocks.append(dataBlock) 
                dataBlock = []
            else : 
                dataBlock.append(line.replace('\r\n', ''))
    return dataBlocks         

SweepSti = 'Data/SweepSinus'
Sweep = np.loadtxt(SweepSti)

"""
Afstande beregnes med interpolationsmetoden, forbedret gauss-\
metode og fasemetoden
"""

afstande = []
afstandeFase = []
afstandeGauss = []
for i0,i in enumerate(fileToOpen):
    log = LoadLogFile(i)
    Opt = 'Data/PraecisionsForsoeg/PF' + str(i0 +1) +'/Optagelse-'
    for k in log:
        OptagelseSti = Opt + k[0].split()[1]
        Optagelse = np.loadtxt(OptagelseSti)
        y = Optagelse[:,0]
        x_pc = Optagelse[:,1]
#        h = Udfold(y,x_pc)
        h = KrydsKoreller(y,x_pc)
        gl, diff = FindReflektion(h)
        
        ny = Interpolation(h,gl)
        nyrefl = Interpolation(h,gl+diff)
        mm = Punkt2Afstand(nyrefl-ny)/2 * 1000
        afstande.append(mm)
        
        ny = ImprovedGaussian(h,gl,1000)
        nyrefl = ImprovedGaussian(h,gl+diff,1000)
        mm = Punkt2Afstand(nyrefl-ny)/2 * 1000
        afstandeGauss.append(mm)
        
        ny = PhaseLinReg(np.fft.fft(x_pc),np.fft.fft(y),gl)
        nyrefl = PhaseLinReg(np.fft.fft(x_pc),np.fft.fft(y),gl+diff)
        mm = Punkt2Afstand(nyrefl-ny)/2 * 1000
        afstandeFase.append(mm)
    

"""
Figurer oprettes og gemmes til brug i rapport
"""

latexify(fig_width=7,fig_height=1.5)
colors = ['r','g','b']
tegn = ['ro','g>','b<']
gsA = []

plt.figure(0)
plt.xlim(495,500)
plt.xlabel('Millimeter')
for c,i in enumerate([afstande[:10],afstande[10:20],afstande[20:30]]):
    for k in i:
        plt.plot(k,0,tegn[c])
    gs = np.mean(i)
    gsA = np.append(gsA,gs)
    plt.axvline(gs,color=colors[c])
plt.savefig('Interpolation.pdf')

plt.figure(1)
plt.xlim(495,500)
plt.xlabel('Millimeter')
for c,i in enumerate([afstandeGauss[:10],afstandeGauss[10:20],afstandeGauss[20:30]]):
    for k in i:
        plt.plot(k,0,tegn[c])
    gs = np.mean(i)
    gsA = np.append(gsA,gs)
    plt.axvline(gs,color=colors[c])
plt.savefig('GaussMetoden.pdf')
    
plt.figure(2)
plt.xlabel('Millimeter')
for c,i in enumerate([afstandeFase[:10],afstandeFase[10:20],afstandeFase[20:30]]):
    for k in i:
        plt.plot(k,0,tegn[c])
    gs = np.mean(i)
    gsA = np.append(gsA,gs)
    plt.axvline(gs,color=colors[c])
plt.savefig('FaseMetoden.pdf')

plt.show()  

"""
Forskellen i afstand ved de forskellige metoder printes.
"""

print "Interpolation:"
print "%.2f , %.2f" %(gsA[0]-gsA[1] , gsA[0]-gsA[2])
print "ImprovedGaussian:"
print "%.2f , %.2f" %(gsA[3]-gsA[4] , gsA[3]-gsA[5])
print "FaseMetoden:"
print "%.2f , %.2f" %(gsA[6]-gsA[7] , gsA[6]-gsA[8])