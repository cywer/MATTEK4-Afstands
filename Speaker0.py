# -*- coding: utf-8 -*-
"""
Created on Tue Mar 03 09:35:09 2015

@author: Nneling
"""
from __future__ import division
import matplotlib.pyplot as plt
import numpy as np
from scipy import signal
import Funktions_Signal as fs

def LoadLogFile(fileToOpen) : 

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

fileToOpen = 'Data/360/FindRumImpulsPolar360.txt'
dataBlocks = LoadLogFile(fileToOpen)

Speaker0F = None
for x in dataBlocks[:3]:
    time = x[0].replace("time: ",'')
    Optagelse = np.loadtxt("Data/360/Optagelse-" + time)
    xpc = Optagelse[:,1]
    y = Optagelse[:,0]
    Xpc = np.fft.fft(xpc)
    Y = np.fft.fft(y)
    ImportpulsFrekvens = Y/Xpc
    if Speaker0F == None:
        Speaker0F = np.zeros_like(ImportpulsFrekvens)
    Speaker0F += ImportpulsFrekvens
Speaker0F /= 3
freq = np.fft.fftfreq(len(y),1/44100)
Speaker0 = np.fft.ifft(Speaker0F)
hann = signal.hann(2*(260 - np.argmax(Speaker0)) + 1)
diff = np.argmax(Speaker0) - np.argmax(hann)
shiftHann = np.concatenate([np.zeros(diff), hann, np.zeros(len(Speaker0) - diff - len(hann))], axis=0)
print len(shiftHann), len(Speaker0)


def Clip() : 
    plt.plot(shiftHann * np.max(Speaker0))
    plt.xlim(100, 300)
    plt.xlabel('Samples')
    plt.ylabel('Amplitude')
    plt.plot(Speaker0)

##Clip()

clipSpeaker0 = Speaker0 * shiftHann
clipSpeaker0F = np.fft.fft(clipSpeaker0)

Speaker0F = np.fft.fft(Speaker0)
Delay = np.zeros_like(Speaker0)
Delay[np.argmax(Speaker0)] = 1
DelayF = np.fft.fft(Delay)
clipSpeaker0F = clipSpeaker0F/DelayF
hoejtalerFrekvens = clipSpeaker0F

plt.plot(np.abs(clipSpeaker0F))

#remove preringing
##b, a = signal.butter(1, 15000 / 22050, 'low') ## 15000 er hvor knk frekvens er. 
##hfilt = signal.filtfilt(b, a, clipSpeaker0)
##hfilt = np.concatenate([hfilt[192:], np.zeros(192)], axis=0)
##Hfilt = np.fft.rfft(hfilt)
##plt.plot(np.fft.rfftfreq(len(hfilt),1/44100), np.abs(Hfilt), 'r')
##plt.plot(hfilt, 'r')

np.savetxt('Speaker0ClipUdenDelay.txt', clipSpeaker0F.view(float))

plt.show()
