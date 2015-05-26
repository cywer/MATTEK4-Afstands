# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 22:41:56 2015

@author: Cywer
"""
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
import Funktions_Signal as fs
import Sweeps as sweeps

 
def Save_Testdata():
    string='Data/Tids_Simulering/tidsforsinkelse_'
    np.savetxt(string+'test_f_'+str(SNR_dB)+'_dB.csv',test_f)
    np.savetxt(string+'test_k_'+str(SNR_dB)+'_dB.csv',test_k)

def print_test(test):
    print 'Delay, argmax, inter, PhLin'
    for i in np.vsplit(test,len(test)):
        print ("{:6.2f} {:6.2f} {:6.2f} {:6.2f} ".format(
                i[0][0],i[0][1],i[0][2],i[0][5]))
    
def test_sum(test):
    print 'Metode, Genemsnits fejl, LeastSqures'
    metodes=['Argmax','Interp','Cos','ImGaus','Phase']#,'PhaseN']
    for i in range(1,len(metodes)+1):
        print ('{:6} {:f} {:6.1f} {:f} {:10.1f} {:f}'.format(
        metodes[i-1]
        ,np.mean(np.abs(test[:,0]-test[:,i]))
        ,1/(np.mean(np.abs(test[:,0]-test[:,i]))/np.mean(np.abs(test[:,0]-test[:,1])))
        ,np.mean((test[:,0]-test[:,i])**2)
        ,1/(np.mean((test[:,0]-test[:,i])**2)/np.mean((test[:,0]-test[:,1])**2) )
        ,np.std((test[:,0]-test[:,i]))
        ))
    
def test_sum_latex():
    
    def Load_Testdata(load_dB,mode=''):
        string='Data/Tids_Simulering/tidsforsinkelse_'
        test = np.loadtxt(string+mode+str(load_dB)+'_dB.csv')
        print 'Loading test results'
        return test
        
    metodes=['Argmax','Interpolation','Cosinus fitting','Improved gaussian','Fasemetoden']#,'PhaseN']
    tabel=np.zeros((len(metodes),len(dB))) 
    for i1,v in enumerate(dB):
        test=Load_Testdata(v,mode='test_f_')
        for i in range(1,len(metodes)+1):
#            tabel[i-1,i1]=np.mean(np.abs(test[:,0]-test[:,i])**2)
            tabel[i-1,i1]=1/(np.mean(np.abs(test[:,0]-test[:,i])**2)/np.mean(np.abs(test[:,0]-test[:,1])**2))
    
    for i in range(1,len(metodes)+1):
        print ('\\hline')
        print ('{:6} & {:6.0f} & {:6.0f} & {:6.0f} & {:6.0f} & {:6.0f} & {:6.0f} \\\\'.format(
        metodes[i-1]
        ,tabel[i-1,0]
        ,tabel[i-1,1]
        ,tabel[i-1,2]
        ,tabel[i-1,3]
        ,tabel[i-1,4]
        ,tabel[i-1,5]
        ))

# Genererer uforsinket sweep
Sweep_x,SampRate = sweeps.SinusSweepSubsample(1,44100//2,1,Delay=0)
Sweep_x = np.hstack([np.zeros(44100),Sweep_x,np.zeros(44100*1)])
Sweep_X = np.fft.fft(Sweep_x)

n_test=100
test_range = np.linspace(234.,235.,n_test,endpoint=False)

dB=['inf',40,30,20,10,0] # De ønskede mængder af SNR
SNR_dB=dB[1]


test_f=np.zeros((n_test,6))
test_k=test_f.copy()
for i,v in enumerate(test_range):
    # Forsinket sweep generes
    Sweep_y0,SampRate = sweeps.SinusSweepSubsample(1,44100//2,1,Delay=v-int(v))
    Sweep_y = np.hstack([np.zeros(44100+int(v)),Sweep_y0,np.zeros(44100-int(v))])

    # Støj tilføres
    if SNR_dB is not 'inf':
        stoj=fs.White_noise(Sweep_y,SNR_dB)
        Sweep_y=Sweep_y+stoj
    Sweep_Y = np.fft.fft(Sweep_y)
    
    # Impulsresponcet findes
    f = np.real(np.fft.ifft(Sweep_Y/Sweep_X))

    #Find Indekset til den direkte impulse
    argmaxf,_ = fs.FindReflektion(f)   
    
    # Gemmer test resultater for affoldning
    test_f[i,:] =(v, argmaxf
                ,fs.Interpolation(f,argmaxf)
                ,fs.Cosfitting(f,argmaxf)
                ,fs.ImprovedGaussian(f,argmaxf,0.1)
                ,fs.PhaseLinReg(Sweep_X,Sweep_Y,argmaxf)
#                ,PhaseDelayImpulse(f,argmaxf)
                )
       
    # Correlationen findes    
    k = np.real(np.fft.ifft(Sweep_Y*np.conj(Sweep_X)))

    #Find Indekset til den direkte impulse
    argmaxk,_ = fs.FindReflektion(k)   

    test_k[i,:] =(v, argmaxk
            ,fs.Interpolation(k,argmaxk)
            ,fs.Cosfitting(k,argmaxk)
            ,fs.ImprovedGaussian(k,argmaxk,1000)
            ,fs.PhaseLinReg(Sweep_X,Sweep_Y,argmaxk)
#            ,PhaseDelay(Sweep_X,Sweep_Y,argmaxk)
            )   




#print_test(test_f)
#print
#print_test(test_k)
#print

print 'SNR_db',SNR_dB
print 'Test Foldning'
test_sum(test_f)
print 
print 'Test Krydskorrelation'
test_sum(test_k)

#test_sum_latex()