# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 10:31:49 2015

@author: cywer
"""
import time
import subprocess

"""
Gentager den samme måling N gange.

Køres i terminal. Kalder FindRumImpuls.py N gange, således at et sweep 
afspilles, og rummets svar optages. Hvert resultat gemmes desuden i en ny
datafil. 
"""


#time.sleep(5)

N=10 # Lav forsøget N gange
for i in range(N):
        subprocess.call([['python'],['FindRumImpuls.py']], shell=True)
