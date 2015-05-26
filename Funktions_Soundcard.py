# -*- coding: utf-8 -*-
"""
Created on Mon Mar  9 09:56:27 2015

@author: cywer
"""
#http://docs.scipy.org/doc/scipy/reference/tutorial/signal.html
#https://github.com/bastibe/PySoundFile
#https://github.com/bastibe/PySoundCard

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import pysoundcard 
from pysoundcard import Stream, continue_flag
import time

def pysoundcard_test():
    """
    Retunerer microfon og højtaler dictionarys 
    """
    return pysoundcard.default_input_device(),pysoundcard.default_output_device()

def Afspil(Lyd):
     """
    Afspil en lyd med 
    
    Parameters
    -------
    Lyd : 1D ndarray of floats
        Digital lyd gemt i et array. 
    """
    s = Stream()
    s.start()
    s.write(Lyd)
    s.stop()
    
def Optag(Tid):
    """
    Optag lyd fra mikrofonen, og gem det i array.  
    
    Parameters
    -------
    Tid : float
        Tid over hvor længe der skal optages.
    Returns
    -------
    Optagelse : 1D ndarray of floats
        Digital lyd gemt i array. 
    """
    s = Stream()
    s.start()
    Optagelse = s.read(44100*Tid)
    s.stop()
    return Optagelse
    
if __name__ == '__main__':
    Optagelse = Optag(3)
    OptagelseSti = "-".join(('Data/Optagelse','1'))
    np.savetxt(OptagelseSti,Optagelse)  
    plt.plot(Optagelse)
