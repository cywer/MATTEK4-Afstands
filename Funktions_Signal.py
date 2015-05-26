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
from scipy import signal
import Sweeps as sweeps
"""
TO DO
- Funktion forklaringer
- Funktions test ( reproducerbart gør koden de den skal)
"""

def FindImpuls(Ind,Ud):
    """
    Beregner Impulsrespons, ved at affolde output (Ud) med input (Ind) via DTF.
    
    Parameters
    -------
    Ind : 1D ndarray of floats
        Input signal i tidsdomæne.
    Ud : 1D ndarray of floats
        Output signal i tidsdomæne.
    Returns
    -------
    output_signal : 1D ndarray of floats
        Det affoldede signal i tidsdomænet.
    """
    N = np.shape(Ind)[0]
    IndFrek = np.fft.fft(Ind[:N])
    UdFrek = np.fft.fft(Ud[:N])
    ImpulsBFrek = UdFrek/IndFrek
    ImpulsB = np.fft.ifft(ImpulsBFrek)
    return ImpulsB
    
def Udfold(y,xpc,HOJ='Speaker0ClipUdenDelay.txt',mode=None):       
    """
    Beregner rummets impuls ved at affolde output (y) med input (xpc) og
    højtaleren (HOJ).

    Parameters
    -------
    y : 1D ndarray of floats
        Output signal i tidsdomæne.
    xpc : 1D ndarray of floats
        Signalet xpc er computerens forsinkelse tidsdomæne.
    HOJ : string
        Filnavnet på højtaleren frekvensrespons, der skal affoldes med.
    Returns
    -------
    rum : 1D ndarray of floats
        Rummets impulsrespons i tidsdomænet.
    """
    Y=np.fft.fft(y)
    Xpc=np.fft.fft(xpc)
    H0=np.loadtxt(HOJ).view(complex)
    
#    plt.plot(np.fft.ifft(H0))
    
    freq=np.fft.fftfreq(len(y),1./44100.)
    
    RUM=Y/Xpc/H0
    RUM[np.where(np.abs(freq)>20000)]=0
    RUM[np.where(np.abs(freq)<20)]=0
    
    rum=np.real(np.fft.ifft(RUM))
    
    if mode == 'full':
        return rum , freq , H0
    return rum

def KrydsKoreller(y,xpc,HOJ='Speaker0ClipUdenDelay.txt',mode=None):       
    """
    Beregner rummets impuls ved at affolde output (y) med input (xpc) og
    højtaleren (HOJ).

    Parameters
    -------
    y : 1D ndarray of floats
        Output signal i tidsdomæne.
    xpc : 1D ndarray of floats
        Signalet xpc er computerens forsinkelse tidsdomæne.
    HOJ : string
        Filnavnet på højtaleren frekvensrespons, der skal affoldes med.
    Returns
    -------
    rum : 1D ndarray of floats
        Rummets impulsrespons i tidsdomænet.
    """
    Y=np.fft.fft(y)
    Xpc=np.conj(np.fft.fft(xpc))
    H0=np.loadtxt(HOJ).view(complex)
    freq=np.fft.fftfreq(len(y),1./44100.)
    
    RUM=Y*Xpc
    RUM[np.where(np.abs(freq)>20000)]=0
    RUM[np.where(np.abs(freq)<20)]=0
    
    rum=np.real(np.fft.ifft(RUM))
    
    if mode == 'full':
        return rum , freq , H0
    return rum
    
def dB(array):
    """
    Addditive White Gaussian Noise (AWGN) Channel.
    
    Parameters
    -------
    input_signal : 1D ndarray of floats
        Input signal to the channel.
    snr_dB : float
        Output SNR required in dB.
    Returns
    -------
    output_signal : 1D ndarray of floats
        Output signal from the channel with the specified SNR.
    """
    return 20*np.log10(array)

def SNR(signal, stoej, mode='') :
    """
    Beregner signalstøjforholdet med root-mean-square metoden. 
    OBS. der er forskel på om sammenligne signal/støj og 
    signal+støj / støj, især ved lave signalstøjforhold.
    
    Parameters
    -------
    input_signal : 1D ndarray of floats
        Array med signalet
    stoej : 1D ndarray of floats
        Array kun med støj
    Returns
    -------
    SNR : float
        Signalstøjforholdet
    """
    rmsSignal = np.sqrt(np.mean(signal **2))
    rmsNoise = np.sqrt(np.mean(stoej **2))
    if mode == 'db' :
        SNRdb = 20 * np.log10(rmsSignal / rmsNoise)
        return SNRdb
    elif mode == '':
        SNR = (rmsSignal / rmsNoise)**2
        return SNR
    else:
        raise 'Invalid mode'

def White_noise(input_signal, snr_dB):
    """
    Generate appropiate White Gaussian Noise.
    
    Parameters
    -------
    input_signal : 1D ndarray of floats
        Pure input signal.
    snr_dB : float
        Final output SNR required in dB.
    Returns
    -------
    noise : 1D ndarray of floats
        Noise to be added to the input_signal, 
        which results in the specified SNR.
    """

    NonzeroSignal = input_signal[np.nonzero(input_signal)]
    avg_energy = np.sqrt(np.mean(NonzeroSignal**2))
    
    snr_linear = 10**(snr_dB/10.0)
    noise_std = np.sqrt(avg_energy**2/(snr_linear))

    if input_signal.dtype is complex:
#        noise = np.sqrt(noise_variance) * np.random.randn(len(input_signal)) * (1+1j)
        raise 'Not implemented'
    else:
        noise = noise_std * np.random.randn(len(input_signal))

    return noise

def Punkt2Afstand(Punkt,freq=44100.,vel=343.):
    """
    Funktion finder Afstanden mellem et givet 
    antal målepunkter.
    
    Parameters
    -------
    Punkt : float
        Et antal punkter.
    freq : float
        Samplings frekvensen.
    vel : float
        Lydens hastighed i meter per sekund.
    Returns
    -------
    afstand : float
        Afstanden i meter.
    """
    return Punkt * vel/freq
    
def Afstand2Punkt(afstand,freq=44100,vel=343.):
    """
    Funktionen finder et antal målepunkter mellem en gvet afstand.
    
    Parameters
    -------
    afstand : float
        En afstand i meter.
    freq : float
        Samplings frekvensen.
    vel : float
        Lydens hastighed i meter per sekund.
    Returns
    -------
    punkter : float
        Et antal punkter imellem afstand.
    """
    return afstand * freq / vel 

def Afstand2Vag(r, dm):
    """
    Funktionen finder afstanden mellem et målepunkt og væg
    udfra en refleksions måling og afstand mellem måleudstyret.
    Se afsnit 'Relativ afstand'
    
    Parameters
    -------
    r : float
        Afstanden mellem højtaleren akustisk center og mikrofonen i meter.
    dm : float
        Antal målingspunter mellem den direktelyd og refelktionen.
    Returns
    -------
    afstand : float
        Afstanden til vægen
    """
    d = Afstand2Punkt(dm)
    dv = np.sqrt(r**2-(d/2)**2)
    return Punkt2Afstand(dv)

    
def FindReflektion(h,endp=None):
    """
    Finder sample nummer til reflektionen af en lyd, som ligger over 7 samples
    væk fra den direkte lyd. 

    Parameters
    -------
    h : 1D ndarray of floats
        Impulsresponsen for et rum.
    Returns
    -------
    sample : int
        Sample af reflektionen. 
    """
    gl = np.argmax(h)
    diff = np.argmax(h[gl+7:10000]) + 7
    return gl, diff

def Interpolation(array,argmax):
    """
    Bruger interpolation til at forøge præcisionen af afstandsbestemmelse.
    For formel, se afsnit Interpolation i rapporten. 

    Parameters
    -------
    array : 1D ndarray of floats
        Signalet estimationen skal laves på.
    argmax : float
        Indekset på samplet som ønskes estimeret.
    Returns
    -------
    subsample : float
        Subsample estimationen af det ønkede delay.
    """
    x_1,x0,x1 = array[argmax-1], array[argmax],array[argmax+1]
    t = x_1 - x1
    n=2*(x_1 - 2*x0 + x1)
    return argmax + t/n

def Cosfitting(array,argmax):
    """
    Bruger interpolation til at forøge præcisionen af afstandsbestemmelse.
    For formel, se (Wiens og Bradley 2009). 

    Parameters
    -------
    array : 1D ndarray of floats
        Signalet estimationen skal laves på.
    argmax : float
        Indekset på samplet som ønskes estimeret.
    Returns
    -------
    subsample : float
        Subsample estimationen af det ønkede delay.
    """
    w=np.arccos((array[argmax-1]+array[argmax+1])/(2*array[argmax]))
    theta=np.arctan((array[argmax-1]-array[argmax+1])/(2*array[argmax]*np.sin(w)))
    return argmax - theta/w

def ImprovedGaussian(array,argmax,alpha):
    """
    Bruger interpolation til at forøge præcisionen af afstandsbestemmelse.
    For formel, se (Wiens og Bradley 2009). 

    Parameters
    -------
    array : 1D ndarray of floats
        Signalet estimationen skal laves på.
    argmax : float
        Indekset på samplet som ønskes estimeret.
    alpha : float
        Konstant til at skrue på.
    Returns
    -------
    subsample : float
        Subsample estimationen af det ønkede delay.
    """
    b = -np.min((array[argmax-1], array[argmax+1])) + alpha
    x_1,x0,x1 = np.log(np.array([array[argmax-1], array[argmax],array[argmax+1]])+b)

    t = x_1 - x1
    n=2*(x_1 - 2*x0 + x1)
    t_d = argmax + t/n
    return t_d


def PhaseLinReg(X0,X1,argmax):
    """
    Bruger en lineær approksimation af fasen fra frekvensresponset(DTF) 
    til at bestemme et signal. For formel, se (Wiens og Bradley 2009). 

    Parameters
    -------
    X0 : 1D ndarray of floats
        Det uforsinkede signal i frekvensdomænet.
    X1 : 1D ndarray of floats
        Det forsinkede signal i frekvensdomænet.        
    argmax : float
        Indekset på den nærmeste int af forsninkelsen som ønskes estimeret.
    Returns
    -------
    subsample : float
        Subsample estimationen af det ønkede delay.
    """
    Phase = X1 *np.exp(-2j*np.pi* np.arange(len(X1))*-argmax /len(X1))
    ETFE  = Phase/X0

    x = np.fft.fftfreq(len(X0),0.5/np.pi)
    x=np.fft.fftshift(x)    
    
    ang =np.angle(ETFE)
    ang = np.fft.fftshift(ang)
    
    freq_lim = np.pi*0.9
    index=np.where(np.abs(x)<freq_lim)
    
    poly = np.polyfit(x[index],ang[index],1)    
    
    return argmax - poly[0]


#==============================================================================
# Plotting funktions
#==============================================================================

def get_filenames(path='./Data/'):    
    """
    Get at list of files in folder
    
    Parameters
    -------
    path : string
        Relative path to folder. An Example: './Data/360/'
    Returns
    -------
    filenames : list of strings
        List af filenames in the folder.
    """
    from os import walk
    
    for (dirpath, dirnames, filenames) in walk(path):
        break
    return filenames


def latexify(fig_scale=1,fig_width=None, fig_height=None, columns=2):
    import matplotlib
    SPINE_COLOR = 'gray'
    """Set up matplotlib's RC params for LaTeX plotting.
    Call this before plotting a figure.

    Parameters
    ----------
    fig_width : float, optional, inches
    fig_height : float,  optional, inches
    columns : {1, 2}
    """

    # code adapted from http://www.scipy.org/Cookbook/Matplotlib/LaTeX_Examples

    # Width and max height in inches for IEEE journals taken from
    # computer.org/cms/Computer.org/Journal%20templates/transactions_art_guide.pdf

    assert(columns in [1,2])

    if fig_width is None:
        fig_width = 3.39*fig_scale if columns==1 else 6*fig_scale # width in inches

    if fig_height is None:
        golden_mean = 0.6180339887498949#(sqrt(5)-1.0)/2.0# Aesthetic ratio
        fig_height = fig_width*golden_mean # height in inches

    MAX_HEIGHT_INCHES = 8.0
    if fig_height > MAX_HEIGHT_INCHES:
        print("WARNING: fig_height too large:" + fig_height + 
              "so will reduce to" + MAX_HEIGHT_INCHES + "inches.")
        fig_height = MAX_HEIGHT_INCHES

    params = {'backend': 'ps',
              'text.latex.preamble': ['\usepackage{gensymb}'],
              'axes.labelsize': 8, # fontsize for x and y labels (was 10)
              'axes.titlesize': 8,
              'text.fontsize': 8, # was 10
              'legend.fontsize': 8, # was 10
              'xtick.labelsize': 8,
              'ytick.labelsize': 8,
              'text.usetex': True,
              'figure.figsize': [fig_width,fig_height],
              'font.family': 'serif'
    }

    matplotlib.rcParams.update(params)

def format_axes(ax):

    for spine in ['top', 'right']:
        ax.spines[spine].set_visible(False)

    for spine in ['left', 'bottom']:
        ax.spines[spine].set_color(SPINE_COLOR)
        ax.spines[spine].set_linewidth(0.5)

    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')

    for axis in [ax.xaxis, ax.yaxis]:
        axis.set_tick_params(direction='out', color=SPINE_COLOR)

    return ax

def plotsweep_test(d1,d2):
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
    
    ax4 = plt.subplot(gridspec2[gsx//2:,:],sharex=ax3)
    plt.plot(d2[0], d2[2])

    fig.suptitle('An overall title', size=20)
    plt.setp(ax1, title='Test')
    plt.setp(ax2, title='Test')
    plt.setp(ax3, title='Test')
    plt.setp(ax4, title='Test')

    
def PlotSweep(Sweep,Varighed):
    latexify(fig_scale=1)
    Tid = np.linspace(0,Varighed,len(Sweep))
    fig, ax = plt.subplots(dpi=200)
    ax.set_xlabel('Tid')
    ax.set_ylabel('Amplitude')
    ax.plot(Tid,Sweep)
    format_axes(ax)
    plt.savefig('../Grafik/Sweep.pdf')
    plt.show()

def plotsignal(signal,f_s,unwrap=False,newfig=True):
    freq=np.fft.rfftfreq(len(signal),1./f_s)
    F=np.fft.rfft(signal)
    
    plt.figure()
    gsx,gsy = (3, 1)
    gridspec = plt.GridSpec(gsx, gsy)
    gridspec.update(hspace=0.7)
    
    ax1 = plt.subplot(gridspec[0,:])
#    plt.setp(ax1.get_xticklabels(), visible=False)
    plt.stem(signal)
    ax1.set_xticks(range(len(signal)))
    plt.xlabel("Time (sec)")
    plt.ylabel("x")
        
    ax2 = plt.subplot(gridspec[1,:])
    plt.plot(freq, np.abs(F))
    plt.xlabel("Freq (Hz)")
    plt.ylabel("Amplituden")

    ax3 = plt.subplot(gridspec[2,:],sharex=ax2)
    ax3.set_xticks(freq)
    if unwrap:
        plt.plot(freq, np.unwrap(np.angle(F)))
    else:
        plt.plot(freq, (np.angle(F)))
        plt.ylim((-np.pi,np.pi))
#    plt.setp(ax3.get_xticklabels(), visible=False)
    
    plt.xlabel("Freq (Hz)")
    plt.ylabel("Angle")

#    fig.suptitle('An overall title', size=20)
    plt.setp(ax1, title='Tid')
    plt.setp(ax2, title='Amplituden af Frekvens')
    plt.setp(ax3, title='Phasen af Frekvens')
    


if __name__ == '__main__':
    pass
        
    filenames=get_filenames('./Data/360/')
    print len(filenames),filenames[5]
    print [i for i in filenames if 'Opt' not in i]
#    print [i for i in filenames if 'Opt' in i]
