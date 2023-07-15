'''
This code is used for simulate the thermal photon power spectrum create form each dilution fridge tempmerature stage.
Also, it consider the attenuator attenuate the power from upper stage thermal photon. So it can modified how much 
attenuator place in each stage. Let upper thermal photon power spectrum is lower than 10mK thermal photon power spectrum.

Future function:
    1. add the coaxial cable induce heat, current flow induce heat
    2. simulation the photon population at function noise_photon

ref paper:https://epjquantumtechnology.springeropen.com/articles/10.1140/epjqt/s40507-019-0072-0
Author: Jay Chien
'''

import numpy as np
import matplotlib.pyplot as plt
from tqdm.auto import tqdm
from time import time

def coaxial_heat(r_out, r_dielectirc,r_center,  T1, T2, length):
    #T1 was upper stage, T2 was lower stage
    #r_out was thermal conductivities of outer conductor

    heat_flow = (r_out+r_dielectirc+r_center)*(T1-T2)/length
    return heat_flow

def JN_noise(temp, freq):
    """This term return the Johnson-Nyquist noise 

    Parameters
    ----------
    temp: 
        Temperatire at which enviroment
    freq: 
        noise at which frequency
        
        
    """
    k_b = 1.380649e-23      #Bloteman const
    h = 6.626070153e-34     #Plank const
    R = 50                  #Resistor, in microwave line commonly 50ohm
    const = 4*k_b*temp*R
    a = h*freq/(k_b*temp)
    b = np.exp(h*freq/(k_b*temp))-1 #會有inf的問題, np.isinf判斷
    s = const*a/b
    return np.array(s)
    
def db2power(dB):
    """ Convert dB to how power ratial. Convert function is power unit. Also means how 
    much photon(1/A %) can transmit
     
    Parameters
    ----------
    dB: 
        How much dB of attenuator or other resistor
 
    """

    A=10**(dB/10) 
    return 1/A
def power2dB(power):
    """ Convert power ratial to dB

    Parameters
    ----------
    power:
        how much power ratial 
    
    """
    A = 10*np.log10(power)
    return A

    
def BE_dist(temp, freq):
    """ return the bose einstein distribution with temperature and frequecy
    
    Parameters
    ----------
    temp: 
        Temperatire at which enviroment
    freq: 
        Frequency

    """
    k_b = 1.380649e-23
    h = 6.626070153e-34
    n_BE = 1/(np.exp(h*freq/(k_b*temp))-1)
    return n_BE

def noise_photon(freq, att):
    ''' calculation the noise phton at i stage




    '''
    stage = [300, 50, 4, 0.8, 0.1, 10e-3]

    A = np.zeros(len(att))
    for i in range(len(att)):
        A[i] = 10**(att[i]/10)
    noise = np.zeros(len(stage))
    for i, element in enumerate(stage):
        n_300 = BE_dist(freq, 300)
        noise[0]= n_300
        if i != 0:
            noise[i] = noise[i-1]/A[i] + (A[i]-1)/A[i]*BE_dist(freq, element)
    return noise

def calculate_thermalSA(stage_4K, stage_800mK, stage_100mk, stage_10mK, freq):
    """ Calculation the thermal phton power spectrum
    
    
    """
    stage_check = {"4K":stage_4K, "800mK":stage_800mK, "100mK":stage_100mk, "10mK":stage_10mK}# each stage attenuator
    stage = [300,4,0.8,0.1,10e-3]
    #---------------------------------------------------#
    # The original thermal noise power spectrum
    stage_noise=np.vstack([np.zeros(len(freq))]*len(stage))

    for i, temp in enumerate(stage):
        stage_noise[i,:] = JN_noise(temp, freq)

    fig, ax = plt.subplots(figsize=(8,8))
    color_list = ["orange", "blue", "red", "green","black"]

    for i, temp in enumerate(stage):
        ax.plot(f, stage_noise[i],  color = color_list[i], label=f'{temp}K')

    #---------------------------------------------------#
    # After the attenuator thermal noise power spectrum
    stage.remove(10e-3) #mixing chamber no attenaution
    if stage_check["800mK"] == 0:
        #if 800mK stage no attenuator, ignore 800mK
        del stage_check["800mK"]
        stage.remove(0.8)
        stage_noise_att=np.vstack([np.zeros(len(freq))]*(len(stage)))
    att_list=np.zeros(len(stage_check))

    for i in range(len(stage_check)):
        att_list[i] = np.sum(np.array(list(stage_check.values()))[::-1][:len(stage_check)])

    for i, temp in enumerate(stage):
        if i < len(stage_check):
            stage_noise_att[i,:] = JN_noise(temp, f)*db2power(att_list[i])


    for i, temp in enumerate(stage):
            ax.plot(f, stage_noise_att[i],  color = color_list[i], label=f'{temp}K att', linestyle="dashed")


    plt.axvline(1e10, linestyle=":")
    ax.legend()
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(1e6, 1e13)
    ax.set_ylim(1e-30, 1e-17)
    plt.show()
    pass



if __name__=="__main__":
    f = np.linspace(1e6, 1e13, int(1e5))
    # start = time()
    calculate_thermalSA(20,0,20,20,f)  #4K, 800mK, 100mK, 10mK
    # end = time()
    # print(end-start)



