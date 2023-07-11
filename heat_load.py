import numpy as np
import matplotlib.pyplot as plt
from tqdm.auto import tqdm
from time import time

def coaxial_heat(r_out, r_dielectirc,r_center,  T1, T2, length):
    #T1 was upper stage, T2 was lower stage
    #r_out was thermal conductivities of outer conductor

    heat_flow = (r_out+r_dielectirc+r_center)*(T1-T2)/length
    return heat_flow

def JN_noise(T, freq):
    k_b = 1.380649e-23
    h = 6.626070153e-34
    R = 50
    const = 4*k_b*T*R
    a = h*freq/(k_b*T)
    b = np.exp(h*freq/(k_b*T))-1
    s = const*a/b
    return s
    
def dB_calculator(dB):
    # in power unit
    A=10**(dB/10)
    return 1/A

    
def BE_dist(freq, T):
    k_b = 1.380649e-23
    h = 6.626070153e-34
    n_BE = np.exp(h*freq/(k_b*T))-1
    return 1/n_BE

def noise_photon(freq, att):
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

def caculate_thermalSA(stage_4K, stage_800mK, stage_100mk, stage_10mK, f):
    stage_check = {"4K":stage_4K, "800mK":stage_800mK, "100mK":stage_100mk, "10mK":stage_10mK}# each stage attenuator
    
    #---------------------------------------------------#
    # The original thermal noise power spectrum
    kelvin_300k = JN_noise(300, f)
    kelvin_4k = JN_noise(4, f)
    kelvin_800mk = JN_noise(0.8, f)
    kelvin_100mk = JN_noise(0.1, f)
    kelvin_10mk = JN_noise(10e-3, f)

    fig, ax = plt.subplots(figsize=(8,8))
    ax.plot(f, kelvin_300k, "orange", label='300K')
    ax.plot(f, kelvin_4k, "blue", label='4K')
    if stage_check["800mK"] != 0:
        ax.plot(f, kelvin_800mk, "r", label='800mk')
    ax.plot(f, kelvin_100mk, "g", label='100mk')
    ax.plot(f, kelvin_10mk, "black", label='10mk')

    #---------------------------------------------------#
    # After the attenuator thermal noise power spectrum
    kelvin_300k_att = JN_noise(300, f)*dB_calculator(stage_4K+stage_800mK+stage_100mk+stage_10mK)
    kelvin_4k_att = JN_noise(4, f)*dB_calculator(stage_800mK+stage_100mk+stage_10mK)
    kelvin_800mk_att = JN_noise(0.8, f)*dB_calculator(stage_100mk+stage_10mK)
    kelvin_100mk_att = JN_noise(0.1, f)*dB_calculator(stage_10mK)
 


    ax.plot(f, kelvin_300k_att, "orange", label='300K att',linestyle="dashed")
    ax.plot(f, kelvin_4k_att, "blue", label='4K att',linestyle="dashed")
    if stage_check["800mK"] != 0:
        ax.plot(f, kelvin_800mk_att, "r", label='800mk att',linestyle="dashed")
    ax.plot(f, kelvin_100mk_att, "g", label='100mk att',linestyle="dashed")

    plt.axvline(1e10, linestyle=":")
    ax.legend()
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(1e6, 1e13)
    ax.set_ylim(1e-30, 1e-17)
    plt.show()
    return


if __name__=="__main__":
    att=[0, 0, 20, 20, 0 ,20] # attenautor at 50K, 4K, 800mK, 100mK, 10mK
    # f = np.linspace(1e6, 1e13, 100001)
    # start = time()
    # caculate_thermalSA(20,0,20,20,f)
    # end = time()
    # print(end-start)

