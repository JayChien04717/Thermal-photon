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
    b = np.exp(h*freq/(k_b*T))-1 #會有inf的問題, np.isinf判斷
    s = const*a/b
    return np.array(s)
    
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

def calculate_thermalSA(stage_4K, stage_800mK, stage_100mk, stage_10mK, freq):
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
        del stage_check["800mK"]
        stage.remove(0.8)
        stage_noise_att=np.vstack([np.zeros(len(freq))]*(len(stage)))
    att_list=np.zeros(len(stage_check))


    for i in range(len(stage_check)):
        att_list[i] = np.sum(np.array(list(stage_check.values()))[::-1][:len(stage_check)])

    for i, temp in enumerate(stage):
        if i < len(stage_check):
            stage_noise_att[i,:] = JN_noise(temp, f)*dB_calculator(att_list[i])


    for i, temp in enumerate(stage):
            ax.plot(f, stage_noise_att[i],  color = color_list[i], label=f'{temp}K', linestyle="dashed")


    plt.axvline(1e10, linestyle=":")
    ax.legend()
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(1e6, 1e13)
    ax.set_ylim(1e-30, 1e-17)
    plt.show()
    pass


if __name__=="__main__":
    att=[0, 0, 20, 20, 0 ,20] # attenautor at 50K, 4K, 800mK, 100mK, 10mK
    f = np.linspace(1e6, 1e13, int(1e5))
    # start = time()
    calculate_thermalSA(20,0,20,20,f)  #4K, 800mK, 100mK, 10mK
    # end = time()
    # print(end-start)

 
