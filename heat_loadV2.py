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
from time import time
from scipy.constants import h, Boltzmann
from matplotlib.ticker import MultipleLocator, FormatStrFormatter


def coaxial_heat(r_out, r_dielectirc, r_center,  T1, T2, length):
    # T1 was upper stage, T2 was lower stage
    # r_out was thermal conductivities of outer conductor

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
    # k_b = 1.380649e-23      #Bloteman const
    # h = 6.626070153e-34     #Plank const
    R = 50  # Resistor, in microwave line commonly 50ohm
    const = 4*Boltzmann*temp*R
    a = h*freq/(Boltzmann*temp)
    b = np.exp(h*freq/(Boltzmann*temp))-1  # 會有inf的問題, np.isinf判斷
    s = const*a/b
    return np.array(s)


def inversJN(value, freq):
    return (h*freq)/(Boltzmann*np.log(1+(4*50*h*freq/value)))


def db2power(dB):
    """ Convert dB to how power ratial. Convert function is power unit. Also means how 
    much photon(1/A %) can transmit

    Parameters
    ----------
    dB: 
        How much dB of attenuator or other resistor

    """

    A = 10**(dB/10)
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


def noise_photon(freq, att: list):
    ''' calculation the noise phton at i stage


    '''
    stage = [300, 50, 4, 0.8, 0.1, 10e-3]

    A = np.zeros(len(att))
    for i in range(len(att)):
        A[i] = 10**(att[i]/10)
    noise = np.zeros(len(stage))

    for i, element in enumerate(stage):
        n_300 = BE_dist(freq, 300)
        noise[0] = n_300
        if i != 0:
            noise[i] = noise[i-1]/A[i] + (A[i]-1)/A[i]*BE_dist(freq, element)
    return noise


def calculate_thermalSA(stage_50K, stage_4K, stage_800mK, stage_100mk, stage_20mK, freq, interest_freq):
    """ Calculation the thermal phton power spectrum


    """
    stage_check = {
        "50K": stage_50K,
        "4K": stage_4K,
        "800mK": stage_800mK,
        "100mK": stage_100mk,
        "20mK": stage_20mK}  # each stage attenuator
    stage = [300, 50, 4, 0.8, 0.1, 10e-3]
    # ---------------------------------------------------#
    # The original thermal noise power spectrum
    stage_noise = np.vstack([np.zeros(len(freq))]*len(stage))

    fig, ax = plt.subplots(figsize=(12, 8))
    color_list = ["b", "g", "c", "m", "y", "k"]

    for i, temp in enumerate(stage):
        stage_noise[i, :] = JN_noise(temp, freq)
        ax.plot(f, stage_noise[i],  color=color_list[i], label=f'{temp}K')

    # ---------------------------------------------------#
    # After the attenuator thermal noise power spectrum

    stage_noise_att = np.vstack([np.zeros(len(freq))]*(len(stage)))
    att_list = np.zeros(len(stage_check))

    for i in range(len(stage_check)):
        att_list[i] = sum(list(stage_check.values())[i:])

    stage.remove(stage[-1])
    for i, temp in enumerate(stage):
        if i < len(stage_check):
            stage_noise_att[i, :] = JN_noise(temp, f)*db2power(att_list[i])
        ax.plot(freq, stage_noise_att[i],  color=color_list[i],
                label=f'{temp}K att', linestyle=(0, (3, 5, 1, 5, 1, 5)))
    total_temperature = sum(stage_noise_att)+stage_noise[-1]
    ax.plot(freq, total_temperature, label='effect temperature', lw=3, color='r')

    efft_temp = inversJN(total_temperature[np.where(
        freq == interest_freq)], interest_freq).item()

    ax.plot(freq, JN_noise(efft_temp, freq),
            label=f'JN noise at {efft_temp*1e3:.0f}mK', ls='-.', lw=3)
    # ax.fill_between(freq, JN_noise(efft_temp-0.01, freq),
    #                 JN_noise(efft_temp+0.01, freq), alpha=0.3)
    ax.axvline(interest_freq, linestyle="--", lw=2)
    plt.title(
        f'Effective temperature at {interest_freq/1e9:.2f} GHz = {efft_temp*1e3:.0f}mK')

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylabel('Power')
    ax.set_xlabel('Frequency [Hz]')
    ax.set_xlim(5e7, 1e13)
    ax.set_ylim(1e-30, 1e-17)
    plt.legend(bbox_to_anchor=(1.02, 1), borderaxespad=0)
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    f = np.arange(0, 1e13, 5e7)
    # 50K, 4K, 800mK, 100mK, 20mK
    calculate_thermalSA(11, 11, 11, 11, 0, f, 0.3e9)
    # a = noise_photon(6e9, [0, 11, 11, 11, 11, 21])
    # print(a)
