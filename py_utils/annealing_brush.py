from math import sqrt
from math import pi
from numpy import exp
import numpy as np
import math
import scipy.integrate as integrate
from scipy.optimize import brentq
import pandas as pd

def annealing_anion(
        N: int = 300, 
        S: float = 100,
        Cs: float = 0.001,
        lb: float = 1,
        a: float = 1, 

        delta_pK: float = -0.3, #Отступ от pK белка
        pK: float = 5,

        #Параметры белка
        f_plus: float = 0.5, #доля заряженных групп на поверхности
        pK_plus: float = 5,
        pK_minus: float = 5,
        file_name = 'annealing_brush_temp.pro',
        # way = '/home/tpopova/prj/polymer_brush/Free_energy_SCF/SCF_scripts/output/'
        ):
    #                                   CONST по pH
    #IEP of protein

    def find_IEP(Pk_, Pk_plus):

        X_iep = (10**(Pk_plus-Pk_))/2 * (2*f_plus-1)/(1-f_plus) \
                + np.sqrt(((10**(Pk_plus-Pk_))**2)/4 * ((2*f_plus-1)/(1-f_plus))**2
                           + (10**(Pk_plus-Pk_) * f_plus)/(1-f_plus))
        
        return math.log10(X_iep) + Pk_

    pH_iep = find_IEP(pK_minus, pK_plus)

    #pH of solution

    pH_b: float = pK - delta_pK
    
    #delta for protein
    
    delta_pH_b: float = pH_b - pH_iep
    
    #pH_sfbox
    
    pH_sfbox: float = (1 + 10**(pH_b))**(-1)

    #                             ПОСТОЯННЫЕ ПАРАМЕТРЫ ЩЕТКИ

    #степень ионизации на бесконечности

    alpha_b: float = (1 + 10**(delta_pK))**(-1)

    #Обратная длина Дебая

    K: float = sqrt(8 * pi * lb * Cs)

    #Характеристическая толщина щетки

    H_0: float = sqrt(8 / (3 * pi**2)) * N * sqrt(alpha_b)/a

    #dzeta bulk
    dzeta_b: float = (2 * pi * lb * alpha_b * N * H_0)/S

    #Находим alpha_h и толщину щетки h

    def dzeta_func(alpha_h):

        h = ((alpha_h - alpha_b) * K * H_0 * alpha_h) / (2 * alpha_b * np.sqrt((alpha_b - alpha_b * alpha_h) * (alpha_h - alpha_h * alpha_b)))
                                                        
        return -1*(alpha_b**2 * integrate.quad(lambda t: (1 - (1 - alpha_h) * (1 + 2*alpha_b*t**2) * np.exp(alpha_b * (h**2 - t**2)))/(1 \
            - (1 -alpha_h) * np.exp(alpha_b * (h**2 - t**2)))**3, 0, h)[0] + ((K * H_0)**2 * alpha_b)/(4) * integrate.quad(lambda t:\
            (alpha_b/(1 - alpha_b) * ((1 - alpha_h) * np.exp(alpha_b * (h**2 - t**2)))/(1 - (1 - alpha_h)*np.exp(alpha_b*(h**2 - t**2)))**2)\
            - ((1 - alpha_b)/alpha_b * (np.exp(-1 * alpha_b*(h**2 - t**2)))/(1 - alpha_h)), 0, h)[0])
    
    start = 0.0001

    while True:
        try:
            result = brentq(lambda alpha_h: dzeta_func(alpha_h) - dzeta_b, alpha_b - start, alpha_b)
            break
        except:
            start += 0.001

    alpha_H: float = result # type: ignore

    h_find = -1 * (((alpha_H - alpha_b) * K * H_0 * alpha_H) / (2 * alpha_b * np.sqrt((alpha_b - alpha_b * alpha_H) * (alpha_H - alpha_H * alpha_b))))

    H = h_find * H_0

    tlambda = (alpha_H * H_0**2)/(alpha_b * H)

    #Координата

    z_in_range = np.linspace(0, H, num=500)
    z_out_range = np.linspace(H, H + 30, num=200)

    def c_p(z):
        const1 = alpha_b/(2 * np.pi * lb  * H_0**2)
        const2 = (K**2 * H_0**2)/(4 * alpha_b)
        first = (1 - (1 - alpha_H) * (1 + 2*alpha_b*z**2/H_0**2) * np.exp(alpha_b * (H**2 - z**2)/H_0**2))/(1 - \
                (1 - alpha_H) * np.exp(alpha_b * (H**2 - z**2)/H_0**2))**3
        second = alpha_b/(1 - alpha_b) * ((1 - alpha_H) * np.exp(alpha_b * (H**2 - z**2)/H_0**2))/(1 -\
                (1 - alpha_H) * np.exp(alpha_b * (H**2 - z**2)/H_0**2))**2
        third = (1 - alpha_b)/alpha_b * (np.exp(-1 * alpha_b * (H**2 - z**2)/H_0**2))/(1 - alpha_H)
        return const1 * (first + const2 * (second - third))
    
    c_polymer = c_p(z_in_range)

    #alpha z

    def alpha_from_z(z):
        return 1 - (1 - alpha_H) * np.exp((alpha_b * (H**2 - z**2)) / (H_0**2))
    
    alpha_z = alpha_from_z(z_in_range)

    alpha_z_mean = (integrate.quad(lambda z: alpha_from_z(z) * c_p(z), 0, H)[0])/(N/S)
    
    #                              ЭЛЕКТРОСТАТИЧЕСКИЙ ПОТЕНЦИАЛ

    #на границе щетки

    psi_in_H = np.log((alpha_H * (1 - alpha_b))/(alpha_b * (1 - alpha_H)))

    # внутри щетки

    def psi_in_z(z):
        return -1 * alpha_b * (H**2 - z**2)/(H_0**2) \
            + np.log((1 - alpha_b)/(alpha_b * (1 - alpha_H)) * (1 - (1 - alpha_H) * np.exp(alpha_b * (H**2 - z**2)/(H_0**2))))
    
    psi_in = psi_in_z(z_in_range)

    #снаружи щетки

    def psi_out_z(z):

        minus_one = K * tlambda + np.sqrt((K * tlambda)**2 + 1) - 1
        plus_one = K * tlambda - np.sqrt((K * tlambda)**2 + 1) + 1
        exp = np.exp(-1 * K * (z - H))

        return -2 * np.log((minus_one + plus_one * exp)/(minus_one - plus_one * exp))
    
    psi_out = psi_out_z(z_out_range)

    # на поверхности прививки

    def psi_in_0():
        return -1 * alpha_b * (H**2)/(H_0**2) \
            + np.log((1 - alpha_b)/(alpha_b * (1 - alpha_H)) * (1 - (1 - alpha_H) * np.exp(alpha_b * (H**2)/(H_0**2))))
    
    psi_in_zero = psi_in_0()

    #               СТЕПЕНЬ ИОНИЗАЦИИ В БУФЕРЕ

    def alpha_buf_plus(pH_b, pK_plus_protein):
        return (1 + 10**(pH_b - pK_plus_protein))**(-1)

    def alpha_buf_minus(pK_minus_protein, pH_b):
        return (1 + 10**(pK_minus_protein - pH_b))**(-1)
    
    #             СТЕПЕНЬ ИОНИЗАЦИИ В ЩЕТКЕ

    def alpha_plus_in():
        return (1 + (1 - alpha_buf_plus(pH_b, pK_plus))/(alpha_buf_plus(pH_b, pK_plus))
                * np.exp(psi_in_z(z_in_range)))**(-1)

    def alpha_plus_out():
        return (1 + (1 - alpha_buf_plus(pH_b, pK_plus))/(alpha_buf_plus(pH_b, pK_plus))
                * np.exp(psi_out_z(z_out_range)))**(-1)

    def alpha_minus_in():
        return (1 + (1 - alpha_buf_minus(pK_minus, pH_b))/(alpha_buf_minus(pK_minus, pH_b))
                * np.exp(-1 * psi_in_z(z_in_range)))**(-1)

    def alpha_minus_out():
        return (1 + (1 - alpha_buf_minus(pK_minus, pH_b))/(alpha_buf_minus(pK_minus, pH_b))
                * np.exp(-1 * psi_out_z(z_out_range)))**(-1)
    
    alpha_pl_in = alpha_plus_in()
    alpha_pl_out = alpha_plus_out()
    alpha_min_in = alpha_minus_in()
    alpha_min_out = alpha_minus_out()

    # Свободная энергия

    def alpha_buf_plus_fast():
        return (1 + 10**(delta_pH_b + pH_iep - pK_plus))**(-1)

    def alpha_buf_minus_fast():
        return (1 + 10**(pK_minus - delta_pH_b - pH_iep))**(-1)

    def delta_F_ion(psi, alpha_buf_plus_fast, alpha_buf_minus_fast):
        return f_plus * np.log(np.exp(psi)/(alpha_buf_plus_fast + (1 - alpha_buf_plus_fast) * np.exp(psi))) +\
               (1 - f_plus) * np.log(np.exp(-1 * psi)/(alpha_buf_minus_fast + (1 - alpha_buf_minus_fast) * np.exp(-1 * psi)))

    f_ion_in = delta_F_ion(psi_in_z(z_in_range),
                              alpha_buf_plus_fast(),
                               alpha_buf_minus_fast())
    f_ion_out = delta_F_ion(psi_out_z(z_out_range),
                                alpha_buf_plus_fast(),
                               alpha_buf_minus_fast())
    
    #Заряд

    def Q_exp(alpha_buf_plus_fast, alpha_buf_minus_fast, psi):
        return f_plus * (1 + (1 - alpha_buf_plus_fast)/alpha_buf_plus_fast * np.exp(psi))**(-1)\
               - (1 - f_plus) * (1 + (1 - alpha_buf_minus_fast)/alpha_buf_minus_fast * np.exp(-1 * psi))**(-1)

    charge_in = Q_exp(alpha_buf_plus_fast(),
                        alpha_buf_minus_fast(),
                       psi_in_z(z_in_range))

    charge_out = Q_exp(alpha_buf_plus_fast(),
                       alpha_buf_minus_fast(),
                       psi_out_z(z_out_range))
    
    # #Плотность заряда rho (z)
    
    def rho(z):
        memb = -1 * alpha_b/(2*pi*lb*(H_0)**2)
        numerator = 1 - (1-alpha_H) * (1 + 2*alpha_b*z**2/H_0**2)*exp(alpha_b*(H**2-z**2)/H_0**2)
        denominator = (1 - (1-alpha_H) * exp(alpha_b*(H**2-z**2)/H_0**2))**2
        return memb * (numerator)/(denominator)
    
    rho_range = [rho(z_in_range[i]) for i in range(len(z_in_range))]

    polymer_dens_anneal = [(-1 * rho_range[i] - Cs * np.exp(psi_in[i]) + Cs * np.exp(-psi_in[i]))/(alpha_z[i]) for i in range(len(z_in_range))]
    
    #Свободная энергия через SCF

    #NAMICS
    # file = f'{way}{file_name}'
    file = f'{file_name}'
    parse_SCF_psi = pd.read_csv(file, sep='\t')['sys_noname_psi']
    parse_SCF_phi = pd.read_csv(file, sep='\t')['mol_brush_phi']
    
    f_ion_SCF = delta_F_ion(parse_SCF_psi,
                              alpha_buf_plus_fast(),
                               alpha_buf_minus_fast())
    charge_SCF = Q_exp(alpha_buf_plus_fast(),
                        alpha_buf_minus_fast(),
                       parse_SCF_psi)
    
    #Средний заряд белка
    
    Q_mean = charge_out[-1]
    
    return H, alpha_H, tlambda, K, alpha_z_mean,alpha_z,alpha_b, \
            delta_pK, pH_b, pH_iep,pH_sfbox, delta_pH_b, \
            z_in_range, z_out_range,psi_in, psi_out, f_ion_in, f_ion_out, charge_in, charge_out,\
            polymer_dens_anneal, \
            f_ion_SCF, charge_SCF, parse_SCF_psi, parse_SCF_phi, Q_mean