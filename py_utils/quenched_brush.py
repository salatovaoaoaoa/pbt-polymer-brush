import math
from math import sqrt
import numpy as np
import scipy.optimize
from numpy import pi
import scipy.integrate as integrate
from numpy import trapz
from scipy.optimize import brentq
import pandas as pd


def quecnhed_br(
        N: int = 400,
        S: float = 100,
        alpha: float = 0.5,
        Cs: float = 0.001,
        lb: float = 1,
        a: float = 1,
        
        PK_MINUS: float = 7.0,
        PK_PLUS: float = 7.0,
        f_plus: float = 0.5,
        
        pH_B: float = 5.0,

        file_names = 'annealing_brush_temp.pro',
        way = '/home/tpopova/prj/polymer_brush/SCF_AT/SCF_scripts/output/'
):
    #                                   CONST по pH
    #IEP of protein

    def find_IEP(Pk_, Pk_plus):

        X_iep = (10**(Pk_plus-Pk_))/2 * (2*f_plus-1)/(1-f_plus) \
                + np.sqrt(((10**(Pk_plus-Pk_))**2)/4 * ((2*f_plus-1)/(1-f_plus))**2
                           + (10**(Pk_plus-Pk_) * f_plus)/(1-f_plus))
        
        return math.log10(X_iep) + Pk_

    pH_iep_in_quen = find_IEP(PK_MINUS, PK_PLUS)
    
    #delta for protein
    
    d_pH_b: float = pH_B - pH_iep_in_quen
    
    #разница между pK

    Delta = (PK_PLUS - PK_MINUS)/2

    #                             ПОСТОЯННЫЕ ПАРАМЕТРЫ ЩЕТКИ
    # Обратная длина Дебая
    K: float = sqrt(8 * pi * lb * Cs)

    # Характеристическая длина
    H_0: float = sqrt(8 / (3 * math.pi**2)) * N * sqrt(alpha)/a

    dzeta: float = 2 * sqrt(8 / 3) * alpha ** 1.5 * lb * N ** 2 * a / S

    def dzeta_from(h):
        return h + 1 / 4 * ((np.sqrt((K * H_0) ** 2 + h ** 2)) + h) ** 2 \
               * np.exp(h ** 2) * integrate.quad(lambda t: np.exp(-t ** 2), 0, h)[0]\
               - 1 / 4 * ((np.sqrt((K * H_0) ** 2 + h ** 2)) - h) ** 2 * np.exp(-h ** 2)\
               * integrate.quad(lambda t: np.exp(t ** 2), 0, h)[0]

    h_from_dzeta_m = brentq(lambda h: dzeta_from(h) - dzeta, 0, 10)

    H_q = H_0 * h_from_dzeta_m # type: ignore

    Lambda = H_0 ** 2 / H_q

    z_in_range_q = np.linspace(0, H_q, num=500)
    z_out_range_q = np.linspace(H_q, H_q + 70, num=100)

    # Потенциал

    # 1_out_the_brush (z>H)
    def psi_out(z):
        first = K * Lambda + np.sqrt((K * Lambda) ** 2 + 1) - 1
        second = (K * Lambda - np.sqrt((K * Lambda) ** 2 + 1) + 1) * np.exp(-K * (z - H_q))
        return -2 * np.log((first + second) / (first - second))

    # 2_on_the_boundary (z=H)
    def psi_in_h():
        return 2 * np.log((np.sqrt((K * Lambda) ** 2 + 1) - 1) / (K * Lambda))

    # 3_in_the_brush (0<z<H)
    def psi_in_z(z):
        return (z ** 2 - H_q ** 2) / (H_0 ** 2) + psi_in_h()

########################################################
    y_in = psi_in_z(z_in_range_q)
    y_out = psi_out(z_out_range_q)
########################################################
    #через \delta pH_b

    def alpha_buf_plus_exp(delta_pH_b, pH_iep, Pk_plus):
        return (1 + 10**(delta_pH_b + pH_iep - Pk_plus))**(-1)

    def alpha_buf_minus_exp(delta_pH_b, pH_iep, Pk_):
        return (1 + 10**(Pk_ - delta_pH_b - pH_iep))**(-1)

    def delta_F_ion_exp(psi, alpha_buf_plus_exp, alpha_buf_minus_exp):
        return f_plus * np.log(np.exp(psi)/(alpha_buf_plus_exp + (1 - alpha_buf_plus_exp) * np.exp(psi))) +\
               (1 - f_plus) * np.log(np.exp(-1 * psi)/(alpha_buf_minus_exp + (1 - alpha_buf_minus_exp) * np.exp(-1 * psi)))

    y_exp_in = delta_F_ion_exp(psi_in_z(z_in_range_q),
                               alpha_buf_plus_exp(d_pH_b, pH_iep_in_quen, PK_PLUS),
                               alpha_buf_minus_exp(d_pH_b, pH_iep_in_quen, PK_MINUS))
    y_exp_out = delta_F_ion_exp(psi_out(z_out_range_q),
                                alpha_buf_plus_exp(d_pH_b, pH_iep_in_quen, PK_PLUS),
                                alpha_buf_minus_exp(d_pH_b, pH_iep_in_quen, PK_MINUS))

    def Q_exp(alpha_buf_plus_e, alpha_buf_minus_e, psi):
        return f_plus * (1 + (1 - alpha_buf_plus_e)/alpha_buf_plus_e * np.exp(psi))**(-1)\
               - (1 - f_plus) * (1 + (1 - alpha_buf_minus_e)/alpha_buf_minus_e * np.exp(-1 * psi))**(-1)

    y_exp_q_in = Q_exp(alpha_buf_plus_exp(d_pH_b, pH_iep_in_quen, PK_PLUS),
                               alpha_buf_minus_exp(d_pH_b, pH_iep_in_quen, PK_MINUS),
                       psi_in_z(z_in_range_q))

    y_exp_q_out = Q_exp(alpha_buf_plus_exp(d_pH_b, pH_iep_in_quen, PK_PLUS),
                       alpha_buf_minus_exp(d_pH_b, pH_iep_in_quen, PK_MINUS),
                       psi_out(z_out_range_q))

    #Свободная энергия в максимуме

    def f_ion_max():
        return 0.5 * np.log((10**(-Delta) + 10**(Delta) + 10**(-d_pH_b) + 10**(d_pH_b))/(10**(-Delta) + 10**(Delta) + 2))
    
    F_ion_max = f_ion_max()
    
    #Плотность полимера
    def c_polymer(psi):
        rho = -1 * (2*pi*lb*H_0**2)**(-1)
        c_plus = Cs * np.exp(-psi)
        c_minus = Cs * np.exp(psi)
        
        return (-rho + c_plus - c_minus)/alpha
    
    c_p = c_polymer(psi_in_z(z_in_range_q))
    
    #Свободная энергия через SCF

    #NAMICS
    # file = f'{way}{file_names}'
    file = f'{file_names}'
    parse_SCF_psi_quen = pd.read_csv(file, sep='\t')['sys_noname_psi']
    parse_SCF_phi_quen = pd.read_csv(file, sep='\t')['mol_brush_phi']
    
    f_ion_SCF_quen = delta_F_ion_exp(parse_SCF_psi_quen,
                               alpha_buf_plus_exp(d_pH_b, pH_iep_in_quen, PK_PLUS),
                               alpha_buf_minus_exp(d_pH_b, pH_iep_in_quen, PK_MINUS))
    
    Q_SCF_quen = Q_exp(alpha_buf_plus_exp(d_pH_b, pH_iep_in_quen, PK_PLUS),
                               alpha_buf_minus_exp(d_pH_b, pH_iep_in_quen, PK_MINUS),
                      parse_SCF_psi_quen)

    return H_q, Lambda, \
        d_pH_b, pH_iep_in_quen, \
        z_in_range_q, z_out_range_q, y_in, y_out, y_exp_in, y_exp_out, y_exp_q_in, y_exp_q_out, \
            f_ion_SCF_quen, Q_SCF_quen, parse_SCF_psi_quen, parse_SCF_phi_quen
