import sys
import os
sys.path.append('/home/tpopova/prj/polymer_brush')

from math import sqrt
from math import exp
from math import pi

import numpy as np 

import scipy.integrate as integrate 
from scipy.optimize import root

from scipy.special import i0 
from scipy.special import i1 

from scipy.special import k0 
from scipy.special import k1 

def pore_final(S : float  = 100, #площадь поверхности поры на цепь
              alpha : float = 0.5,
              Cs : float = 0.001,
              lb : float = 1,
              D : float = 300,
              N : int = 300):
    
    # inverse Debye length
    K: float = sqrt(8 * pi * lb * Cs)
    
    # Characteristic length
    H_0: float = sqrt(8 / (3 * pi**2)) * N * sqrt(alpha) 

    #dzeta b
    zeta_b: float = (2 * pi * D * lb * alpha * N)/S
    
    #Theta
    #The relation of the lateral surface area of the cylinder to the unit length l
    l_t = S/(2*pi*D)

    #theta
    theta = N/l_t

    #We calculate the unlimited thickness of the brush:

    def find_h(h):

        t_lambda = H_0/h

        rho = h * (D/H_0 - h)

        c_plus = (K*H_0/2)**2 * np.exp(h**2 + 2/(K*t_lambda) * i0((D-h*H_0)*K)/i1((D-h*H_0)*K))\
            * (D/H_0 * (integrate.quad(lambda t: np.exp(-t ** 2), 0, h)[0]) - (1 - np.exp(-1 * h**2))/2)

        c_minus = (K*H_0/2)**2 * np.exp(-1 * h**2 - 2/(K*t_lambda) * i0((D-h*H_0)*K)/i1((D-h*H_0)*K))\
            * (D/H_0 * (integrate.quad(lambda t: np.exp(t ** 2), 0, h)[0]) - (np.exp(h**2) - 1)/2)
        
        return rho + c_plus - c_minus - zeta_b
   
    h_solution = root(lambda h: find_h(h), 0.01, method='lm')

    H = h_solution.x * H_0

    t_lambda_answ = H_0/h_solution.x 
    
    #Electrostatic potential

    def psi_out(r, tLambda, H_):
        return -2/(K*tLambda) * (i0(r*K))/(i1((D-H_)*K))
    
    def psi_in(r, tLambda, H_):
        first = ((D-r)**2 - H_**2)/H_0**2
        psi_H = -2/(K*tLambda) * (i0((D-H_)*K))/(i1((D-H_)*K))
        return first + psi_H
    
    #Charge density

    def rho_in(r):
        return -1/(2*pi*lb*H_0**2) * (2*r - D)/r
    
    #Polymer density

    def c_p(r, psi):
        
        rho = rho_in(r)

        c_plus = Cs * np.exp(-1 * psi)

        c_minus = Cs * np.exp(psi)

        return (-1 * rho + c_plus - c_minus)/alpha
    
    if H < D:
    
        #We calculate the ranking of the coordinate according to the unlimited H:

        r_in_range = np.linspace(D, D-H[0], num = 500)
        r_out_range = np.linspace(D-H[0], 0, num = 500)

        #ELECTROSTATIC POTENTIAL THROUGH UNLIMITED THICKNESS

        psi_in_range = psi_in(r_in_range, t_lambda_answ, H)
        psi_out_range = psi_out(r_out_range, t_lambda_answ, H)
        
        c_pol = c_p(r_in_range, psi_in(r_in_range, t_lambda_answ, H))
        
        rho = rho_in(r_in_range)
        
        return H, r_in_range, r_out_range, psi_in_range, psi_out_range, c_pol, rho, theta
    
    else:
        
        def find_h_cut(h_cut):

            t_lambda_cut = H_0/h_cut

            rho_cut = h_cut * (D/H_0 - h_cut)

            c_plus_cut = (K*H_0/2)**2 * np.exp(h_cut**2 + 2/(K*t_lambda_cut) * k0((D-h_cut*H_0)*K)/k1((D-h_cut*H_0)*K))\
                * (D/H_0 * (integrate.quad(lambda t: np.exp(-t ** 2), 0, h_cut)[0]) - (1 - np.exp(-1 * h_cut**2))/2)

            c_minus_cut = (K*H_0/2)**2 * np.exp(-1 * h_cut**2 - 2/(K*t_lambda_cut) * k0((D-h_cut*H_0)*K)/k1((D-h_cut*H_0)*K))\
                * (D/H_0 * (integrate.quad(lambda t: np.exp(t ** 2), 0, h_cut)[0]) - (np.exp(h_cut**2) - 1)/2)
        
            return rho_cut + c_plus_cut - c_minus_cut - zeta_b
   
    h_solution_cut = root(lambda h_cut: find_h_cut(h_cut), 0.01, method='lm').x

    H_cut = h_solution_cut * H_0

    t_lambda_answ_cut = H_0/h_solution_cut
    
    #We calculate the ranking of the coordinate according to the unlimited H:

    r_in_range_cut = np.linspace(D, D-H_cut[0], num = 600)
    r_out_range_cut = np.linspace(D-H_cut[0], 0, num = 300)

    #ELECTROSTATIC POTENTIAL THROUGH UNLIMITED THICKNESS

    psi_in_range_cut = psi_in(r_in_range_cut, t_lambda_answ_cut, H_cut)
    psi_out_range_cut = psi_out(r_out_range_cut, t_lambda_answ_cut, H_cut)
    
    c_pol_cut = c_p(r_in_range_cut, psi_in(r_in_range_cut, t_lambda_answ_cut, H_cut))
    
    rho_cut = rho_in(r_in_range_cut)
        
    return H_cut, r_in_range_cut, r_out_range_cut, psi_in_range_cut, psi_out_range_cut, c_pol_cut, rho_cut, theta
    