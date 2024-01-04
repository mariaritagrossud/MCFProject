"""
File con le funzioni usate nel progetto 
_________________________________________
"""

import numpy as np
import math as mt

#Parametri uguali per tutte le simulazioni
c = 3*10**8 #m/s
h = 6.63*10**(-34) #J*s
k_b = 1.38*10**(-23) #J/K
n_refr = 1.00029
density = 2.504*(10**(25)) #molecole/m**3
R_T = 6378*10**3 #m
S_z = 8000 #m


def D(lam,T):

    """
    Funzione densità di fotoni per lunghezza d'onda e unità di superficie e tempo
    
    D(lamda,T) =2*c/((lamda^4)*(e^(alfa/lambda)-1))

    con costanti: 
    c = velocità della luce nel vuoto
    alfa = h*c/k_B*T con
           h = costante di Planck
           k_B = costante di Boltzmann 
           T = temperatura della Stella
    
    input:
    lam = array di lunghezze d'onda espresse in nanometri
    T = temperatura della Stella di cui si sta studiando lo spettro

    output:
    Array dei valori della funzione nelle lunghezze d'onda in ingresso

    """

    #Lunghezza d'onda in entrata in nanometri, converto in metri
    lam_fromnano = lam*10**(-9)
    alfa = h*c/(k_b*T)
    expon = alfa/lam_fromnano
    den=np.power(lam_fromnano, 4)*(np.exp(expon)-1)
    return 2*c/den

def D_norm(lam,T,k_norm):

    """
    Funzione densità di fotoni per lunghezza d'onda e unità di superficie e tempo normalizzata per avere
    significato di probabilità

    D_norm(lambda,T,k_norm) = D(lambda,T)*k_norm

    con k_norm la costante di normalizzazione, ossia il reciproco dell'integrale della distribuzione 
    nell'interavallo di lunghezze d'onda considerate

    """
    return D(lam,T)*k_norm

def D_norm_opposite(lam,T,k_norm):
    """
    Funzione opposta della densità di fotoni per lunghezza d'onda e unità di superficie e tempo normalizzata
    usata per trovare il massimo di D_norm

    """
    return -D_norm(lam,T,k_norm)

def n_obs(lam,T,angle):

    """
    Funzione densità di fotoni per lunghezza d'onda e unità di superficie e tempo considerando scattering di 
    Rayleigh
    n_obs(lam,T,angle) = D(lam,T)*e^(beta(lam)*S(angle))

    con:
    beta(lam) = (8*π^3*(n^2-1)^2)/(3*N*lam^4)
        n = indice di rifrazione dell'aria
        N = densità di molecole nell'aria
    
    S(angle) = ((R_T*cos(angle))**2+2*R_T*S_z+S_z**2)^(0.5)-R_T*.cos(angle)
        R_T = raggio terrestre
        S_z = spessore della massa d'aria attraversata dai fotoni solari quando il Sole è allo Zenit
    
    input:
    lam = array di lunghezze d'onda espresse in nanometri
    T = temperatura della Stella di cui si sta studiando lo spettro
    angle = angolo del Sole rispetto allo Zenit
    """
    #beta
    lam_fromnano=lam*10**(-9)
    num=8*mt.pi**3*(n_refr**2-1)**2
    den=3*density*(lam_fromnano**4)
    beta= num/den

    #S
    angle_rad=(angle*mt.pi)/180
    S= mt.sqrt((R_T*mt.cos(angle_rad))**2+2*R_T*S_z+S_z**2)-R_T*mt.cos(angle_rad)

    return D(lam, T)*np.exp(-beta*S)

def n_obs_norm(lam,T,theta,k_norm):
    """
    Funzione n_obs normalizzata 
    
    """
    return n_obs(lam,T,theta)*k_norm

def n_obs_norm_opposite(lam,T,theta,k_norm):

    return -n_obs_norm(lam,T,theta,k_norm)

#definisco D e n_obs con entrata lunghezze in micrometri poiché nel fit dà problemi

def n_photons_fit(lam,k,T):
    """
    Definisco la funzione analoga a n_obs per il fit per cui si ha bisogno anche di 
    un coefficiente motliplicativo che tenga in conto che la distribuzione non è normalizzata
    """
    lam_fromnano=lam*10**(-9)
    alfa=h*c/(k_b*T)
    den=np.power(lam_fromnano, 4)*(np.exp(alfa/lam_fromnano)-1)
    D_value=2*c/den
    
    theta_xrad=45*mt.pi/180
    S_theta=mt.sqrt((R_T*mt.cos(theta_xrad))**2+2*R_T*S_z+S_z**2)-R_T*mt.cos(theta_xrad)
    beta=8*np.pi**3*(n_refr**2-1)**2/(3*density*lam_fromnano**4)

    return k*D_value*np.exp(-beta*S_theta)