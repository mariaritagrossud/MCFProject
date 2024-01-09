"""
File con le funzioni e la classe usate nel progetto 
___________________________________________________
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
    expon = (alfa/lam_fromnano)
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

    Funzione opposta della densità di fotoni per lunghezza d'onda e unità di superficie e tempo normalizzata.
    Usata per trovare il massimo di D_norm

    """
    return -D_norm(lam,T,k_norm)

def prob_obs(lam,angle):

    """
    Funzione che restituisce la probabilità di osservazione di un fotone ad una data lunghezza d'onda
    quando la stella è ad una certa posizione rispetto allo Zenit

    Input:
    lam: array di lunghezze d'onda
    angle: angolo della stella rispetto allo Zenit

    p_obs(lam,angle) = e^(-beta(lam)*S(angle))

    con:
    - beta(lam) = 8π^3*(n^2-1)^2/(3*N*lam^4)
        n: indice di rifrazione dell'aria
        N: densità di molecole nell'aria

    - S(angle) = ((R_T*cos(angle))^2+2*R_T*S_z+S_z^2)^0.5)-R_T*cos(angle)
        R_T: raggio della Terra
        S_z: spessore massa d'aria allo Zenit
    
    Output:
    Array con le probabilità
    """
    
    lam_fromnano=np.array(lam)*10**(-9)
    num=-8*mt.pi**3*(n_refr**2-1)**2
    den=3*density*(lam_fromnano**4)
    beta= num/den
    angle_rad=angle*mt.pi/180
    S= mt.sqrt((R_T*mt.cos(angle_rad))**2+2*R_T*S_z+S_z**2)-R_T*mt.cos(angle_rad)
    return np.exp(beta*S)


def n_photons_fit(lam,angle,k,T):

    """
    Funzione uguale a N_obs a meno di un coefficiente moltiplicativo che tenga in conto che la distribuzione non è 
    normalizzata
    N_obs(lam,angle,T) = k * D(lam,T) * p_obs(lam,angle)

    """
    lam_fromnano=lam*10**(-9)
    alfa=h*c/(k_b*T)
    den=np.power(lam_fromnano, 4)*(np.exp(alfa/lam_fromnano)-1)
    D_value=2*c/den
    
    theta_xrad=angle*mt.pi/180
    S_theta=mt.sqrt((R_T*mt.cos(theta_xrad))**2+2*R_T*S_z+S_z**2)-R_T*mt.cos(theta_xrad)
    beta=8*np.pi**3*(n_refr**2-1)**2/(3*density*lam_fromnano**4)

    return k*D_value*np.exp(-beta*S_theta)