"""
File con simulazioni Monte Carlo  
_________________________________________
"""
import numpy as np
import math as mt
import matplotlib.pyplot  as plt
from scipy import integrate
from scipy import optimize
import matplotlib.colors
import sys
import argparse

#Importo modulo con funzioni
sys.path.append(" ")
import func as f
import starclass

"""
Gestione delle azioni con argparse
__________________________________
"""
def parse_arguments():
    
    parser = argparse.ArgumentParser(description='Gestione delle simulazioni',
                                     usage      ='--opzione')
    parser.add_argument('-sole', '--opzione1', action='store_true',  help='Simulazione con Sole')
    parser.add_argument('-betelgeuse', '--opzione2', action='store_true',  help='Simulazione con Betelgeuse')
    parser.add_argument('-bellatrix', '--opzione3', action='store_true',  help='Simulazione con Bellatrix')
    parser.add_argument('-alfacrucis', '--opzione4', action='store_true',  help='Simulazione con Alfa Crucis')
    parser.add_argument('-tutte', '--opzione5', action='store_true',  help='Simulazione con tutte le Stelle')
    return  parser.parse_args()


"""
Simulazione eventi per le Stelle in analisi
____________________________________________
"""
def main():

    #Creo le stelle
    sole=starclass.star("Sole",5.75*10**3)
    betelgeuse=starclass.star("Betelgeuse",3*10**3)
    bellatrix=starclass.star("Bellatrix",22*10**3)
    alfa_crucis=starclass.star("Alfa Crucis",28*10**3)
    stelle=[sole,betelgeuse,bellatrix,alfa_crucis]

    #Indago solo nel visibile
    estremo_sx=380 #nm
    estremo_dx=790 #nm
    
    N=50000 

    #Angoli con cui studiare andamento numero fotoni totali
    angles=np.linspace(0,90,100)
    N1=10000

    args = parse_arguments()

    if args.opzione1 == True:

            #Emissione fotoni 
            sim1=sole.no_absorption(N,estremo_sx,estremo_dx)
            photons=sim1[0]
            k_norm=sim1[1]
            sole.no_abs_graph(photons,k_norm)
            
            #Simulazione scattering Rayleigh allo Zenit
            sim2=sole.R_scattering(photons,0)
            sole.R_scattering_graph(photons,sim2,0,k_norm)

            #Scattering scattering Rayleigh all'orizzonte
            sim3=sole.R_scattering(photons,90)
            sole.R_scattering_graph(photons,sim3,90,k_norm)

            #Confronto distribuzioni nel visibile
            sole.compare_distribution(N)

            #Numero di fotoni al variare angolo
            sole.flusso_integrato(N1,angles)
    
    if args.opzione2 == True:

            #Emissione fotoni 
            sim1=betelgeuse.no_absorption(N,estremo_sx,estremo_dx)
            photons=sim1[0]
            k_norm=sim1[1]
            betelgeuse.no_abs_graph(photons,k_norm)
            
            #Simulazione scattering Rayleigh allo Zenit
            sim2=betelgeuse.R_scattering(photons,0)
            betelgeuse.R_scattering_graph(photons,sim2,0,k_norm)

            #Scattering scattering Rayleigh all'orizzonte
            sim3=betelgeuse.R_scattering(photons,90)
            betelgeuse.R_scattering_graph(photons,sim3,90,k_norm)

            #Confronto distribuzioni nel visibile
            betelgeuse.compare_distribution(N)

            #Numero di fotoni al variare angolo
            betelgeuse.flusso_integrato(N1,angles)
    
    if args.opzione3 == True:
            #Emissione fotoni 
            sim1=bellatrix.no_absorption(N,estremo_sx,estremo_dx)
            photons=sim1[0]
            k_norm=sim1[1]
            bellatrix.no_abs_graph(photons,k_norm)
            
            #Simulazione scattering Rayleigh allo Zenit
            sim2=bellatrix.R_scattering(photons,0)
            bellatrix.R_scattering_graph(photons,sim2,0,k_norm)

            #Scattering scattering Rayleigh all'orizzonte
            sim3=bellatrix.R_scattering(photons,90)
            bellatrix.R_scattering_graph(photons,sim3,90,k_norm)

            #Confronto distribuzioni nel visibile
            bellatrix.compare_distribution(N)

            #Numero di fotoni al variare angolo
            bellatrix.flusso_integrato(N1,angles)

    if args.opzione4 == True:
            #Emissione fotoni 
            sim1=alfa_crucis.no_absorption(N,estremo_sx,estremo_dx)
            photons=sim1[0]
            k_norm=sim1[1]
            alfa_crucis.no_abs_graph(photons,k_norm)
            
            #Simulazione scattering Rayleigh allo Zenit
            sim2=alfa_crucis.R_scattering(photons,0)
            alfa_crucis.R_scattering_graph(photons,sim2,0,k_norm)

            #Scattering scattering Rayleigh all'orizzonte
            sim3=alfa_crucis.R_scattering(photons,90)
            alfa_crucis.R_scattering_graph(photons,sim3,90,k_norm)

            #Confronto distribuzioni nel visibile
            alfa_crucis.compare_distribution(N)

            #Numero di fotoni al variare angolo
            alfa_crucis.flusso_integrato(N1,angles)

    if args.opzione5 ==True:

        for item in stelle:
            print("--------------------------------------------------------------------------------------------")
            print("Simulazione distribuzione senza assorbimento -  {:} ".format(item.name))
            sim1=item.no_absorption(N,estremo_sx,estremo_dx)
            photons=sim1[0]
            k_norm=sim1[1]
            item.no_abs_graph(photons,k_norm)
            print("--------------------------------------------------------------------------------------------")
            print("Simulazione distribuzione con scattering Rayleigh allo zenit -  {:} ".format(item.name))
            sim2=item.R_scattering(photons,0)
            item.R_scattering_graph(photons,sim2,0,k_norm)
            print("--------------------------------------------------------------------------------------------")
            print("Simulazione distribuzione con scattering Rayleigh all'orizzonte -  {:} ".format(item.name))
            sim3=item.R_scattering(photons,90)
            item.R_scattering_graph(photons,sim3,90,k_norm)
            print("--------------------------------------------------------------------------------------------")
            print("Confronto simulazione diverse distribuzioni -  {:} ".format(item.name))
            item.compare_distribution(N)
            print("--------------------------------------------------------------------------------------------")
            print("Studio andamento numero fotoni totali osservati al variare angolo  -  {:} ".format(item.name))
            item.flusso_integrato(N1,angles)
        

if __name__ == "__main__":

    main()