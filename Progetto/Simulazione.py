"""
File con simulazioni Monte Carlo  
_________________________________________
"""
import numpy as np
import math as mt
import matplotlib.pyplot  as plt
from scipy import integrate
from scipy import optimize
from scipy import stats
import matplotlib.colors
import sys
import argparse

#Importo modulo con funzioni
sys.path.append(" ")
import func as f


"""
Creo classe con la quale simulare e analizzare la distribuzione radiazione fotoni di varie stelle
__________________________________________________________________________________________________
"""


class star:
    def __init__(self, name, T):
        self.name=name
        self.T=T
    

    def no_absorption(self, N,estremo_sx, estremo_dx):
        """
        Metodo che simula la distribuzione dei fotoni emessi da un corpo nero e che noi vedremmo 
        se non ci fosse lo scattering di rayleigh.
        
        Input:
        N = numero di fotoni da usare per generare la distribuzione
        estremo_sx = più piccola lunghezza d'onda che si vuole generare
        estremo_dx = più grande lunghezza d'onda che si vuole generare

        Output:
        - Array con i fotoni generati dalla distribuzione
        - costate di normalizzazione 
        """
        num=int((estremo_dx-estremo_sx)/0.1) #faccio in modo che i punti siano equispaziati di 0.1 nm
        x=np.linspace(estremo_sx,estremo_dx,num)
        k_norm=1/integrate.simpson(f.D(x,self.T),x) #costante con cui normalizzare distribuzione
        x_max=optimize.minimize(f.D_norm_opposite,x0=(estremo_sx+estremo_dx)/2,args=(self.T,k_norm), tol=1e-9)
        max_value_1=f.D_norm(x_max.x[0],self.T,k_norm)    
        
        #Simulazione con il metodo dell' Hit or Miss
        fotoni=[]
        lam=np.random.uniform(low=estremo_sx,high=estremo_dx,size=N)
        y_value=np.random.uniform(low=0, high=max_value_1, size=N)
        for i in range(N):

             if (y_value[i] <= f.D_norm(lam[i],self.T,k_norm)):
                fotoni.append(lam[i])
        
        return fotoni, k_norm


    def no_abs_graph(self,fotoni,k_norm):
        """
        Metodo con cui si rappresenta graficamente la distribuzione del corpo nero e la si confronta 
        con valori attesi
        """
        #Plot 
        plt.figure(figsize=(10,8))

        n_bins=int(mt.sqrt(len(fotoni)))
        n,bins,p=plt.hist(fotoni,bins=n_bins, color="limegreen", label="Valori simulati")
        
        bincenters=(bins[:-1]+bins[1:])/2
        binwidth=(bins[-1:]-bins[0])/n_bins
        height=f.D_norm(bincenters,self.T,k_norm)
        prob=height*binwidth
        n_exp=len(fotoni)*prob
       
        mask=np.nonzero(n)
        chi2 =  np.sum( (n_exp[mask] - n[mask])**2 /n[mask] ) 
        ndof = len(n[mask])-1
        chi2_red=chi2/ndof

        plt.text(bincenters[1],max(n)*1/4, r'$\chi^2$  : {:.2f}'.format(chi2_red), fontsize=14, color='blue')
        plt.plot(bincenters,n_exp, color="blue", label="Valori attesi", marker=".", linewidth=1)
        plt.title("Stella: {:} \n Distribuzione simulata e attesa dei fotoni emessi dal corpo nero. ".format(self.name), 
                   color="mediumblue",fontsize=13)
        plt.ylabel(r'$ \text{Numero di fotoni osservati per unità  di superficie e tempo} \  \left( \frac{1}{m^2 \cdot t }\right)  $'
                    ,fontsize=9.5)
        plt.xlabel(r'$ \text{Lunghezza d\'onda (nm)} $')
        plt.legend(frameon=False)
        plt.show()

    
    def R_scattering(self, fotoni, angle) :

        """
        Metodo che simula lo scattering di Ryleigh dei fotoni emessi dal corpo nero una volta
        arrivati in atmosfera
        
        Input:
        fotoni = fotoni che sono emessi dal corpo nero e che arrivano nell'atmosfera
        angle = angolo della stella rispetto allo Zenit

        output:
        Array con i fotoni osservati, non diffusi

        """

        prob=f.prob_obs(fotoni,angle)
        n_obs=[]
        for i in range(len(fotoni)):
            y_value=np.random.uniform()
            if y_value < prob[i]: #lo scattering non è avvenuto, osservo fotone
                n_obs.append(fotoni[i])

        return n_obs
    
    def R_scattering_graph(self,fotoni,n_obs, angle, k_norm):

        """
        Metodo con cui si rappresenta graficamente la distribuzione dei fotoni osservati a seguito 
        scattering di Rayleigh e la si confronta con valori attesi
        
        """
        
        plt.figure(figsize=(10,8))
        n_bins=int(mt.sqrt(len(n_obs)))
        n,bins,p=plt.hist(n_obs,bins=n_bins, color="limegreen", label="Valori simulati")
     
        #Valori attesi

        bincenters=(bins[:-1]+bins[1:])/2 
        binwidth=(bins[-1:]-bins[0])/n_bins
        height=f.D_norm(bincenters,self.T,k_norm)
        prob_0=height*binwidth 
        n_0_exp=len(fotoni)*prob_0 #numero di fotoni iniziali attesi
        
        n_exp=n_0_exp*f.prob_obs(bincenters,angle)
        
        mask=np.nonzero(n)

        chi2 =  np.sum( (n_exp[mask] - n[mask])**2 /n[mask] ) 
        ndof = len(n[mask])-1
        chi2_red=chi2/ndof
    
        plt.plot(bincenters,n_exp, color="blue", label="Valori attesi", marker=".", linewidth=1)
        plt.text(bincenters[1],max(n)*1/4, r'$\chi^2$  : {:.2f}'.format(chi2_red), fontsize=14, color='blue')
        plt.title("Stella: {:}, angolo: {:}°. \n Distribuzione simulata e attesa dei fotoni osservati considerando scattering di Ryleigh.".format(self.name,angle),
                   color="mediumblue",fontsize=13)
        plt.ylabel(r'$ \text{Numero di fotoni osservati per unità  di superficie e tempo} \  \left( \frac{1}{m^2 \cdot t }\right)  $',fontsize=9.5)
        plt.xlabel(r'$ \text{Lunghezza d\'onda (nm)} $')
        plt.legend(frameon=False)
        plt.show()

    
    def compare_distribution(self,N):
        """
        Metodo che confronta nella zona del visibile la distribuzione dei fotoni emessi dalla stella
        con quelle che considerano lo scattering di Rayleigh quando la stella è allo Zenit e all'orizzonte
        """

        fotoni_emessi=self.no_absorption(N,380,790)[0]
        fotoni_zenit=self.R_scattering(fotoni_emessi,0)
        fotoni_orizzonte=self.R_scattering(fotoni_emessi,90)

        fig, axs = plt.subplots(2, 1, gridspec_kw={'height_ratios': [4 ,0.5]},sharex=True,figsize=(10,8))
        fig.subplots_adjust(hspace=0)
        n_bins=100
        
        axs[0].set_title("Confronto diverse distribuzioni simulate di fotoni osservati - Stella: {:}".format(self.name), fontsize=12)
        axs[0].hist(fotoni_emessi,bins=n_bins, color="blueviolet", label="Senza assorbimento")
        axs[0].hist(fotoni_zenit,bins=n_bins, color="cyan", label="Con scattering e Sole allo Zenit ")
        axs[0].hist(fotoni_orizzonte,bins=n_bins, color="greenyellow",label="Con scattering e Sole all'orizzonte")
        axs[0].set_ylabel(r'$ \text{Numero di fotoni per  unità  di superficie e tempo} \  \left( \frac{1}{m^2 \cdot t }\right)  $',fontsize=9.5)
        axs[0].legend(fontsize=8, frameon=False)
        norm=plt.Normalize(380,790)
        cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["purple","blue","limegreen","yellow","orange","red"])
        plt.colorbar(plt.cm.ScalarMappable( cmap=cmap, norm=norm),cax=axs[1], orientation='horizontal',label="something")
        axs[1].set_xlabel(r'$ \text{Lunghezza d\'onda (nm)} $',fontsize=10) 
        plt.show()

    
    
    def flusso_integrato(self,N,angle):

        """
        Metodo che calcola il numero di fotoni totali osservato in funzione dell'angolo della stella 
        dallo Zenit

        Input:
        - N numero di fotoni che si vogliono utilizzare per la simulazione
        - angle: array di angoli da studiare

        Output:
        -Array degli integrali, ossia numero fotoni al variare dell'angolo

    

        """
        integrali=[]
        #Emissione fotoni
        fotoni=self.no_absorption(N,380,790)[0]

        #Simulo scattering a vari angoli
        for i in range (len(angle)):

            #Simulo all'angolo i-esimo
            n_obs=self.R_scattering(fotoni,angle[i])
            n_obs_tot=np.sum(n_obs)
            integrali.append(n_obs_tot)

        plt.figure(figsize=(10,8))
        plt.plot(angle,integrali, marker=".", linewidth=0, color="mediumseagreen")
        plt.title(r'Andamento del numero di fotoni totali osservati in funzione angolo rispetto allo zenit. Stella: {:}'.format(self.name), fontsize=11)
        plt.ylabel(r'$ \text{Numero di fotoni totali per  unità  di superficie e tempo} \  \left( \frac{1}{m^2 \cdot t }\right)  $',fontsize=9.5)
        plt.xlabel("Angolo (°)")
        plt.show()

        return integrali

"""
Gestione delle azioni con argparse
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

    #Creo alcune stelle
    sole=star("Sole",5.75*10**3)
    betelgeuse=star("Betelgeuse",3*10**3)
    bellatrix=star("Bellatrix",22*10**3)
    alfa_crucis=star("Alfa Crucis",28*10**3)
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

            #Confronto distribuzioni
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

            #Confronto distribuzioni
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

            #Confronto distribuzioni
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

            #Confronto distribuzioni
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