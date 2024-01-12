import numpy as np
import math as mt
import matplotlib.pyplot  as plt
from scipy import integrate
from scipy import optimize
import matplotlib.colors
import sys

#Importo modulo con funzioni
sys.path.append(" ")
import func as f


"""
Creo classe con la quale simulare e analizzare la radiazione fotoni da parte di varie stelle
____________________________________________________________________________________________
"""


class star:

    def __init__(self, name, T):
        self.name=name
        self.T=T
    

    def no_absorption(self, N,estremo_sx, estremo_dx):
        """
        Metodo che simula la distribuzione dei fotoni emessi da un corpo nero e che noi vedremmo 
        se non ci fosse lo scattering di Rayleigh.
        
        Input:
        N = numero di fotoni da usare per generare la distribuzione
        estremo_sx = più piccola lunghezza d'onda che si vuole generare
        estremo_dx = più grande lunghezza d'onda che si vuole generare

        Output:
        - Array con i fotoni generati dalla distribuzione
        - costate di normalizzazione 
        """
        num = int((estremo_dx-estremo_sx)/0.1) #faccio in modo che i punti siano equispaziati di circa 0.1 nm
        x = np.linspace(estremo_sx,estremo_dx,num)
        k_norm = 1/integrate.simpson(f.D(x,self.T),x) #costante con cui normalizzare distribuzione
        x_max = optimize.minimize(f.D_norm_opposite,x0=(estremo_sx+estremo_dx)/2,args=(self.T,k_norm), tol=1e-9)
        max_value_1 = f.D_norm(x_max.x[0],self.T,k_norm)    
        
        #Simulazione con il metodo dell' Hit or Miss
        fotoni = []
        lam = np.random.uniform(low=estremo_sx,high=estremo_dx,size=N)
        y_value = np.random.uniform(low=0, high=max_value_1, size=N)
        for i in range(N):

             if (y_value[i] <= f.D_norm(lam[i],self.T,k_norm)):
                fotoni.append(lam[i])
        
        return fotoni, k_norm


    def no_abs_graph(self,fotoni,k_norm):
        """
        Metodo con cui si rappresenta graficamente la distribuzione del corpo nero e la si confronta 
        con valori attesi
        """

        plt.figure(figsize=(10,8))

        n_bins = 100
        n,bins,p = plt.hist(fotoni,bins=n_bins, color="limegreen", label="Valori simulati")
        
        #Valori attesi
        bincenters = (bins[:-1]+bins[1:])/2
        binwidth =(bins[-1:]-bins[0])/n_bins
        height = f.D_norm(bincenters,self.T,k_norm)
        prob = height*binwidth
        n_exp = len(fotoni)*prob
       
        #Chi quadro
        #Prendo solo bin con un numero di eventi maggiore o uguale a 4
        mask = n >=4
        chi2 =  np.sum( (n_exp[mask] - n[mask])**2 /n[mask] ) 
        ndof = len(n[mask])-1
        chi2_red = chi2/ndof

        plt.text(bincenters[1],max(n)*1/4, r'$\chi_r^2$  : {:.2f}'.format(chi2_red), fontsize=14, color='blue')
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

        Output:
        Array con i fotoni osservati, non diffusi

        """

        prob = f.prob_obs(fotoni,angle)
        n_obs = []
        for i in range(len(fotoni)):
            y_value = np.random.uniform()
            if y_value < prob[i]: #lo scattering non è avvenuto, osservo fotone
                n_obs.append(fotoni[i])

        return n_obs
    
    def R_scattering_graph(self,fotoni,n_obs, angle, k_norm):

        """
        Metodo con cui si rappresenta graficamente la distribuzione dei fotoni osservati a seguito 
        scattering di Rayleigh e la si confronta con valori attesi

        """
        
        plt.figure(figsize=(10,8))
        n_bins =int(mt.sqrt(len(n_obs)))

        n,bins,p = plt.hist(n_obs,bins=n_bins, color="limegreen", label="Valori simulati")
     
        #Valori attesi
        bincenters = (bins[:-1]+bins[1:])/2 
        binwidth = (bins[-1:]-bins[0])/n_bins
        height = f.D_norm(bincenters,self.T,k_norm)
        prob_0 = height*binwidth 
        n_0_exp = len(fotoni)*prob_0 #numero di fotoni iniziali attesi
        
        n_exp = n_0_exp*f.prob_obs(bincenters,angle)
        
        #Prendo solo bin con un numero di eventi maggiore o uguale a 4
        mask = n >=4
        chi2 =  np.sum( (n_exp[mask] - n[mask])**2 /n[mask] ) 
        ndof = len(n[mask])-1
        chi2_red=chi2/ndof
    
        plt.plot(bincenters,n_exp, color="blue", label="Valori attesi", marker=".", linewidth=1)
        plt.text(bincenters[1],max(n)*1/4, r'$\chi_r^2$  : {:.2f}'.format(chi2_red), fontsize=14, color='blue')
        plt.title("Stella: {:}, angolo: {:}°. \n Distribuzione simulata e attesa dei fotoni osservati considerando scattering di Ryleigh.".format(self.name,angle),
                   color="mediumblue",fontsize=13)
        plt.ylabel(r'$ \text{Numero di fotoni osservati per unità  di superficie e tempo} \  \left( \frac{1}{m^2 \cdot t }\right)  $',fontsize=9.5)
        plt.xlabel(r'$ \text{Lunghezza d\'onda (nm)} $')
        plt.legend(frameon=False)
        plt.show()

    
    def compare_distribution(self,N):
        """
        Metodo che confronta nella zona del visibile la distribuzione dei fotoni emessi dalla stella
        con quelle dei fotoni osservati che considerano lo scattering di Rayleigh quando 
        la stella è allo Zenit e all'orizzonte
        """

        fotoni_emessi = self.no_absorption(N,380,790)[0]
        fotoni_zenit = self.R_scattering(fotoni_emessi,0)
        fotoni_orizzonte = self.R_scattering(fotoni_emessi,90)

        fig, axs = plt.subplots(2, 1, gridspec_kw={'height_ratios': [4 ,0.5]},sharex=True,figsize=(10,8))
        fig.subplots_adjust(hspace=0)
        n_bins = 100
        axs[0].set_title("Stella: {:}. \n Confronto diverse distribuzioni simulate di fotoni osservati ".format(self.name), fontsize=12)
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
        -Array degli integrali, ossia numero totale fotoni al variare dell'angolo

        """

        integrali=[]
        #Emissione fotoni considerando solo quelli nel visibile
        fotoni=self.no_absorption(N,380,790)[0]

        #Simulo scattering a vari angoli
        for i in range (len(angle)):

            #Simulo per angolo i-esimo
            n_obs=self.R_scattering(fotoni,angle[i])
            n_obs_tot=len(n_obs) 
            integrali.append(n_obs_tot)

        plt.figure(figsize=(10,8))
        plt.plot(angle,integrali, marker=".", linewidth=0, color="mediumseagreen")
        plt.title("Stella: {:}. \nAndamento del numero di fotoni totali osservati in funzione angolo rispetto allo Zenit.".format(self.name), fontsize=11)
        plt.ylabel(r'$ \text{Numero di fotoni totali per  unità  di superficie e tempo} \  \left( \frac{1}{m^2 \cdot t }\right)  $',fontsize=9.5)
        plt.xlabel("Angolo (°)")
        plt.show()

        return integrali

