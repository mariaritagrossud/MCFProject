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
        Metodo che simula la distribuzione dei fotoni solari osservati in assenza di scattering di Rayleigh e
        plotta il risultato confrontandolo con quello atteso
        
        input:
        N = numero di fotoni da usare per generare la distribuzione
        estremo_sx = più piccola lunghezza d'onda che si vuole generare
        estremo_dx = più grande lunghezza d'onda che si vuole generare

        output:
        - Array con i fotoni generati dalla distribuzione
        - costate di normalizzazione da usare nel metodo successivo
        """
        num=int((estremo_dx-estremo_sx)/0.1) #faccio in modo che i punti siano equispaziati di 0.1 nm
        x=np.linspace(estremo_sx,estremo_dx,num)
        k_norm=1/integrate.simpson(f.D(x,self.T),x) #costante con cui normalizzare distribuzione
        x_max=optimize.minimize(f.D_norm_opposite,x0=(estremo_sx+estremo_dx)/2,args=(self.T,k_norm), tol=1e-9)
        max_value_1=f.D_norm(x_max.x[0],self.T,k_norm)    
        
        #Simulazione
        fotoni=[]
        for i in range(N):
            lam=np.random.uniform(low=estremo_sx,high=estremo_dx,size=None)
            y_value=np.random.uniform(low=0, high=max_value_1, size=None)
        
            if (y_value <= f.D_norm(lam,self.T,k_norm)):
                fotoni.append(lam)
        
        return fotoni, k_norm


    def no_abs_graph(self,fotoni,k_norm):
        """

        Metodo con cui si rappresenta graficamente la distribuzione e la si confronta con valori attesi

        """
        #Plot 
        plt.figure(figsize=(10,8))

        n_bins=150
        n,bins,p=plt.hist(fotoni,bins=n_bins, color="cyan",alpha=0.8, label="Valori simulati")
        
        bincenters=(bins[:-1]+bins[1:])/2
        binwidth=(bins[len(bins)-1]-bins[0])/n_bins
        height=f.D_norm(bincenters,self.T,k_norm)
        prob=height*binwidth
        n_exp=len(fotoni)*prob
       
        mask=np.nonzero(n)
        chi2 =  np.sum( (n_exp[mask] - n[mask])**2 /n[mask] ) 
        ndof = len(n[mask])-1
        chi2_red=chi2/ndof

        plt.text(bincenters[1],30, r'$\chi^2$  : {:.2f}'.format(chi2_red), fontsize=14, color='blue')
        plt.plot(bincenters,n_exp, color="blue", label="Valori attesi", marker=".", linewidth=1)
        plt.title("Distribuzione fotoni senza assorbimento simulata e attesa - Stella: {:}".format(self.name))
        plt.ylabel(r'$ \text{Numero di fotoni osservati per unità  di superficie e tempo} \  \left( \frac{1}{m^2 \cdot t }\right)  $'
                    ,fontsize=9.5)
        plt.xlabel(r'$ \text{Lunghezza d\'onda (nm)} $')
        plt.legend(frameon=False)
        plt.show()
    
    def R_scattering(self, N,estremo_sx, estremo_dx, angle) :
        """
        Metodo che simula la distribuzione dei fotoni solari osservati con scattering di Rayleigh
        
        input:
        N = numero di fotoni da usare per generare la distribuzione
        estremo_sx = più piccola lunghezza d'onda che si vuole generare
        estremo_dx = più grande lunghezza d'onda che si vuole generare
        angle = angolo della stella rispetto allo Zenit

        output:
        Array con i fotoni generati dalla distribuzione
        """
        num=int((estremo_dx-estremo_sx)/0.1)
        x=np.linspace(estremo_sx,estremo_dx,num)
        k_norm=1/integrate.simpson(f.n_obs(x,self.T,angle),x)
        x_max=optimize.minimize(f.n_obs_norm_opposite,x0=(estremo_sx+estremo_dx)/2,args=(self.T, angle,k_norm), tol=1e-9)
        max_value=f.n_obs_norm(x_max.x[0],self.T,angle,k_norm)    

        fotoni=[]
        for i in range(N):
            lam=np.random.uniform(low=estremo_sx,high=estremo_dx,size=None)
            y_value=np.random.uniform(low=0, high=max_value, size=None)
            if (y_value <= f.n_obs_norm(lam,self.T,angle,k_norm)):
                fotoni.append(lam)

        return fotoni, k_norm
    
    def R_scattering_graph(self,fotoni,angle, k_norm):
         #Plot
                
        plt.figure(figsize=(10,8))
        n_bins=150
        n,bins,p=plt.hist(fotoni,bins=n_bins, color="cyan",alpha=0.9, label="Valori simulati")
        bincenters=(bins[:-1]+bins[1:])/2
        binwidth=(bins[len(bins)-1]-bins[0])/n_bins
        height=f.n_obs_norm(bincenters,self.T,angle,k_norm)
        prob=height*binwidth
        n_exp=len(fotoni)*prob

        mask=np.nonzero(n)
        chi2 =  np.sum( (n_exp[mask] - n[mask])**2 /n[mask] ) 
        ndof = len(n[mask])-1
        chi2_red=chi2/ndof

        plt.plot(bincenters,n_exp, color="blue", label="Valori attesi", marker=".", linewidth=1)
        plt.text(bincenters[1],30, r'$\chi^2$  : {:.2f}'.format(chi2_red), fontsize=14, color='blue')
        plt.title("Distribuzione fotoni con scattering Rayleigh simulata e attesa con angolo {:}° - Stella: {:}".format(angle,self.name))
        plt.ylabel(r'$ \text{Numero di fotoni osservati per unità  di superficie e tempo} \  \left( \frac{1}{m^2 \cdot t }\right)  $',fontsize=9.5)
        plt.xlabel(r'$ \text{Lunghezza d\'onda (nm)} $')
        plt.legend(frameon=False)
        plt.show()

    
    def compare_distribution(self,N):
        """
        Metodo che confronta nella zona del visibile la distribuzione senza assorbimento con quelle che considerano 
        lo scattering di Rayleigh con stella allo Zenit e all'orizzonte
        """

        fotoni_solari=self.no_absorption(N,380,790)[0]
        fotoni_zenit=self.R_scattering(N,380,790,0)[0]
        fotoni_orizzonte=self.R_scattering(N,380,790,90)[0]

        fig, axs = plt.subplots(2, 1, gridspec_kw={'height_ratios': [4 ,0.5]},sharex=True,figsize=(10,8))
        fig.subplots_adjust(hspace=0)
        n_bins=224
        axs[0].set_title("Confronto diverse distribuzioni simulate di fotoni osservati - Stella: {:}".format(self.name), fontsize=12)
        axs[0].hist(fotoni_solari,bins=n_bins, color="blueviolet", label="Senza assorbimento")
        axs[0].hist(fotoni_zenit,bins=n_bins, color="cyan",alpha=0.8, label="Con scattering e Sole allo Zenit ")
        axs[0].hist(fotoni_orizzonte,bins=n_bins, color="greenyellow",label="Con scattering e Sole all'orizzonte")
        axs[0].set_ylabel(r'$ \text{Numero di fotoni per  unità  di superficie e tempo} \  \left( \frac{1}{m^2 \cdot t }\right)  $',fontsize=9.5)
        axs[0].legend(fontsize=8, frameon=False)
        norm=plt.Normalize(380,790)
        cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["purple","blue","limegreen","yellow","orange","red"])
        plt.colorbar(plt.cm.ScalarMappable( cmap=cmap, norm=norm),cax=axs[1], orientation='horizontal',label="something")
        axs[1].set_xlabel(r'$ \text{Lunghezza d\'onda (nm)} $',fontsize=10) 
        plt.show()

        return fotoni_solari, fotoni_zenit, fotoni_orizzonte
    
    def flusso_integrato(self,N,angle):

        """
        Metodo che calcola il numero di fotoni totali osservato in funzione dell'angolo della stella 
        dallo Zenit

        """
        integrali=[]
        #Simulo distribuzione fotoni con quell'angolo e calcolo integrale

        for i in range (len(angle)):

            #Simulo distribuzione fotoni con quell'angolo
            fotoni=self.R_scattering(N,380,790,angle[i])[0]
            n_bins=150
            n,bins=np.histogram(fotoni,bins=n_bins)
            bincenters=(bins[:-1]+bins[1:])/2
            integrale=integrate.simpson(n,bincenters)
            integrali.append(integrale)

        plt.figure(figsize=(9,7))
        plt.plot(angle,integrali, marker=".", linewidth=0, color="mediumseagreen")
        plt.title(r'Andamento del flusso integrato in funzione angolo rispetto allo zenit. Stella: {:}'.format(self.name), fontsize=11)
        plt.ylabel(r'$ \text{Numero di fotoni totali per  unità  di superficie e tempo} \  \left( \frac{1}{m^2 \cdot t }\right)  $',fontsize=9.5)
        plt.xlabel("Angolo (°)")
        plt.show()

        return integrali


"""
Simulazione eventi per le Stelle in analisi
____________________________________________

"""
sole=star("Sole",5.75*10**3)
betelgeuse=star("Betelgeuse",3*10**3)
bellatrix=star("Bellatrix",22*10**3)
alfa_crucis=star("Alfa Crucis",28*10**3)
stelle=[sole,betelgeuse,bellatrix,alfa_crucis]
#indago solo nel visibile
estremo_sx=380 #nm
estremo_dx=790 #nm
N=50000
for item in stelle:
    print("--------------------------------------------------------------------------------------------")
    print("Simulazione distribuzione senza assorbimento -  {:} ".format(item.name))
    sim1=item.no_absorption(N,estremo_sx,estremo_dx)
    item.no_abs_graph(sim1[0],sim1[1])
    print("--------------------------------------------------------------------------------------------")
    print("Simulazione distribuzione con scattering Rayleigh allo zenit -  {:} ".format(item.name))
    sim2=item.R_scattering(N,estremo_sx,estremo_dx,0)
    item.R_scattering_graph(sim2[0],0,sim2[1])
    print("--------------------------------------------------------------------------------------------")
    print("Simulazione distribuzione con scattering Rayleigh all'orizzonte -  {:} ".format(item.name))
    sim3=item.R_scattering(N,estremo_sx,estremo_dx,90)
    item.R_scattering_graph(sim3[0],90,sim3[1])
    print("--------------------------------------------------------------------------------------------")
    print("Confronto simulazione diverse distribuzioni -  {:} ".format(item.name))
    item.compare_distribution(N)
    print("--------------------------------------------------------------------------------------------")
    print("Studio andamento flusso integrato al variare angolo  -  {:} ".format(item.name))
    angle=np.linspace(0,90,200)
    sim4=item.flusso_integrato(1000,angle)
    
